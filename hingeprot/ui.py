from __future__ import annotations

import os, re, glob, datetime, base64, uuid, shutil, subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import requests
import ipywidgets as W
from IPython.display import display, clear_output, HTML, FileLink

LAST_UI_STATE: dict | None = None
LAST_INPUTS: dict | None = None

def get_last_inputs() -> dict | None:
    return LAST_INPUTS

# ----------------------------- shell helpers -----------------------------
def _sh(cmd: str, cwd: str | None = None, timeout: int | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "-lc", cmd], cwd=cwd, capture_output=True, text=True, timeout=timeout)

def _ldconfig_has_libg2c() -> bool:
    return _sh(r'ldconfig -p | grep -E "libg2c\.so\.0" >/dev/null 2>&1').returncode == 0

def _ensure_libg2c() -> None:
    if _ldconfig_has_libg2c():
        return
    os.makedirs("/content", exist_ok=True)
    os.chdir("/content")
    deb1 = "gcc-3.4-base_3.4.6-6ubuntu3_amd64.deb"
    deb2 = "libg2c0_3.4.6-6ubuntu3_amd64.deb"
    url1 = f"https://old-releases.ubuntu.com/ubuntu/pool/universe/g/gcc-3.4/{deb1}"
    url2 = f"https://old-releases.ubuntu.com/ubuntu/pool/universe/g/gcc-3.4/{deb2}"
    r = _sh(f"wget -q {url1} -O {deb1}")
    if r.returncode != 0:
        raise RuntimeError(f"wget failed for {url1}\\n{r.stderr}")
    r = _sh(f"wget -q {url2} -O {deb2}")
    if r.returncode != 0:
        raise RuntimeError(f"wget failed for {url2}\\n{r.stderr}")
    _sh(f"dpkg -i {deb1} {deb2} || true")
    r = _sh("apt-get -y -qq -f install")
    if r.returncode != 0:
        raise RuntimeError(f"apt-get -f install failed:\\n{r.stderr}")
    _sh("ldconfig")
    if not _ldconfig_has_libg2c():
        raise RuntimeError("libg2c.so.0 still not found after installation.")

def _ensure_repo(fresh: bool = False) -> str:
    root = "/content/hingeprot_fortran"
    hp = os.path.join(root, "hingeprot")
    url = "https://github.com/enesemretas/hingeprot_fortran.git"
    here = os.path.abspath(__file__) if "__file__" in globals() else os.path.abspath(os.getcwd())
    running_inside = here.startswith(os.path.abspath(root) + os.sep)
    if fresh and not running_inside:
        shutil.rmtree(root, ignore_errors=True)
    if not os.path.isdir(hp):
        os.makedirs("/content", exist_ok=True)
        os.chdir("/content")
        r = _sh(f"git clone -q {url}")
        if r.returncode != 0:
            raise RuntimeError(f"git clone failed:\\n{r.stderr}")
    if not os.path.isdir(hp):
        raise RuntimeError("Repo clone incomplete: missing /content/hingeprot_fortran/hingeprot")
    return hp

def _write_runHingeProt_pl(hingeprot_dir: str, gnm_cut: float, anm_cut: float) -> str:
    gnm_i = int(round(float(gnm_cut)))
    anm_i = int(round(float(anm_cut)))
    pl_path = os.path.join(hingeprot_dir, "runHingeProt.pl")
    content = f"""#!/usr/bin/perl -w
use strict;
use File::Copy;
use FindBin;
my $home = "$FindBin::Bin";
if ($#ARGV != 1) {{ print "runHingeProt.pl <PDB_file> <chain ids>\\n"; exit; }}
my $pdb = $ARGV[0];
my $pdbCode = $ARGV[0];
my $chains = $ARGV[1];
my $dirname = "$pdbCode.$chains";
mkdir $dirname or print "cannot create $dirname\\n";
chdir $dirname or die "cannot change to $dirname\\n";
if (!-e "../$pdb") {{ die "cannot find file $pdb\\n"; }}
copy("../$pdb","$pdb");
`$home/getChain.Linux $chains $pdb > pdb`;
`$home/part1 {gnm_i} {anm_i}; $home/useblz; $home/part2`;
rename("gnm1anmvector", "$pdbCode.1vector");
rename("gnm2anmvector", "$pdbCode.2vector");
rename("hinges", "$pdbCode.hinge");
rename("newcoordinat.mds", "$pdbCode.new");
`$home/processHinges $pdbCode.new $pdbCode.hinge $pdbCode.1vector $pdbCode.2vector 15 14.0`;
rename("$pdbCode.new.moved1.pdb", "$pdbCode.mode1.pdb");
rename("$pdbCode.new.moved2.pdb", "$pdbCode.mode2.pdb");
"""
    with open(pl_path, "w", encoding="utf-8") as f:
        f.write(content)
    _sh("chmod +x ./runHingeProt.pl || true", cwd=hingeprot_dir)
    _sh("find . -maxdepth 1 -type f -exec chmod +x {} \\; || true", cwd=hingeprot_dir)
    return pl_path

def _safe_html(text: str) -> str:
    return (text or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def _find_hinges_file(out_dir: str, pdb_filename: str) -> str | None:
    candidates = [
        os.path.join(out_dir, f"{pdb_filename}.new.hinges"),
        os.path.join(out_dir, f"{pdb_filename}.new.hinge"),
        os.path.join(out_dir, f"{pdb_filename}.hinges"),
        os.path.join(out_dir, f"{pdb_filename}.hinge"),
        os.path.join(out_dir, "hinges"),
        os.path.join(out_dir, "hinge"),
    ]
    for p in candidates:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None

# ------------------------- parse *.new.hinges report -------------------------
_NEW_MODE_RE = re.compile(r"^---->\\s*Slowest\\s+mode\\s+(\\d+)\\s*:", re.IGNORECASE)
_NEW_NPART_RE = re.compile(r"#\\s*of\\s*rigid\\s*parts\\s*:\\s*(\\d+)", re.IGNORECASE)
_NEW_HINGE_RE = re.compile(r"^Hinge\\s+residues\\s*:\\s*(.*)$", re.IGNORECASE)

def _find_new_hinges_report(out_dir: str, pdb_filename: str) -> str | None:
    base = os.path.basename(pdb_filename)
    base_noext = os.path.splitext(base)[0]
    candidates = [
        os.path.join(out_dir, f"{base}.new.hinges"),
        os.path.join(out_dir, f"{base_noext}.new.hinges"),
        os.path.join(out_dir, f"{base_noext}.pdb.new.hinges"),
        os.path.join(out_dir, f"{base}.new.hinge"),
        os.path.join(out_dir, f"{base_noext}.new.hinge"),
    ]
    for p in candidates:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    found = sorted(Path(out_dir).glob("*.new.hinges"), key=lambda x: x.stat().st_size, reverse=True)
    return str(found[0]) if found else None

def parse_new_hinges_report(report_path: Path) -> dict[int, dict[str, object]]:
    report: dict[int, dict[str, object]] = {}
    mode: int | None = None
    with report_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            m = _NEW_MODE_RE.match(s)
            if m:
                mode = int(m.group(1))
                report.setdefault(mode, {"n_parts": None, "parts": [], "hinge_tokens": [], "hinge_token_set": set()})
                continue
            if mode is None:
                continue
            m = _NEW_NPART_RE.search(s)
            if m:
                report[mode]["n_parts"] = int(m.group(1))
                continue
            if s.startswith("Part "):
                mnum = re.match(r"^Part\\s+(\\d+)\\s+", s)
                if not mnum:
                    continue
                pno = int(mnum.group(1))
                rest = s[mnum.end():]
                idx = rest.lower().find("direction")
                residues_str = rest[:idx].rstrip() if idx != -1 else rest.rstrip()
                report[mode]["parts"].append((pno, residues_str))
                continue
            mh = _NEW_HINGE_RE.match(s)
            if mh:
                tail = (mh.group(1) or "").strip().replace(",", " ")
                tokens = re.findall(r"-?\\d+[A-Za-z]?", tail)
                report[mode]["hinge_tokens"] = tokens
                report[mode]["hinge_token_set"] = set(tokens)
                continue
    for m in report:
        parts = list(report[m].get("parts", []) or [])
        parts.sort(key=lambda x: x[0])
        report[m]["parts"] = parts
    return report

# ------------------------- Short Flexible Fragments -------------------------
def _leading_int_str(resid: str) -> str | None:
    m = re.match(r"^(-?\\d+)", (resid or "").strip())
    return m.group(1) if m else None

def _resnum_int(resid: str) -> int | None:
    s = _leading_int_str(resid)
    if s is None:
        return None
    try:
        return int(s)
    except Exception:
        return None

def _fmt_resid_with_chain(resid: object, ch: str) -> str:
    s = str(resid).strip()
    if not s:
        return s
    if s[-1].isalpha() and s[-1].upper() == ch.upper():
        return s
    return f"{s}{ch}"

def _parse_pdb_like_chain_resid(line: str) -> Optional[Tuple[str, str]]:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    if len(line) >= 27:
        chain = line[21].strip()
        resnum = line[22:26].strip()
        icode = line[26].strip()
        resid = (resnum + icode).strip()
        if chain and resid:
            return chain, resid
    parts = line.split()
    if len(parts) >= 6:
        chain = parts[4].strip()
        resid = parts[5].strip()
        if chain and resid:
            return chain, resid
    return None

def read_residue_order_from_pdb(pdb_path: Path) -> Dict[str, List[str]]:
    residues_by_chain: Dict[str, List[str]] = {}
    last: Optional[Tuple[str, str]] = None
    with pdb_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parsed = _parse_pdb_like_chain_resid(line)
            if not parsed:
                continue
            chain, resid = parsed
            key = (chain, resid)
            if key == last:
                continue
            residues_by_chain.setdefault(chain, []).append(resid)
            last = key
    return residues_by_chain

def parse_hinge_file(hinge_path: Path) -> Dict[int, Dict[str, List[Tuple[int, str]]]]:
    modes: Dict[int, Dict[str, List[Tuple[int, str]]]] = {}
    mode: Optional[int] = None
    def _mode_from_header(s: str) -> Optional[int]:
        s_low = s.lower()
        if "1st" in s_low: return 1
        if "2nd" in s_low: return 2
        m = re.search(r"(\\d+)", s_low)
        return int(m.group(1)) if m else None
    with hinge_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("---->"):
                mode = _mode_from_header(line)
                continue
            if mode is None:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                seq_idx = int(float(parts[0]))
            except Exception:
                continue
            resid_token_raw = parts[-2].strip().rstrip(",;")
            chain = parts[-1].strip()[:1]
            if not resid_token_raw or not chain:
                continue
            modes.setdefault(mode, {}).setdefault(chain, []).append((seq_idx, resid_token_raw))
    for m in modes:
        for ch in modes[m]:
            modes[m][ch].sort(key=lambda x: x[0])
    return modes

def compute_short_flexible_fragments(
    residues_by_chain: dict[str, list[str]],
    hinge_modes: dict[int, dict[str, list[tuple[int, str]]]],
    report: dict[int, dict[str, object]],
    min_len: int = 14,
) -> dict[int, list[str]]:
    out: dict[int, list[str]] = {}
    for mode, ch_map in hinge_modes.items():
        hinge_set = set(report.get(mode, {}).get("hinge_token_set", set()) or set())
        mode_fragments: list[str] = []
        for ch, entries in (ch_map or {}).items():
            if not entries:
                continue
            start_label = residues_by_chain.get(ch, [None])[0] if residues_by_chain.get(ch) else None
            start_num = _resnum_int(start_label) if start_label else None
            hinges = []
            for seq_idx, resid in entries:
                resid = str(resid).strip()
                tok = _fmt_resid_with_chain(resid, ch)
                rn = _resnum_int(resid)
                if tok in hinge_set or resid in hinge_set or (rn is not None and str(rn) in hinge_set):
                    continue
                hinges.append({"seq": int(seq_idx), "resid": resid, "resnum": _resnum_int(resid)})
            hinges.sort(key=lambda x: x["seq"])
            removed_int: list[tuple[int, int]] = []
            removed_str: list[tuple[str, str]] = []
            def gap(i: int, j: int) -> int:
                ri, rj = hinges[i].get("resnum"), hinges[j].get("resnum")
                if isinstance(ri, int) and isinstance(rj, int) and ri != rj:
                    return abs(rj - ri)
                return abs(int(hinges[j]["seq"]) - int(hinges[i]["seq"]))
            def head_len() -> int:
                if not hinges: return 0
                r1 = hinges[0].get("resnum")
                if isinstance(start_num, int) and isinstance(r1, int):
                    return abs(r1 - start_num) + 1
                return abs(int(hinges[0]["seq"]) - 1) + 1
            while True:
                if not hinges: break
                if len(hinges) >= 2 and gap(0, 1) < min_len:
                    h1, h2 = hinges[0], hinges[1]
                    r1, r2 = h1.get("resnum"), h2.get("resnum")
                    if isinstance(r1, int) and isinstance(r2, int):
                        a = r1 + 1; b = r2
                        if a <= b: removed_int.append((a, b))
                    else:
                        removed_str.append((str(h1.get("resid")), str(h2.get("resid"))))
                    del hinges[1]; del hinges[0]
                    continue
                if head_len() < min_len:
                    h1 = hinges[0]; r1 = h1.get("resnum")
                    if isinstance(start_num, int) and isinstance(r1, int):
                        a = start_num; b = r1
                        if a <= b: removed_int.append((a, b))
                    else:
                        if start_label: removed_str.append((str(start_label), str(h1.get("resid"))))
                    del hinges[0]
                    continue
                if len(hinges) < 2: break
                removed_any = False
                for i in range(len(hinges) - 1):
                    if gap(i, i + 1) < min_len:
                        h1, h2 = hinges[i], hinges[i + 1]
                        r1, r2 = h1.get("resnum"), h2.get("resnum")
                        if isinstance(r1, int) and isinstance(r2, int):
                            a = r1 + 1; b = r2
                            if a <= b: removed_int.append((a, b))
                        else:
                            removed_str.append((str(h1.get("resid")), str(h2.get("resid"))))
                        del hinges[i + 1]; del hinges[i]
                        removed_any = True
                        break
                if removed_any: continue
                break
            removed_int.sort()
            merged: list[tuple[int, int]] = []
            for a, b in removed_int:
                if not merged: merged.append((a, b))
                else:
                    pa, pb = merged[-1]
                    if a <= pb + 1: merged[-1] = (pa, max(pb, b))
                    else: merged.append((a, b))
            for a, b in merged:
                mode_fragments.append(f"{_fmt_resid_with_chain(a, ch)}-{_fmt_resid_with_chain(b, ch)}")
            for a, b in removed_str:
                mode_fragments.append(f"{_fmt_resid_with_chain(a, ch)}-{_fmt_resid_with_chain(b, ch)}")
        seen=set(); cleaned=[]
        for x in mode_fragments:
            if x not in seen:
                seen.add(x); cleaned.append(x)
        out[mode]=cleaned
    return out

# ------------------------- py3Dmol helpers -------------------------
def _ensure_py3dmol():
    try:
        import py3Dmol  # type: ignore
        return py3Dmol
    except Exception:
        r = _sh("python3 -m pip -q install py3Dmol")
        if r.returncode != 0:
            raise RuntimeError("py3Dmol install failed. Try: !pip -q install py3Dmol")
        import py3Dmol  # type: ignore
        return py3Dmol

def _html_with_unique_divid(raw_html: str) -> str:
    m = re.search(r'id="([^"]+)"', raw_html)
    if not m:
        return raw_html
    old = m.group(1)
    new = f"hp3d_{uuid.uuid4().hex}"
    return raw_html.replace(old, new)

def _bfactor_minmax(pdb_text: str) -> tuple[float, float]:
    mn = None; mx = None
    for ln in pdb_text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")): continue
        if len(ln) < 66: continue
        try: b=float(ln[60:66])
        except Exception: continue
        mn = b if mn is None else min(mn, b)
        mx = b if mx is None else max(mx, b)
    if mn is None or mx is None: return (0.0, 100.0)
    if abs(mx - mn) < 1e-9: return (mn - 1.0, mx + 1.0)
    return (mn, mx)

def _standardize_pdb_for_3dmol(pdb_text: str) -> str:
    out=[]
    for ln in (pdb_text or "").splitlines():
        s = ln.rstrip("\\n")
        if s.startswith(("MODEL","ENDMDL")):
            out.append(s); continue
        if s.startswith("END"): continue
        if not s.startswith(("ATOM","HETATM")):
            out.append(s); continue
        parts = s.split()
        if len(parts) < 11:
            out.append(s); continue
        rec = parts[0]
        try:
            serial=int(float(parts[1])); atom=parts[2]; resn=parts[3]
            chain=(parts[4] or "A")[0]; resi=int(float(parts[5]))
            x=float(parts[6]); y=float(parts[7]); z=float(parts[8])
            occ=float(parts[9]); b=float(parts[10])
        except Exception:
            out.append(s); continue
        atom_name=(atom or "CA").strip()
        atom_field=atom_name.rjust(4)[:4]
        resn_field=(resn or "UNK")[:3].rjust(3)
        altloc=" "; icode=" "
        elem=(atom_name[0] if atom_name else "C").upper()[:2].rjust(2)
        out.append(
            f"{rec:<6}{serial:>5} {atom_field}{altloc}{resn_field} {chain:1}{resi:>4}{icode}   "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{b:>6.2f}          {elem}"
        )
    return "\\n".join(out).rstrip() + "\\n"

def _to_multimodel_pdb_string(pdb_text: str) -> tuple[str, int]:
    txt = _standardize_pdb_for_3dmol(pdb_text)
    lines = txt.splitlines()
    if not any(ln.startswith("MODEL") for ln in lines):
        return "MODEL        1\\n" + txt + "ENDMDL\\nEND\\n", 1
    models=[]; cur=None
    for ln in lines:
        if ln.startswith("MODEL"):
            cur=[]; continue
        if ln.startswith("ENDMDL"):
            if cur is not None and any(x.startswith(("ATOM","HETATM")) for x in cur):
                models.append(cur)
            cur=None; continue
        if ln.startswith("END"): continue
        if cur is not None:
            cur.append(ln)
    if not models:
        return "MODEL        1\\n" + txt + "ENDMDL\\nEND\\n", 1
    out=[]
    for i, mlines in enumerate(models, start=1):
        out.append(f"MODEL        {i}")
        out.extend(mlines)
        out.append("ENDMDL")
    out.append("END")
    return "\\n".join(out) + "\\n", len(models)

def _render_mode_into_output(out_widget: W.Output, mode_pdb_path: str) -> None:
    py3Dmol = _ensure_py3dmol()
    with out_widget:
        clear_output(wait=True)
        pdb_text = Path(mode_pdb_path).read_text(encoding="utf-8", errors="ignore")
        if not pdb_text.strip():
            display(HTML("<div style='color:#dc2626;font-weight:800;'>Empty mode file.</div>"))
            return
        mm_text, nmodels = _to_multimodel_pdb_string(pdb_text)
        bmin, bmax = _bfactor_minmax(mm_text)
        v = py3Dmol.view(width=560, height=280)
        v.setBackgroundColor("white")
        if nmodels >= 2:
            v.addModelsAsFrames(mm_text, "pdb")
        else:
            v.addModel(mm_text, "pdb")
        cs = {"prop":"b", "gradient":"roygb", "min":float(bmin), "max":float(bmax)}
        v.setStyle({}, {"trace":{"colorscheme":cs}, "cartoon":{"colorscheme":cs}})
        v.zoomTo()
        raw = _html_with_unique_divid(v._make_html())
        display(HTML(raw))

# ------------------------- report widget -------------------------
def rigidparts_report_widget_from_report(pdb_label: str, report: dict[int, dict[str, object]], short_frags_by_mode: dict[int, list[str]], out_dir: str, download_fn) -> W.VBox:
    def _locate_mode_pdb(out_dir_: str, pdb_label_: str, mode_: int) -> str | None:
        cands=[os.path.join(out_dir_, f"{pdb_label_}.mode{mode_}.pdb"), os.path.join(out_dir_, f"{os.path.splitext(pdb_label_)[0]}.mode{mode_}.pdb")]
        for p in cands:
            if os.path.exists(p) and os.path.getsize(p)>0: return p
        hits=sorted(glob.glob(os.path.join(out_dir_, f"*mode{mode_}.pdb")), key=os.path.getsize, reverse=True)
        for p in hits:
            if os.path.getsize(p)>0: return p
        return None

    blocks=[]
    mode_paths={}
    for m in (1,2):
        p=_locate_mode_pdb(out_dir, pdb_label, m)
        if p: mode_paths[m]=p

    if mode_paths:
        title=W.HTML("<div style='text-align:center;font-family:Arial;font-weight:900;color:#111827;margin:6px 0 2px 0;'>Mode Viewer</div>")
        opts=[(f"Mode {m}", str(m)) for m in sorted(mode_paths)]
        mode_sel=W.ToggleButtons(options=opts, value=opts[0][1], layout=W.Layout(width="320px"), style={"button_width":"150px"})
        mode_out=W.Output(layout=W.Layout(width="580px", border="1px solid #e5e7eb", border_radius="12px", padding="6px"))
        def _do_render(*_):
            _render_mode_into_output(mode_out, mode_paths[int(mode_sel.value)])
        mode_sel.observe(lambda ch: _do_render(), names="value")
        _do_render()
        blocks.append(W.VBox([title, W.HBox([mode_sel], layout=W.Layout(justify_content="center")), mode_out], layout=W.Layout(width="100%", gap="6px", align_items="center")))

    def _css_cell() -> str:
        return "padding:6px 8px;border-bottom:1px solid #e5e7eb;font-size:13px;"

    for mode in sorted(report):
        blocks.append(W.HTML(f"<div style='width:100%; text-align:center; color:#dc2626; font-weight:900;'>----&gt; Slowest Mode {mode}: {pdb_label}</div>"))
        if mode in (1,2):
            fpath=_locate_mode_pdb(out_dir, pdb_label, mode)
            if fpath:
                btn=W.Button(description=f"Download mode{mode}.pdb", icon="download", layout=W.Layout(width="220px"))
                btn.on_click(lambda _b, p=fpath: download_fn(p))
                blocks.append(W.HBox([btn], layout=W.Layout(justify_content="center")))
        n_parts=report[mode].get("n_parts")
        parts=list(report[mode].get("parts", []) or [])
        rows=[]
        for pno, residues_str in parts:
            rows.append(f"<tr><td style='{_css_cell()}'>{pno}</td><td style='{_css_cell()}'>{residues_str}</td></tr>")
        hinge_tokens=list(report[mode].get("hinge_tokens", []) or [])
        hinge_line=" ".join(hinge_tokens) if hinge_tokens else "-"
        frags=short_frags_by_mode.get(mode, []) or []
        if frags:
            items=[f"<div style='margin:2px 0;'>{k}. {frag}</div>" for k, frag in enumerate(frags, start=1)]
            short_html="<div style='margin-top:10px;color:#dc2626;font-weight:900;'>Short Flexible Fragments:</div>"+ "".join(items)
        else:
            short_html="<div style='margin-top:10px;color:#dc2626;font-weight:900;'>Short Flexible Fragments:</div><div style='margin:2px 0;'>-</div>"
        body="<div style='font-family:Arial;'>"
        if isinstance(n_parts, int):
            body+=f"<div style='margin:4px 0 8px 0;font-weight:800;'># of rigid parts: {n_parts}</div>"
        body+=(
            f"<table style='width:100%;border-collapse:collapse;'>"
            f"<thead><tr><th style='text-align:left;{_css_cell()}border-bottom:2px solid #e5e7eb;'>Rigid Part No</th>"
            f"<th style='text-align:left;{_css_cell()}border-bottom:2px solid #e5e7eb;'>Residues</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table>"
            f"<div style='margin-top:6px;color:#1d4ed8;font-weight:900; text-align:center;'>Hinge residues: {hinge_line}</div>"
            f"{short_html}</div>"
        )
        blocks.append(W.HTML(body))
    return W.VBox(blocks, layout=W.Layout(width="100%", gap="6px"))

# ----------------------------- UI -----------------------------
def launch(runs_root: str = "/content/hingeprot_runs"):
    from google.colab import output
    output.enable_custom_widget_manager()
    os.makedirs(runs_root, exist_ok=True)

    def _fetch_pdb_text(pdb_code: str) -> str:
        code = pdb_code.strip().upper()
        if not re.fullmatch(r"[0-9A-Z]{4}", code):
            raise ValueError("PDB code must be 4 characters (e.g., 4CLN).")
        url = f"https://files.rcsb.org/download/{code}.pdb"
        r = requests.get(url, timeout=30)
        if r.status_code != 200 or len(r.text) < 200:
            raise RuntimeError(f"Failed to fetch PDB {code} (HTTP {r.status_code}).")
        return r.text

    def _detect_chains_from_text(pdb_text: str) -> list[str]:
        chains=set()
        for line in pdb_text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                ch=line[21].strip()
                if ch: chains.add(ch)
        return sorted(chains)

    def _list_or_custom_float(label: str, options, default_value: float, minv: float, maxv: float):
        opts = sorted(set(float(x) for x in options) | {float(default_value)})
        lbl=W.Label(label, layout=W.Layout(width="120px"))
        toggle=W.ToggleButtons(options=[("List","list"),("Custom","custom")], value="list", layout=W.Layout(width="180px"), style={"button_width":"80px"})
        dropdown=W.Dropdown(options=opts, value=float(default_value), layout=W.Layout(width="240px"))
        fbox=W.BoundedFloatText(value=float(default_value), min=minv, max=maxv, step=0.1, layout=W.Layout(width="240px"))
        value_box=W.Box([dropdown], layout=W.Layout(align_items="center"))
        def _on_toggle(ch):
            value_box.children = [dropdown] if ch["new"]=="list" else [fbox]
        toggle.observe(_on_toggle, names="value")
        def get_value() -> float:
            return float(dropdown.value) if toggle.value=="list" else float(fbox.value)
        row=W.HBox([lbl, toggle, value_box], layout=W.Layout(align_items="center", gap="12px"))
        return row, get_value

    _CHAIN_PALETTE=["red","blue","green","orange","purple","cyan","magenta","yellow","teal","brown","pink","lime","navy","gold"]
    def _assign_chain_colors(chains: list[str]) -> dict[str,str]:
        return {ch:_CHAIN_PALETTE[i % len(_CHAIN_PALETTE)] for i, ch in enumerate(chains)}

    css = W.HTML("<style>.hp-card{border:1px solid #e5e7eb;border-radius:14px;padding:14px 16px;background:#fff;box-sizing:border-box}.hp-pre{white-space:pre-wrap;font-family:ui-monospace,monospace;font-size:13px;line-height:1.35;background:#0b1020;color:#e5e7eb;padding:12px;border-radius:12px;border:1px solid #1f2937}</style>")
    nav = W.ToggleButtons(options=[("Web Server","web"),("Help","help")], value="web", style={"button_width":"140px"})

    input_mode=W.ToggleButtons(options=[("Enter PDB code","code"),("Upload PDB file","upload")], value="code", style={"button_width":"170px"}, layout=W.Layout(width="420px"))
    pdb_code=W.Text(value="", description="PDB code:", placeholder="e.g., 4cln", style={"description_width":"80px"}, layout=W.Layout(width="420px"))
    btn_choose_file=W.Button(description="Choose file", icon="upload", layout=W.Layout(width="180px"))
    upload_prog=W.IntProgress(value=0, min=0, max=100, layout=W.Layout(width="160px")); upload_prog.bar_style=""
    file_lbl=W.Label("No file chosen")
    btn_load=W.Button(description="Load / Detect Chains", button_style="info", icon="search", layout=W.Layout(width="260px"))

    all_chains=W.Checkbox(value=False, description="All Chains", indent=False, layout=W.Layout(width="120px"))
    chains_wrap=W.Box([], layout=W.Layout(display="flex", flex_flow="row wrap", gap="10px", border="1px solid #e5e7eb", border_radius="12px", padding="8px 10px"))

    gnm_row, get_gnm_cut=_list_or_custom_float("GNM cutoff (Å):", [7,8,9,10,11,12,13,20], 10.0, 1.0, 100.0)
    anm_row, get_anm_cut=_list_or_custom_float("ANM cutoff (Å):", [10,13,15,18,20,23,36], 18.0, 1.0, 100.0)

    progress=W.IntProgress(value=0, min=0, max=4, description="Progress:", bar_style="")
    btn_run_fortran=W.Button(description="Run HingeProt", button_style="primary", icon="play", layout=W.Layout(width="320px"))
    btn_clear=W.Button(description="Clear", button_style="warning", icon="trash", layout=W.Layout(width="180px"))

    table_box=W.VBox([], layout=W.Layout(width="100%", gap="8px"))
    status_box=W.HTML('<div class="hp-pre">Load a PDB to detect chains.</div>')
    output_info=W.HTML("", layout=W.Layout(width="100%"))

    def _set_status(t: str):
        status_box.value=f'<div class="hp-pre">{_safe_html(t)}</div>'

    def _download_file(path: str):
        path=os.path.abspath(path)
        if (not os.path.exists(path)) or os.path.getsize(path)==0:
            _set_status(f"ERROR: file not found or empty:\\n{path}")
            return
        try:
            from google.colab import files  # type: ignore
            files.download(path)
        except Exception:
            display(FileLink(path))

    # viewer
    CARD_W=620; CARD_PAD=16; OUT_PAD=6; OUT_BORDER=1
    INNER_W=CARD_W-2*CARD_PAD
    VIEW_W=INNER_W-2*(OUT_PAD+OUT_BORDER); VIEW_H=280
    viewer_out=W.Output(layout=W.Layout(width="100%", height=f"{VIEW_H}px", border="1px solid #e5e7eb", border_radius="12px", padding=f"{OUT_PAD}px", overflow="hidden"))
    viewer_card=W.VBox([W.HTML("<b>3D Viewer</b>"), viewer_out], layout=W.Layout(width=f"{CARD_W}px", gap="10px"))
    viewer_card.add_class("hp-card")

    def _viewer_placeholder(msg: str="Load a PDB to preview it here."):
        with viewer_out:
            clear_output(wait=True)
            print(msg)
    _viewer_placeholder()

    state={"pdb_text":None,"pdb_filename":None,"pdb_path":None,"run_dir":None,"pdb_tag":None,"upload_name":None,"upload_bytes":None,
           "detected_chains":[],"chain_cbs":{},"chain_colors":{},"manual_selection":(),"_syncing":False,"last_out_dir":None}
    global LAST_UI_STATE
    LAST_UI_STATE=state

    def _selected_chains() -> list[str]:
        detected=state.get("detected_chains", [])
        return [ch for ch in detected if ch in state["chain_cbs"] and state["chain_cbs"][ch].value]

    def _set_selection(sel: list[str]):
        detected=state.get("detected_chains", [])
        sel=[c for c in sel if c in detected]
        state["_syncing"]=True
        try:
            for ch, cb in state["chain_cbs"].items():
                cb.value=(ch in sel)
        finally:
            state["_syncing"]=False

    def _update_all_checkbox_from_selection():
        if state["_syncing"]: return
        detected=state.get("detected_chains", [])
        if not detected: return
        sel=_selected_chains()
        all_now=(len(sel)==len(detected)) and (len(detected)>0)
        state["_syncing"]=True
        try: all_chains.value=all_now
        finally: state["_syncing"]=False
        if not all_now:
            state["manual_selection"]=tuple(sel)

    def _refresh_viewer():
        if not state.get("pdb_text"):
            _viewer_placeholder()
            return
        py3Dmol=_ensure_py3dmol()
        pdb_text=state["pdb_text"]
        detected=state.get("detected_chains", [])
        selected=list(detected) if (detected and all_chains.value) else _selected_chains()
        chain_colors=state.get("chain_colors", {}) or {}
        v=py3Dmol.view(width=VIEW_W, height=VIEW_H)
        v.setBackgroundColor("white")
        v.addModel(pdb_text, "pdb")
        v.setStyle({}, {"cartoon":{"color":"lightgray"}})
        for ch in selected:
            v.setStyle({"chain":ch}, {"cartoon":{"color":chain_colors.get(ch,"red")}})
        v.setStyle({"hetflag": True, "not":{"resn":["HOH","WAT","DOD"]}}, {"stick":{}, "sphere":{"scale":0.25}})
        v.zoomTo()
        raw=_html_with_unique_divid(v._make_html())
        with viewer_out:
            clear_output(wait=True)
            display(HTML(raw))

    def _sync_input_visibility(*_):
        pdb_code.layout.display="" if input_mode.value=="code" else "none"
        btn_choose_file.layout.display="" if input_mode.value=="upload" else "none"
        upload_prog.layout.display="" if input_mode.value=="upload" else "none"
        file_lbl.layout.display="" if input_mode.value=="upload" else "none"
    _sync_input_visibility()
    input_mode.observe(lambda ch: _sync_input_visibility(), names="value")

    cb_name=f"hingeprot_uploader_{uuid.uuid4().hex}"
    cb_prog=f"hingeprot_uploadprog_{uuid.uuid4().hex}"

    def _js_upload_progress_callback(payload):
        try:
            pct=int(payload.get("pct", 0)); pct=max(0,min(100,pct))
            upload_prog.value=pct
            upload_prog.bar_style="info" if pct<100 else "success"
        except Exception:
            pass

    def _js_upload_callback(payload):
        try:
            name=payload.get("name","upload.pdb")
            data_b64=payload.get("data_b64","")
            if not data_b64:
                _set_status("Upload callback received empty data.")
                return
            data=base64.b64decode(data_b64.encode("utf-8"))
            state["upload_name"]=name
            state["upload_bytes"]=data
            file_lbl.value=name
            upload_prog.value=100; upload_prog.bar_style="success"
            _set_status(f"Uploaded file: {name} ({len(data)} bytes)\\nNow click 'Load / Detect Chains'.")
        except Exception as e:
            _set_status(f"Upload callback error: {e}")

    output.register_callback(cb_prog, _js_upload_progress_callback)
    output.register_callback(cb_name, _js_upload_callback)

    def on_choose_file(_):
        upload_prog.value=0; upload_prog.bar_style="info"
        js=f"""
        (async () => {{
          const input = document.createElement('input');
          input.type='file'; input.accept='.pdb,.ent'; input.style.display='none';
          document.body.appendChild(input);
          input.onchange = async () => {{
            const file = input.files && input.files[0];
            document.body.removeChild(input);
            if (!file) return;
            const reader = new FileReader();
            reader.onloadstart = async () => {{ try {{ await google.colab.kernel.invokeFunction("{cb_prog}", [{{pct:0}}], {{}}); }} catch(e){{}} }};
            reader.onprogress = async (e) => {{
              try {{
                if (e.lengthComputable) {{
                  const pct = Math.round((e.loaded/e.total)*100);
                  await google.colab.kernel.invokeFunction("{cb_prog}", [{{pct:pct}}], {{}});
                }}
              }} catch(e){{}}
            }};
            reader.onloadend = async () => {{ try {{ await google.colab.kernel.invokeFunction("{cb_prog}", [{{pct:100}}], {{}}); }} catch(e){{}} }};
            reader.onload = async () => {{
              const b64 = (reader.result || "").split(",")[1] || "";
              await google.colab.kernel.invokeFunction("{cb_name}", [{{name:file.name, data_b64:b64}}], {{}});
            }};
            reader.readAsDataURL(file);
          }};
          input.click();
        }})();
        """
        output.eval_js(js)

    btn_choose_file.on_click(on_choose_file)

    def _on_chain_cb_change(_):
        _update_all_checkbox_from_selection()
        _refresh_viewer()

    def _rebuild_chain_checkboxes(chains: list[str], default_selected: list[str]):
        state["chain_cbs"]={}
        items=[all_chains]
        for ch in chains:
            cb=W.Checkbox(value=(ch in default_selected), description=ch, indent=False, layout=W.Layout(width="48px"))
            cb.observe(_on_chain_cb_change, names="value")
            state["chain_cbs"][ch]=cb
            items.append(cb)
        chains_wrap.children=tuple(items)

    def _on_all_chains_toggle(ch):
        if state["_syncing"]:
            return
        detected=state.get("detected_chains", [])
        if not detected:
            return
        if ch["new"] is True:
            sel=_selected_chains()
            if len(sel) != len(detected):
                state["manual_selection"]=tuple(sel)
            _set_selection(detected)
        else:
            prev=list(state.get("manual_selection") or [])
            prev=[c for c in prev if c in detected]
            if not prev: prev=[detected[0]]
            _set_selection(prev)
        _update_all_checkbox_from_selection()
        _refresh_viewer()

    all_chains.observe(_on_all_chains_toggle, names="value")

    def on_load_clicked(_):
        progress.value=0; table_box.children=(); progress.bar_style="info"; output_info.value=""
        try:
            ts=datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            if input_mode.value=="upload":
                if state["upload_bytes"] is None:
                    raise ValueError("Please click 'Choose file' and upload a PDB first.")
                pdb_text=state["upload_bytes"].decode("utf-8", errors="ignore")
                upname=(state.get("upload_name") or "upload.pdb").strip()
                if not re.search(r"\\.(pdb|ent)$", upname, flags=re.I):
                    upname += ".pdb"
                pdb_filename=os.path.basename(upname)
                base=os.path.splitext(pdb_filename)[0]
                tag=re.sub(r"[^0-9A-Za-z]+","",base).upper() or "UPLOAD"
            else:
                code=pdb_code.value.strip()
                if not code:
                    raise ValueError("Please enter a PDB code (e.g., 3lzg).")
                code4=code.upper()
                pdb_text=_fetch_pdb_text(code4)
                pdb_filename=f"{code4.lower()}.pdb"
                tag=code4
            run_dir=os.path.join(runs_root, f"{tag}_run_{ts}")
            os.makedirs(run_dir, exist_ok=True)
            state["run_dir"]=run_dir; state["pdb_text"]=pdb_text; state["pdb_filename"]=pdb_filename; state["pdb_tag"]=tag
            pdb_path=os.path.join(run_dir, pdb_filename)
            Path(pdb_path).write_text(pdb_text, encoding="utf-8")
            state["pdb_path"]=pdb_path
            chs=_detect_chains_from_text(pdb_text)
            if not chs:
                raise RuntimeError("No chains detected in the PDB.")
            state["detected_chains"]=chs
            state["chain_colors"]=_assign_chain_colors(chs)
            default_sel=[chs[0]]
            state["manual_selection"]=tuple(default_sel)
            _rebuild_chain_checkboxes(chs, default_sel)
            state["_syncing"]=True
            try: all_chains.value=False
            finally: state["_syncing"]=False
            progress.value=1
            _set_status(f"Loaded PDB (ID={tag})\\nRun folder: {run_dir}\\nDetected chains: {', '.join(chs)}\\n")
            _refresh_viewer()
        except Exception as e:
            progress.bar_style="danger"
            _set_status(f"ERROR: {e}")

    def _capture_inputs() -> dict:
        if not state.get("pdb_text"):
            raise RuntimeError("Please click 'Load / Detect Chains' first.")
        detected=state.get("detected_chains", [])
        if not detected:
            raise RuntimeError("No detected chains. Load again.")
        chain_list=detected if all_chains.value else _selected_chains()
        if not chain_list:
            raise RuntimeError("Please select at least one chain (or tick All chains).")
        return {
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
            "pdb_tag": state.get("pdb_tag"),
            "pdb_filename": state.get("pdb_filename"),
            "run_dir_runsroot": state.get("run_dir"),
            "chains_str": "".join(chain_list),
            "gnm_cutoff_A": float(get_gnm_cut()),
            "anm_cutoff_A": float(get_anm_cut()),
        }

    def on_run_fortran_clicked(_):
        progress.bar_style="info"
        try:
            captured=_capture_inputs()
            global LAST_INPUTS
            LAST_INPUTS=captured
            table_box.children=(W.HTML("<div style='font-family:Arial; color:#6b7280;'>Running HingeProt…</div>"),)
            progress.value=1
            _ensure_libg2c()
            progress.value=2
            hp_dir=_ensure_repo(fresh=False)
            progress.value=3
            _write_runHingeProt_pl(hp_dir, captured["gnm_cutoff_A"], captured["anm_cutoff_A"])
            pdb_filename=captured["pdb_filename"]
            (Path(hp_dir)/pdb_filename).write_text(state["pdb_text"] or "", encoding="utf-8")
            chains_str=captured["chains_str"]
            out_dir_repo=os.path.join(hp_dir, f"{pdb_filename}.{chains_str}")
            if os.path.isdir(out_dir_repo):
                shutil.rmtree(out_dir_repo, ignore_errors=True)
            r=_sh(f"perl ./runHingeProt.pl {pdb_filename} {chains_str}", cwd=hp_dir, timeout=1200)
            if r.returncode != 0:
                raise RuntimeError(f"runHingeProt.pl failed (return code {r.returncode}).\\n{r.stderr}")
            run_dir=captured["run_dir_runsroot"]
            if not run_dir or not os.path.isdir(run_dir):
                raise RuntimeError("Run folder not found. Please 'Load / Detect Chains' again.")
            dest_out_dir=os.path.join(run_dir, os.path.basename(out_dir_repo))
            if os.path.isdir(dest_out_dir):
                shutil.rmtree(dest_out_dir, ignore_errors=True)
            if os.path.isdir(out_dir_repo):
                shutil.move(out_dir_repo, dest_out_dir)
            else:
                raise RuntimeError(f"Expected output folder not found: {out_dir_repo}")
            state["last_out_dir"]=dest_out_dir
            output_info.value=f"<div style='font-family:Arial;font-weight:800;margin:2px 0 6px 0;'>{_safe_html(pdb_filename)} for Chain(s): {_safe_html(', '.join(list(chains_str)))}</div>"
            pdb_chain_path=Path(dest_out_dir)/"pdb"
            if not pdb_chain_path.exists():
                pdb_chain_path=Path(dest_out_dir)/pdb_filename
            report_file=_find_new_hinges_report(dest_out_dir, pdb_filename)
            if not report_file:
                raise RuntimeError(f"Report file not found: expected '{pdb_filename}.new.hinges' in {dest_out_dir}")
            report=parse_new_hinges_report(Path(report_file))
            hinge_path=Path(dest_out_dir)/f"{pdb_filename}.hinge"
            if (not hinge_path.exists()) or (hinge_path.stat().st_size==0):
                alt=_find_hinges_file(dest_out_dir, pdb_filename)
                if alt: hinge_path=Path(alt)
                else: raise RuntimeError(f".hinge file not found: {hinge_path}")
            residues_by_chain=read_residue_order_from_pdb(pdb_chain_path)
            hinge_modes=parse_hinge_file(hinge_path)
            short_frags=compute_short_flexible_fragments(residues_by_chain, hinge_modes, report, min_len=14)
            table_box.children=(rigidparts_report_widget_from_report(pdb_filename, report, short_frags, dest_out_dir, _download_file),)
            progress.value=4; progress.bar_style="success"
        except Exception as e:
            progress.bar_style="danger"
            table_box.children=(W.HTML("<div style='color:#dc2626;font-weight:800;'>ERROR</div>"),)
            _set_status(f"ERROR: {e}")

    def on_clear_clicked(_):
        pdb_code.value=""; table_box.children=(); input_mode.value="code"
        state["upload_name"]=None; state["upload_bytes"]=None
        file_lbl.value="No file chosen"; upload_prog.value=0; upload_prog.bar_style=""
        state.update({"pdb_text":None,"pdb_filename":None,"pdb_path":None,"run_dir":None,"pdb_tag":None,"last_out_dir":None,
                      "detected_chains":[], "chain_cbs":{}, "chain_colors":{}, "manual_selection":(), "_syncing":False})
        all_chains.value=False; chains_wrap.children=()
        output_info.value=""
        global LAST_INPUTS
        LAST_INPUTS=None
        progress.value=0; progress.bar_style=""
        _set_status("Cleared. Load a PDB to detect chains.")
        _viewer_placeholder()

    btn_load.on_click(on_load_clicked)
    btn_run_fortran.on_click(on_run_fortran_clicked)
    btn_clear.on_click(on_clear_clicked)

    form_card=W.VBox([
        W.HBox([W.HTML("<b>Input:</b>"), input_mode]),
        pdb_code,
        W.HBox([btn_choose_file, upload_prog, file_lbl]),
        btn_load,
        W.HTML("<hr>"),
        W.HBox([W.HTML("<b>Select Chains:</b>"), chains_wrap], layout=W.Layout(align_items="center", gap="12px")),
        W.VBox([gnm_row, anm_row], layout=W.Layout(gap="8px")),
        progress,
        W.HBox([btn_run_fortran, btn_clear]),
    ], layout=W.Layout(width="620px", gap="10px"))
    form_card.add_class("hp-card")

    output_card=W.VBox([output_info, table_box, status_box], layout=W.Layout(width="100%", gap="8px"))
    output_card.add_class("hp-card")

    web_page=W.VBox([W.HBox([form_card, viewer_card], layout=W.Layout(display="flex", flex_flow="row wrap", gap="14px")), output_card], layout=W.Layout(width="100%", gap="10px"))

    help_page=W.HTML("<div style='font-family:Arial;line-height:1.4;'><h3>Help</h3><ul><li>Load a PDB (code or upload), select chain(s), then run HingeProt.</li><li>Mode viewer normalizes ATOM lines so CA-only files render reliably.</li></ul></div>")

    main_view=W.VBox([web_page], layout=W.Layout(width="100%"))
    def _switch(_):
        main_view.children=[web_page] if nav.value=="web" else [help_page]
    nav.observe(_switch, names="value"); _switch(None)

    display(css, nav, main_view)
    return None
