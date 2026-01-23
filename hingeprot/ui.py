from __future__ import annotations

import os
import re
import datetime
import base64
import uuid
import shutil
import subprocess

import requests
import ipywidgets as W
from IPython.display import display, clear_output, HTML, FileLink
from pathlib import Path
from typing import Dict, List, Tuple, Optional



LAST_UI_STATE: dict | None = None
LAST_INPUTS: dict | None = None


def get_last_inputs() -> dict | None:
    return LAST_INPUTS


# ----------------------------- shell helpers -----------------------------
def _sh(cmd: str, cwd: str | None = None, timeout: int | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _ldconfig_has_libg2c() -> bool:
    r = _sh(r'ldconfig -p | grep -E "libg2c\.so\.0" >/dev/null 2>&1')
    return r.returncode == 0


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
        raise RuntimeError(f"wget failed for {url1}\n{r.stderr}")

    r = _sh(f"wget -q {url2} -O {deb2}")
    if r.returncode != 0:
        raise RuntimeError(f"wget failed for {url2}\n{r.stderr}")

    _sh(f"dpkg -i {deb1} {deb2} || true")
    r = _sh("apt-get -y -qq -f install")
    if r.returncode != 0:
        raise RuntimeError(f"apt-get -f install failed:\n{r.stderr}")

    _sh("ldconfig")
    if not _ldconfig_has_libg2c():
        raise RuntimeError("libg2c.so.0 still not found after installation.")


def _ensure_repo(fresh: bool = False) -> str:
    root = "/content/hingeprot_fortran"
    hp = os.path.join(root, "hingeprot")
    url = "https://github.com/enesemretas/hingeprot_fortran.git"

    # SAFE: __file__ yoksa cwd kullan
    here = os.path.abspath(__file__) if "__file__" in globals() else os.path.abspath(os.getcwd())
    running_inside = here.startswith(os.path.abspath(root) + os.sep)

    if fresh and not running_inside:
        shutil.rmtree(root, ignore_errors=True)

    if not os.path.isdir(hp):
        os.makedirs("/content", exist_ok=True)
        os.chdir("/content")
        r = _sh(f"git clone -q {url}")
        if r.returncode != 0:
            raise RuntimeError(f"git clone failed:\n{r.stderr}")

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

if ($#ARGV != 1) {{
  print "runHingeProt.pl <PDB_file> <chain ids>\\n";
  exit;
}}

my $pdb = $ARGV[0];
my $pdbCode = $ARGV[0];
my $chains = $ARGV[1];

my $dirname = "$pdbCode.$chains";

mkdir $dirname or print "cannot create $dirname\\n";
chdir $dirname or die "cannot change to $dirname\\n";

if (!-e "../$pdb") {{
  die "cannot find file $pdb\\n";
}}

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


def _read_text_file(path: str, max_lines: int = 900) -> str:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    if len(lines) > max_lines:
        head = lines[: max_lines // 2]
        tail = lines[-max_lines // 2 :]
        lines = head + ["", "[... truncated ...]", ""] + tail
    return "\n".join(lines)

def _read_optional_file(path: str) -> str:
    try:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            return _read_text_file(path, max_lines=20000)
    except Exception:
        pass
    return ""


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

# ------------------------- NEW: parse .new.hinges report -------------------------

_NEW_MODE_RE = re.compile(r"^---->\s*Slowest\s+mode\s+(\d+)\s*:", re.IGNORECASE)
_NEW_NPART_RE = re.compile(r"#\s*of\s*rigid\s*parts\s*:\s*(\d+)", re.IGNORECASE)
_NEW_HINGE_RE = re.compile(r"^Hinge\s+residues\s*:\s*(.*)$", re.IGNORECASE)

def _find_new_hinges_report(out_dir: str, pdb_filename: str) -> str | None:
    base = os.path.basename(pdb_filename)
    base_noext = os.path.splitext(base)[0]

    candidates = [
        os.path.join(out_dir, f"{base}.new.hinges"),         # 3lzg.pdb.new.hinges
        os.path.join(out_dir, f"{base_noext}.new.hinges"),   # 3lzg.new.hinges
        os.path.join(out_dir, f"{base_noext}.pdb.new.hinges"),
        os.path.join(out_dir, f"{base}.new.hinge"),
        os.path.join(out_dir, f"{base_noext}.new.hinge"),
    ]
    for p in candidates:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p

    # fallback: dizinde ne varsa bul (en büyük dosyayı seç)
    try:
        found = sorted(
            Path(out_dir).glob("*.new.hinges"),
            key=lambda x: x.stat().st_size,
            reverse=True,
        )
        if found:
            return str(found[0])
    except Exception:
        pass

    return None



def parse_new_hinges_report(report_path: Path) -> dict[int, dict[str, object]]:
    """
    Returns:
      report[mode] = {
        "n_parts": int|None,
        "parts": list[tuple[int, str]]  # (part_no, residues_str) residues_str = direction öncesi
        "hinge_tokens": list[str]       # e.g. ["45A","305A","60B","87B"]
        "hinge_token_set": set[str]     # aynı tokenların set'i (hızlı lookup)
      }
    """
    report: dict[int, dict[str, object]] = {}
    mode: int | None = None

    with report_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            s = line.strip()
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
                # "Part 2 A:46-55,...  direction 0.0 ..."
                # İstenen: direction'a kadar aynen
                # 1) part no
                mnum = re.match(r"^Part\s+(\d+)\s+", s)
                if not mnum:
                    continue
                pno = int(mnum.group(1))
                rest = s[mnum.end():]  # after "Part k "
                # 2) cut at "direction"
                idx = rest.lower().find("direction")
                residues_str = rest[:idx].rstrip() if idx != -1 else rest.rstrip()
                # aynı part no tekrar gelirse overwrite yerine append; rapor genelde tekil
                report[mode]["parts"].append((pno, residues_str))
                continue

            mh = _NEW_HINGE_RE.match(s)
            if mh:
                tail = (mh.group(1) or "").strip()
                tail = tail.replace(",", " ")
                tokens = re.findall(r"-?\d+[A-Za-z]?", tail)
                report[mode]["hinge_tokens"] = tokens
                report[mode]["hinge_token_set"] = set(tokens)
                continue

    # parts sorted
    for m in report:
        parts = list(report[m].get("parts", []))
        parts.sort(key=lambda x: x[0])
        report[m]["parts"] = parts

    return report


# ------------------------- NEW: short flexible fragments -------------------------

def _leading_int_str(resid: str) -> str | None:
    t = (resid or "").strip()
    m = re.match(r"^(-?\d+)", t)
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
    # zaten zincir harfiyle bitiyorsa tekrar ekleme
    if s[-1].isalpha() and s[-1].upper() == ch.upper():
        return s
    return f"{s}{ch}"

def compute_short_flexible_fragments(
    residues_by_chain: dict[str, list[str]],
    hinge_modes: dict[int, dict[str, list[tuple[int, str]]]],  # parse_hinge_file çıktın
    report: dict[int, dict[str, object]],
    min_len: int = 14,
) -> dict[int, list[str]]:
    """
    For each mode:
      - take hinges from .hinge
      - remove those that are in report's 'Hinge residues' list
      - apply old 'restart removal' logic (min_len)
      - output fragments with chain IDs at ends: '46A-55A'
    """
    out: dict[int, list[str]] = {}

    for mode, ch_map in hinge_modes.items():
        hinge_set = set(report.get(mode, {}).get("hinge_token_set", set()) or set())

        mode_fragments: list[str] = []

        for ch, entries in (ch_map or {}).items():
            if not entries:
                continue

            # start residue id for head-pair
            start_label = None
            if ch in residues_by_chain and residues_by_chain[ch]:
                start_label = residues_by_chain[ch][0]
            start_num = _resnum_int(start_label) if start_label else None

            # filter hinges: remove those in Hinge residues
            hinges = []
            for seq_idx, resid in entries:
                resid = str(resid).strip()
                tok = _fmt_resid_with_chain(resid, ch)
                rn = _resnum_int(resid)
                if (
                    tok in hinge_set
                    or str(resid).strip() in hinge_set
                    or (rn is not None and str(rn) in hinge_set)
                ):
                    continue

                hinges.append({
                    "seq": int(seq_idx),
                    "resid": resid,
                    "resnum": _resnum_int(resid),
                })

            hinges.sort(key=lambda x: x["seq"])

            removed_ranges_int: list[tuple[int, int]] = []
            removed_ranges_str: list[tuple[str, str]] = []  # fallback if ints yok

            def gap(i: int, j: int) -> int:
                ri, rj = hinges[i].get("resnum"), hinges[j].get("resnum")
                if isinstance(ri, int) and isinstance(rj, int) and ri != rj:
                    return abs(rj - ri)
                # resnum yoksa / eşitse -> seq farkı
                return abs(int(hinges[j]["seq"]) - int(hinges[i]["seq"]))

            def head_len() -> int:
                if not hinges:
                    return 0
                r1 = hinges[0].get("resnum")
                if isinstance(start_num, int) and isinstance(r1, int):
                    return abs(r1 - start_num) + 1
                # seq index fallback: assume first residue ~1
                return abs(int(hinges[0]["seq"]) - 1) + 1

            # --- restart loop (same spirit as your old code) ---
            while True:
                if not hinges:
                    break

                # Rule: first two hinge too close
                if len(hinges) >= 2 and gap(0, 1) < min_len:
                    h1, h2 = hinges[0], hinges[1]
                    r1, r2 = h1.get("resnum"), h2.get("resnum")

                    if isinstance(r1, int) and isinstance(r2, int):
                        a = r1 + 1
                        b = r2
                        if a <= b:
                            removed_ranges_int.append((a, b))
                    else:
                        # fallback: can't do +1 safely
                        a = str(h1.get("resid"))
                        b = str(h2.get("resid"))
                        removed_ranges_str.append((a, b))

                    del hinges[1]
                    del hinges[0]
                    continue

                # Rule: head segment too short
                if head_len() < min_len:
                    h1 = hinges[0]
                    r1 = h1.get("resnum")
                    if isinstance(start_num, int) and isinstance(r1, int):
                        a = start_num
                        b = r1
                        if a <= b:
                            removed_ranges_int.append((a, b))
                    else:
                        # requirement: if first hinge missing but second exists,
                        # pair with first residue id (label) (we do it here)
                        if start_label:
                            removed_ranges_str.append((str(start_label), str(h1.get("resid"))))
                    del hinges[0]
                    continue

                # scan other pairs
                if len(hinges) < 2:
                    break

                removed_any = False
                for i in range(len(hinges) - 1):
                    if gap(i, i + 1) < min_len:
                        h1, h2 = hinges[i], hinges[i + 1]
                        r1, r2 = h1.get("resnum"), h2.get("resnum")

                        if isinstance(r1, int) and isinstance(r2, int):
                            a = r1 + 1
                            b = r2
                            if a <= b:
                                removed_ranges_int.append((a, b))
                        else:
                            a = str(h1.get("resid"))
                            b = str(h2.get("resid"))
                            removed_ranges_str.append((a, b))

                        del hinges[i + 1]
                        del hinges[i]
                        removed_any = True
                        break

                if removed_any:
                    continue
                break

            # merge int ranges
            removed_ranges_int.sort()
            merged: list[tuple[int, int]] = []
            for a, b in removed_ranges_int:
                if not merged:
                    merged.append((a, b))
                else:
                    pa, pb = merged[-1]
                    if a <= pb + 1:
                        merged[-1] = (pa, max(pb, b))
                    else:
                        merged.append((a, b))

            # render fragments with chain at ends
            for a, b in merged:
                mode_fragments.append(f"{_fmt_resid_with_chain(a, ch)}-{_fmt_resid_with_chain(b, ch)}")
            for a, b in removed_ranges_str:
                # NOTE: a already may include insertion; still append chain at ends
                mode_fragments.append(f"{_fmt_resid_with_chain(a, ch)}-{_fmt_resid_with_chain(b, ch)}")

        # stabilize / dedup (keep order-ish)
        seen = set()
        cleaned = []
        for x in mode_fragments:
            if x not in seen:
                seen.add(x)
                cleaned.append(x)

        out[mode] = cleaned

    return out


# ------------------------- NEW: HTML renderer (mode-based, no chain blocks) -------------------------

def rigidparts_report_html_from_report(
    pdb_label: str,
    report: dict[int, dict[str, object]],
    short_frags_by_mode: dict[int, list[str]],
    out_dir: str,  # NEW
) -> str:
    def _css_cell() -> str:
        return "padding:6px 8px;border-bottom:1px solid #e5e7eb;font-size:13px;"

    blocks: list[str] = []
    for mode in sorted(report.keys()):

        # --- NEW: mode1/mode2 download link (if file exists) ---
        btn_html = ""
        if mode in (1, 2):
            fname = f"{pdb_label}.mode{mode}.pdb"  # e.g. 3lzg.pdb.mode1.pdb
            fpath = os.path.join(out_dir, fname)
            if os.path.exists(fpath) and os.path.getsize(fpath) > 0:
                # Colab serves /content/... as /files/content/...
                href = f"/files{fpath}"
                btn_html = f"<a class='hp-dlbtn' href='{href}' download='{fname}'>Download mode{mode}.pdb</a>"

        blocks.append(
            f"<div class='hp-modebar'>"
            f"  <div class='hp-modetitle' style='color:#dc2626;font-weight:900;'>----&gt; Slowest Mode {mode}: {pdb_label}</div>"
            f"  {btn_html}"
            f"</div>"
        )

        n_parts = report[mode].get("n_parts", None)
        if isinstance(n_parts, int):
            blocks.append(f"<div style='margin:4px 0 8px 0;font-weight:800;'># of rigid parts: {n_parts}</div>")

        parts = list(report[mode].get("parts", []) or [])
        rows = []
        for pno, residues_str in parts:
            rows.append(
                f"<tr>"
                f"<td style='{_css_cell()}'>{pno}</td>"
                f"<td style='{_css_cell()}'>{residues_str}</td>"
                f"</tr>"
            )

        hinge_tokens = list(report[mode].get("hinge_tokens", []) or [])
        hinge_line = " ".join(hinge_tokens) if hinge_tokens else "-"

        frags = short_frags_by_mode.get(mode, []) or []
        if frags:
            items = [f"<div style='margin:2px 0;'>{k}. {frag}</div>" for k, frag in enumerate(frags, start=1)]
            short_html = (
                "<div style='margin-top:10px;color:#dc2626;font-weight:900;'>Short Flexible Fragments:</div>"
                + "".join(items)
            )
        else:
            short_html = (
                "<div style='margin-top:10px;color:#dc2626;font-weight:900;'>Short Flexible Fragments:</div>"
                "<div style='margin:2px 0;'>-</div>"
            )

        blocks.append(
            f"<div style='margin:8px 0 14px 0;'>"
            f"<table style='width:100%;border-collapse:collapse;'>"
            f"<thead><tr>"
            f"<th style='text-align:left;{_css_cell()}border-bottom:2px solid #e5e7eb;'>Rigid Part No</th>"
            f"<th style='text-align:left;{_css_cell()}border-bottom:2px solid #e5e7eb;'>Residues</th>"
            f"</tr></thead>"
            f"<tbody>{''.join(rows) if rows else ''}</tbody>"
            f"</table>"
            f"<div style='margin-top:6px;color:#1d4ed8;font-weight:900; text-align:center;'>Hinge residues: {hinge_line}</div>"
            f"{short_html}"
            f"</div>"
        )

    return "<div style='font-family:Arial, Helvetica, sans-serif;'>" + "".join(blocks) + "</div>"

def rigidparts_report_widget_from_report(
    pdb_label: str,
    report: dict[int, dict[str, object]],
    short_frags_by_mode: dict[int, list[str]],
    out_dir: str,
    download_fn,
) -> W.VBox:

    def _css_cell() -> str:
        return "padding:6px 8px;border-bottom:1px solid #e5e7eb;font-size:13px;"

    blocks: list[W.Widget] = []

    for mode in sorted(report.keys()):


        # --- header text ---
        header_text = W.HTML(
            f"<div style='width:100%; text-align:center; color:#dc2626; font-weight:900;'>"
            f"----&gt; Slowest Mode {mode}: {pdb_label}"
            f"</div>",
            layout=W.Layout(width="100%"),
        )


        # --- inline download button (text'in hemen sonrasında) ---
        btn = None
        if mode in (1, 2):
            fname = f"{pdb_label}.mode{mode}.pdb"  # 4cln.pdb.mode1.pdb
            fpath = os.path.join(out_dir, fname)
            if os.path.exists(fpath) and os.path.getsize(fpath) > 0:
                btn = W.Button(
                    description=f"Download mode{mode}.pdb",
                    button_style="",  # istersen "info"
                    icon="download",
                    layout=W.Layout(width="200px"),
                )
                btn.on_click(lambda _b, p=fpath: download_fn(p))

        blocks.append(header_text)

        if btn is not None:
            blocks.append(
                W.HBox(
                    [btn],
                    layout=W.Layout(width="100%", justify_content="center", align_items="center")
                )
            )


        # --- body html (table + hinge + short fragments) ---
        n_parts = report[mode].get("n_parts", None)
        parts = list(report[mode].get("parts", []) or [])
        rows = []
        for pno, residues_str in parts:
            rows.append(
                f"<tr>"
                f"<td style='{_css_cell()}'>{pno}</td>"
                f"<td style='{_css_cell()}'>{residues_str}</td>"
                f"</tr>"
            )

        hinge_tokens = list(report[mode].get("hinge_tokens", []) or [])
        hinge_line = " ".join(hinge_tokens) if hinge_tokens else "-"

        frags = short_frags_by_mode.get(mode, []) or []
        if frags:
            items = [f"<div style='margin:2px 0;'>{k}. {frag}</div>" for k, frag in enumerate(frags, start=1)]
            short_html = (
                "<div style='margin-top:10px;color:#dc2626;font-weight:900;'>Short Flexible Fragments:</div>"
                + "".join(items)
            )
        else:
            short_html = (
                "<div style='margin-top:10px;color:#dc2626;font-weight:900;'>Short Flexible Fragments:</div>"
                "<div style='margin:2px 0;'>-</div>"
            )

        body_html = "<div style='font-family:Arial, Helvetica, sans-serif;'>"
        if isinstance(n_parts, int):
            body_html += f"<div style='margin:4px 0 8px 0;font-weight:800;'># of rigid parts: {n_parts}</div>"

        body_html += (
            f"<div style='margin:8px 0 14px 0;'>"
            f"<table style='width:100%;border-collapse:collapse;'>"
            f"<thead><tr>"
            f"<th style='text-align:left;{_css_cell()}border-bottom:2px solid #e5e7eb;'>Rigid Part No</th>"
            f"<th style='text-align:left;{_css_cell()}border-bottom:2px solid #e5e7eb;'>Residues</th>"
            f"</tr></thead>"
            f"<tbody>{''.join(rows) if rows else ''}</tbody>"
            f"</table>"
            f"<div style='margin-top:6px;color:#1d4ed8;font-weight:900; text-align:center;'>Hinge residues: {hinge_line}</div>"
            f"{short_html}"
            f"</div>"
            f"</div>"
        )

        blocks.append(W.HTML(body_html))

    return W.VBox(blocks, layout=W.Layout(width="100%", gap="6px"))





# ------------------------- rigid parts report helpers -------------------------

def _parse_pdb_like_chain_resid(line: str) -> Optional[Tuple[str, str]]:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None

    # 1) fixed-column (daha güvenli)
    if len(line) >= 27:
        chain = line[21].strip()
        resnum = line[22:26].strip()
        icode = line[26].strip()
        resid = (resnum + icode).strip()
        if chain and resid:
            return chain, resid

    # 2) fallback split
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
            parsed = _parse_pdb_like_chain_resid(line)  # seninki zaten ATOM/HETATM parse ediyor
            if not parsed:
                continue
            chain, resid = parsed
            key = (chain, resid)
            if key == last:
                continue
            residues_by_chain.setdefault(chain, []).append(resid)
            last = key

    return residues_by_chain

_STRIP_TRAIL_RE = re.compile(r"^(-?\d+)([A-Za-z]+)$")

def _strip_trailing_letters(token: str) -> str:
    """
    '45A' -> '45' (chain zaten ayrı kolonda geldiği için)
    '123' -> '123'
    '45A,' -> '45'  (sondaki virgül/noktalı virgül temizlenir)
    """
    t = (token or "").strip().rstrip(",;")
    m = _STRIP_TRAIL_RE.match(t)
    return m.group(1) if m else t


def parse_hinge_file(hinge_path: Path) -> Dict[int, Dict[str, List[Tuple[int, str]]]]:
    """
    Parses a .hinge file with headers like:
      ----> crosscorrelation : 1st slowest mode
    and data lines like:
      <seq_idx> <resid> <chain>
    Returns: modes[mode][chain] = [(seq_idx, resid), ...]
    """
    modes: Dict[int, Dict[str, List[Tuple[int, str]]]] = {}
    mode: Optional[int] = None

    def _mode_from_header(s: str) -> Optional[int]:
        s_low = s.lower()
        if "1st" in s_low:
            return 1
        if "2nd" in s_low:
            return 2
        m = re.search(r"(\d+)", s_low)
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

            # resid'i aynen tut (insertion code vs bozulmasın)
            resid_token = resid_token_raw

            modes.setdefault(mode, {}).setdefault(chain, []).append((seq_idx, resid_token))


    for m in modes:
        for ch in modes[m]:
            modes[m][ch].sort(key=lambda x: x[0])

    return modes


# ----------------------------- UI -----------------------------
def launch(runs_root: str = "/content/hingeprot_runs"):
    # ---- ENV DETECTION: Colab vs Localhost ----
    IS_COLAB = False
    try:
        from google.colab import output as colab_output  # type: ignore
        colab_output.enable_custom_widget_manager()
        IS_COLAB = True
    except Exception:
        colab_output = None  # local Jupyter

    os.makedirs(runs_root, exist_ok=True)

    # ---------- helpers ----------
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
        chains = set()
        for line in pdb_text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                ch = line[21].strip()
                if ch:
                    chains.add(ch)
        return sorted(chains)

    def _safe_html(text: str) -> str:
        return (text or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
 
    def _linkify_citation_line(line: str) -> str:
        """
        Citation satırı içinde geçen URL / DOI / PubMed PMID'leri tıklanabilir yapar.
        DİKKAT: Bu fonksiyon HTML döndürür (li içine direkt konacak).
        """
        s = (line or "").strip()
        if not s:
            return ""

        # Önce HTML escape (güvenlik)
        s = _safe_html(s)

        # 1) URL'leri link yap (http/https)
        # escaped string içinde URL zaten düz kalır, sadece < > & yok.
        url_pat = re.compile(r"(https?://[^\s<]+)")
        s = url_pat.sub(r'<a href="\1" target="_blank" rel="noopener noreferrer">\1</a>', s)

        # 2) DOI yakala (DOI: 10.... veya düz 10....)
        # Basit ama pratik: 10.<digits>/<stuff>
        doi_pat = re.compile(r"\b(DOI:\s*)?(10\.\d{4,9}/[^\s;,]+)\b", re.IGNORECASE)
        def _doi_repl(m):
            prefix = m.group(1) or ""
            doi = m.group(2)
            # DOI linki
            url = f"https://doi.org/{doi}"
            # prefix (DOI:) görünür kalsın
            return f'{prefix}<a href="{url}" target="_blank" rel="noopener noreferrer">{doi}</a>'
        s = doi_pat.sub(_doi_repl, s)

        # 3) PubMed / PMID yakala (PMID: 12345678 ya da sadece PMID 12345678)
        pmid_pat = re.compile(r"\bPMID[:\s]+(\d{5,10})\b", re.IGNORECASE)
        s = pmid_pat.sub(
            r'PMID: <a href="https://pubmed.ncbi.nlm.nih.gov/\1/" target="_blank" rel="noopener noreferrer">\1</a>',
            s,
        )

        return s

    def _list_or_custom_float(label: str, options, default_value: float, minv: float, maxv: float):
        opts = [float(x) for x in options]
        default_value = float(default_value)
        if default_value not in opts:
            opts = sorted(set(opts + [default_value]))
        else:
            opts = sorted(set(opts))

        lbl = W.Label(label, layout=W.Layout(width="120px"))
        toggle = W.ToggleButtons(
            options=[("List", "list"), ("Custom", "custom")],
            value="list",
            layout=W.Layout(width="180px"),
            style={"button_width": "80px"},
        )
        dropdown = W.Dropdown(options=opts, value=default_value, layout=W.Layout(width="240px"))
        fbox = W.BoundedFloatText(value=default_value, min=minv, max=maxv, step=0.1, layout=W.Layout(width="240px"))
        value_box = W.Box([dropdown], layout=W.Layout(align_items="center"))

        def _on_toggle(ch):
            value_box.children = [dropdown] if ch["new"] == "list" else [fbox]

        toggle.observe(_on_toggle, names="value")

        def get_value() -> float:
            return float(dropdown.value) if toggle.value == "list" else float(fbox.value)

        row = W.HBox([lbl, toggle, value_box], layout=W.Layout(align_items="center", gap="12px"))
        return row, get_value

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
        """
        py3Dmol'un ürettiği HTML içindeki ilk div id'sini bulup
        tüm referanslarda benzersiz bir id ile değiştirir.
        """
        m = re.search(r'id="([^"]+)"', raw_html)
        if not m:
            return raw_html
        old = m.group(1)
        new = f"hp3d_{uuid.uuid4().hex}"
        return raw_html.replace(old, new)

    def _wrap_html_doc(snippet: str) -> str:
        """py3Dmol HTML snippet'ini tam bir HTML dokümanına sarar."""
        return (
            "<!doctype html><html><head><meta charset='utf-8'></head>"
            "<body style='margin:0;overflow:hidden;'>"
            f"{snippet}"
            "</body></html>"
        )

    def _iframe_from_html_doc(doc_html: str, w: int, h: int) -> str:
        """
        /files veya data: kullanmadan iframe içinde HTML'i çalıştırır.
        (srcdoc + base64 => localhost problemlerini bitirir)
        """
        b64 = base64.b64encode(doc_html.encode("utf-8")).decode("ascii")

        # srcdoc attribute'u dışarıda tek tırnakla kullanıldığı için içeride tek tırnak kullanmıyoruz.
        srcdoc = (
            "<!doctype html><html><head><meta charset=\"utf-8\"></head>"
            "<body style=\"margin:0;overflow:hidden;\">"
            "<script>"
            f"const html=atob(\"{b64}\");"
            "document.open();document.write(html);document.close();"
            "</script>"
            "</body></html>"
        )

        return (
            f"<iframe srcdoc='{srcdoc}' "
            f"sandbox='allow-scripts allow-same-origin' "
            f"style='width:{w}px;height:{h}px;border:1px solid #e5e7eb;"
            f"border-radius:12px;overflow:hidden;'></iframe>"
        )

    
    def _wrap_html_doc(snippet: str) -> str:
        return (
            "<!doctype html><html><head><meta charset='utf-8'></head>"
            "<body style='margin:0;overflow:hidden;'>"
            f"{snippet}"
            "</body></html>"
        )
    
    def _write_html_and_get_iframe(raw_html_doc: str, save_dir: str, w: int, h: int) -> str:
        """
        Colab: /files + absolute path works.
        Local: /files + *relative-to-cwd* path works (absolute path usually DOES NOT).
        If file is outside cwd, copy it under cwd/.hp_iframes.
        """
        save_dir = os.path.abspath(save_dir or os.getcwd())
        tmp_dir = os.path.join(save_dir, ".hp_iframes")
        os.makedirs(tmp_dir, exist_ok=True)
    
        html_path = os.path.join(tmp_dir, f"hp3d_{uuid.uuid4().hex}.html")
        Path(html_path).write_text(raw_html_doc, encoding="utf-8")
    
        if IS_COLAB:
            src = f"/files{html_path}"  # absolute ok in Colab
        else:
            # local: need relpath under current notebook root (cwd)
            cwd = os.path.abspath(os.getcwd())
            rel = os.path.relpath(html_path, start=cwd)
            if rel.startswith(".."):
                # copy under cwd to ensure /files can serve it
                safe_dir = os.path.join(cwd, ".hp_iframes")
                os.makedirs(safe_dir, exist_ok=True)
                safe_path = os.path.join(safe_dir, os.path.basename(html_path))
                shutil.copy2(html_path, safe_path)
                rel = os.path.relpath(safe_path, start=cwd)
    
            rel = rel.replace(os.sep, "/")
            src = f"/files/{rel}"
    
        return (
            f"<iframe src='{src}' "
            f"sandbox='allow-scripts allow-same-origin' "
            f"style='width:{w}px;height:{h}px;border:1px solid #e5e7eb;"
            f"border-radius:12px;overflow:hidden;'></iframe>"
        )

    

    def _read_pdb_for_frames(path: str) -> str:
        """MODEL/ENDMDL/END satırlarını temizleyip frame olarak kullanılabilir pdb text döndürür."""
        txt = Path(path).read_text(encoding="utf-8", errors="ignore")
        out_lines = []
        for ln in txt.splitlines():
            if ln.startswith(("MODEL", "ENDMDL", "END")):
                continue
            out_lines.append(ln)
        return "\n".join(out_lines).rstrip() + "\n"

    def _bfactor_minmax(pdb_text: str) -> tuple[float, float]:
        mn = None
        mx = None
        for ln in pdb_text.splitlines():
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            if len(ln) < 66:
                continue
            try:
                b = float(ln[60:66])
            except Exception:
                continue
            mn = b if mn is None else min(mn, b)
            mx = b if mx is None else max(mx, b)
        if mn is None or mx is None:
            return (0.0, 100.0)
        if abs(mx - mn) < 1e-9:
            return (mn - 1.0, mx + 1.0)
        return (mn, mx)

    def _build_multimodel_pdb_string(pdb_text: str) -> tuple[str, int]:
        """
        Her durumda geçerli bir multi-model PDB string'i üretir:
        MODEL        1
        ...atoms...
        ENDMDL
        MODEL        2
        ...atoms...
        ENDMDL
        Döndürür: (multi_model_string, n_models)
        """
        lines = (pdb_text or "").splitlines()
    
        # END, MODEL, ENDMDL gibi satırları kontrol ederek model bloklarını çıkar
        has_model = any(ln.startswith("MODEL") for ln in lines)
        models: list[list[str]] = []
    
        if has_model:
            cur: list[str] = []
            in_model = False
    
            for ln in lines:
                if ln.startswith("MODEL"):
                    # yeni model başlıyor -> önceki bloğu kapat
                    if cur:
                        models.append(cur)
                        cur = []
                    in_model = True
                    continue
    
                if ln.startswith("ENDMDL"):
                    # model bitti
                    models.append(cur)
                    cur = []
                    in_model = False
                    continue
    
                if ln.startswith("END"):
                    continue  # END satırlarını at
    
                # model dışında kalan boş/garip satırları çok büyütmeyelim
                if (not in_model) and (not ln.startswith(("ATOM", "HETATM", "TER", "CONECT"))):
                    continue
    
                cur.append(ln)
    
            if cur:
                models.append(cur)
    
        else:
            # tek modelmiş gibi al; END satırlarını at
            models = [[ln for ln in lines if not ln.startswith("END")]]
    
        # boş modelleri temizle (hiç ATOM/HETATM yoksa at)
        cleaned: list[list[str]] = []
        for m in models:
            if any(ln.startswith(("ATOM", "HETATM")) for ln in m):
                cleaned.append(m)
        if cleaned:
            models = cleaned
    
        # multi-model string üret
        multi = ""
        for i, m in enumerate(models, start=1):
            multi += f"MODEL        {i}\n"
            multi += ("\n".join(m).rstrip() + "\n")
            multi += "ENDMDL\n"
    
        return multi, len(models)

        # ---------- NEW: STEP_*_ANMLD.pdb -> trajectory helpers ----------
    _STEP_RE = re.compile(r"STEP[_-]?(\d+)", re.IGNORECASE)

    def _step_index(p: Path) -> int:
        m = _STEP_RE.search(p.name)
        return int(m.group(1)) if m else 10**9

    def _collect_step_pdbs(out_dir: Path, mode: int) -> list[Path]:
        """
        Mode'a ait STEP frame'lerini olabildiğince sağlam yakala.
        Önce mode'a özel klasör/pattern dener, bulamazsa tüm STEP'lere fallback yapar.
        """
        patterns = [
            f"**/mode{mode}/STEP_*_ANMLD.pdb",
            f"**/MODE{mode}/STEP_*_ANMLD.pdb",
            f"**/*mode{mode}*/STEP_*_ANMLD.pdb",
            f"**/*MODE{mode}*/STEP_*_ANMLD.pdb",
            f"**/*moved{mode}*/STEP_*_ANMLD.pdb",
        ]

        found: list[Path] = []
        for pat in patterns:
            found = list(out_dir.glob(pat))
            if found:
                break

        if not found:
            # son çare: her şeyi topla
            found = list(out_dir.rglob("STEP_*_ANMLD.pdb"))

        # unique + STEP sırasına göre sort
        uniq = sorted({p.resolve() for p in found}, key=lambda x: (_step_index(x), str(x)))
        return uniq

    def _concat_steps_as_models(pdb_files: list[Path], max_frames: int = 200) -> tuple[str, int]:
        """
        STEP dosyalarını MODEL/ENDMDL ile birleştir.
        Çok fazlaysa stride ile seyrekleştir.
        """
        if not pdb_files:
            return "", 0

        files = pdb_files
        if len(files) > max_frames:
            stride = max(1, len(files) // max_frames)
            files = files[::stride]

        chunks: list[str] = []
        for i, f in enumerate(files, start=1):
            chunks.append(f"MODEL        {i}\n")
            # mevcut helper'ını kullanıyoruz: MODEL/ENDMDL/END satırlarını temizler
            chunks.append(_read_pdb_for_frames(str(f)))
            chunks.append("ENDMDL\n")

        return "".join(chunks), len(files)


    def _make_mode_viewer_from_file(mode_pdb_path: str):
        """
        Colab'da: /files + iframe (mevcut stabil yöntem)
        Localhost'ta: data:text/html;base64 + iframe (artık /files/content/... sorunu yok)
        """
        py3Dmol = _ensure_py3dmol()

        MODE_W = 560
        MODE_H = 280

        holder = W.HTML(
            value="<div style='font-family:Arial;color:#6b7280;'>Rendering 3D view…</div>",
            layout=W.Layout(width=f"{MODE_W}px", height=f"{MODE_H}px"),
        )

        try:
            pdb_text = Path(mode_pdb_path).read_text(encoding="utf-8", errors="ignore")
            if not pdb_text.strip():
                raise RuntimeError("Empty PDB text.")

            multi_pdb, nmodels = _build_multimodel_pdb_string(pdb_text)
            bmin, bmax = _bfactor_minmax(multi_pdb)

            view = py3Dmol.view(width=MODE_W, height=MODE_H)
            view.setBackgroundColor("white")
            view.addModelsAsFrames(multi_pdb, "pdb")

            style = {
                "cartoon": {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": float(bmin),
                        "max": float(bmax),
                    }
                }
            }
            view.setStyle({}, style)
            view.zoomTo()

            if nmodels >= 2:
                view.animate({"loop": "backAndForth", "reps": 0, "interval": 180})

            raw = _html_with_unique_divid(view._make_html())
            doc = _wrap_html_doc(raw)
            holder.value = _iframe_from_html_doc(doc, MODE_W, MODE_H)

        except Exception as e:
            holder.value = (
                "<div style='font-family:Arial;color:#dc2626;font-weight:800;'>"
                f"Viewer error: {_safe_html(str(e))}"
                "</div>"
            )

        return W.HBox([holder], layout=W.Layout(width="100%", justify_content="center", align_items="center"))

    def _make_mode_viewer_from_steps(out_dir: str, mode: int, fallback_mode_pdb: str | None = None):
        """
        1) Önce STEP_*_ANMLD.pdb framelerini bul -> animasyon göster
        2) Bulamazsa fallback olarak mode1/mode2 pdb dosyasını göster
        """
        py3Dmol = _ensure_py3dmol()

        MODE_W = 560
        MODE_H = 280

        holder = W.HTML(
            value="<div style='font-family:Arial;color:#6b7280;'>Rendering mode trajectory…</div>",
            layout=W.Layout(width=f"{MODE_W}px", height=f"{MODE_H}px"),
        )

        try:
            od = Path(out_dir)
            if not od.exists():
                raise RuntimeError(f"Output dir not found: {out_dir}")

            step_files = _collect_step_pdbs(od, mode)

            # STEP frameleri gerçekten varsa onları kullan
            if len(step_files) >= 2:
                multi_pdb, nmodels = _concat_steps_as_models(step_files)

                v = py3Dmol.view(width=MODE_W, height=MODE_H)
                v.setBackgroundColor("white")
                v.addModelsAsFrames(multi_pdb, "pdb")
                v.setStyle({"cartoon": {"color": "spectrum"}})
                v.zoomTo()
                v.animate({"loop": "backAndForth", "reps": 0, "interval": 180})

                raw = _html_with_unique_divid(v._make_html())
                doc = _wrap_html_doc(raw)
                holder.value = _iframe_from_html_doc(doc, MODE_W, MODE_H)

            else:
                # frameler yoksa fallback
                if fallback_mode_pdb and os.path.exists(fallback_mode_pdb) and os.path.getsize(fallback_mode_pdb) > 0:
                    return _make_mode_viewer_from_file(fallback_mode_pdb)
                raise RuntimeError("STEP_*_ANMLD.pdb frameleri bulunamadı (ve fallback mode pdb de yok).")

        except Exception as e:
            holder.value = (
                "<div style='font-family:Arial;color:#dc2626;font-weight:800;'>"
                f"Mode viewer error: {_safe_html(str(e))}"
                "</div>"
            )

        return W.HBox([holder], layout=W.Layout(width="100%", justify_content="center", align_items="center"))

    def _make_mode_viewer(mode: int):
        """
        Tek giriş noktası: state'den out_dir + fallback mode.pdb alıp STEP'ten göstermeyi dener.
        """
        out_dir = state.get("mode_steps_dir") or state.get("last_out_dir")
        mf = state.get("mode_files", {}) or {}
        return _make_mode_viewer_from_steps(out_dir, int(mode), mf.get(int(mode)))
    
    
    
    # chain colors (deterministic per detected order)
    _CHAIN_PALETTE = [
        "red", "blue", "green", "orange", "purple", "cyan", "magenta",
        "yellow", "teal", "brown", "pink", "lime", "navy", "gold"
    ]

    def _assign_chain_colors(chains: list[str]) -> dict[str, str]:
        cmap: dict[str, str] = {}
        for i, ch in enumerate(chains):
            cmap[ch] = _CHAIN_PALETTE[i % len(_CHAIN_PALETTE)]
        return cmap

    # ---------- UI ----------
    css = W.HTML(r"""
    <style>
    .hp-card {
      border:1px solid #e5e7eb;
      border-radius:14px;
      padding:14px 16px;
      margin:0;                 /* IMPORTANT: row alignment */
      background:#fff;
      box-sizing:border-box;    /* IMPORTANT: same width actually means same width */
    }
    .hp-headerbar{
      border:1px solid #e5e7eb;
      border-radius:16px;
      padding:12px 18px;
      margin:10px 0 12px 0;
      background:#fff;
      box-shadow: 0 1px 0 rgba(0,0,0,0.03);
      display:flex;
      flex-direction:column;
      align-items:flex-start;
      gap:10px;
      box-sizing:border-box;
    }
    .hp-brand{
      display:flex;
      flex-direction:column;
      align-items:flex-start;
    }
    .hp-logo{
      max-height:78px;
      height:auto;
      display:block;
      margin: 0 0 6px 0;
    }
    .hp-subtitle{
      font-size:16px;
      font-weight:800;
      color:#111827;
      font-family: Arial, Helvetica, sans-serif;
      margin-top:0px;
      line-height:1.2;
      text-align:left;
    }
    .hp-pre{
      white-space:pre-wrap;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
      font-size: 13px;
      line-height: 1.35;
      background:#0b1020;
      color:#e5e7eb;
      padding:12px;
      border-radius:12px;
      border:1px solid #1f2937;
    }

    .hp-about{
      font-family: Arial, Helvetica, sans-serif;
      font-size: 13px;
      line-height: 1.35;
      color: #111827;
      max-width: none;
      width: 100%
    }
    .hp-about h2{
      margin: 6px 0 6px 0;
      font-size: 18px;
      font-weight: 800;
      color: #1d4ed8; /* blue */
    }
    .hp-about h3{
      margin: 14px 0 6px 0;
      font-size: 16px;
      font-weight: 800;
      color: #1d4ed8; /* blue */
    }
    .hp-about p{ margin: 6px 0; }
    .hp-about .hp-ref{
      color:#dc2626; /* red */
      text-decoration: underline;
      font-weight: 700;
    }
    .hp-figrow{
      display:flex;
      justify-content:center;
      align-items:center;
      gap: 18px;
      flex-wrap: wrap;
      margin: 12px 0 10px 0;
    }
    .hp-figbox{
      background:#000;
      padding:6px;
      border-radius:10px;
      border:1px solid #111827;
    }
    .hp-figbox img{
      display:block;
      max-width: 210px;
      width: min(210px, 46vw);
      height: auto;
      border-radius:6px;
    }
    .hp-about .hp-refs{
      margin-top: 10px;
      padding-top: 8px;
      border-top: 1px solid #e5e7eb;
    }
    .hp-about .hp-refs ol{
      margin: 6px 0 0 20px;
      padding: 0;
    }
    .hp-about .hp-refs li{
      margin: 6px 0;
    }

    /* ---------- NEW: References page style (screenshot-like) ---------- */
    .hp-refpage{
      font-family: Arial, Helvetica, sans-serif;
      width: 100%;
      color:#111827;
    }
    .hp-refpage h2{
      margin: 0 0 8px 0;
      font-size: 28px;
      font-weight: 900;
      color: #1d4ed8;
    }
    .hp-refpage .sub{
      margin: 0 0 10px 0;
      font-size: 14px;
      font-weight: 700;
    }
    .hp-refpage .cite-main{
      margin: 0 0 16px 0;
      font-size: 14px;
    }
    .hp-refpage h3{
      margin: 14px 0 8px 0;
      font-size: 22px;
      font-weight: 900;
      color: #1d4ed8;
    }
    .hp-refpage ul{
      margin: 8px 0 0 22px;
      padding: 0;
    }
    .hp-refpage li{
      margin: 6px 0;
      line-height: 1.25;
    }
    
    .hp-refpage a{
      color:#1d4ed8;
      text-decoration: underline;
      font-weight: 700;
    }

    .hp-note{
      margin-top: 14px;
      font-size: 12px;
      color: #6b7280;
    }

    .hp-modebar{
      position:relative;
      display:flex;
      align-items:center;
      justify-content:center;
      gap:10px;
      margin:10px 0 6px 0;
    }
 
    .hp-modetitle{
      width:100%;
      text-align:center;
    }
    .hp-modebar .hp-dlbtn{
      position:absolute;
      right:0;
    }
    
    .hp-dlbtn{
      display:inline-block;
      padding:6px 10px;
      border-radius:10px;
      border:1px solid #e5e7eb;
      background:#ffffff;
      color:#111827;
      font-size:12px;
      font-weight:800;
      text-decoration:none;
      white-space:nowrap;
    }
    .hp-dlbtn:hover{
      border-color:#93c5fd;
      box-shadow:0 1px 0 rgba(0,0,0,0.04);
    }


    </style>
    """)

    logo_url = "https://raw.githubusercontent.com/enesemretas/hingeprot_fortran/main/assets/logo.gif"

    brand = W.HTML(f"""
    <div class="hp-brand">
      <img class="hp-logo" src="{logo_url}" alt="HINGEprot logo">
      <div class="hp-subtitle">An Algorithm For Protein Hinge Prediction Using Elastic Network Models</div>
    </div>
    """)

    nav = W.ToggleButtons(
        options=[
            ("Web Server", "web"),
            ("About HingeProt", "about"),
            ("Help", "help"),
            ("References", "refs"),
        ],
        value="web",
        layout=W.Layout(width="auto"),
        style={"button_width": "140px"},
    )

    # Nav: subtitle'ın ALTINDA
    header_bar = W.VBox(
        [brand, nav],
        layout=W.Layout(width="100%", align_items="flex-start", justify_content="flex-start", gap="8px"),
    )
    header_bar.add_class("hp-headerbar")

    # --- bold Input label ---
    input_label = W.HTML("<b>Input:</b>", layout=W.Layout(width="60px"))
    input_mode = W.ToggleButtons(
        options=[("Enter PDB code", "code"), ("Upload PDB file", "upload")],
        value="code",
        description="",  # handled by input_label
        style={"description_width": "0px", "button_width": "170px"},
        layout=W.Layout(width="420px"),
    )
    input_row = W.HBox([input_label, input_mode], layout=W.Layout(align_items="center", gap="10px"))

    pdb_code = W.Text(
        value="",
        description="PDB code:",
        placeholder="e.g., 4cln",
        style={"description_width": "80px"},
        layout=W.Layout(width="420px"),
    )

    btn_choose_file = W.Button(description="Choose file", icon="upload", layout=W.Layout(width="180px"))
    file_upload = W.FileUpload(accept=".pdb,.ent", multiple=False, layout=W.Layout(width="260px"))  # LOCAL
    upload_prog = W.IntProgress(value=0, min=0, max=100, description="", layout=W.Layout(width="160px"))
    upload_prog.bar_style = ""
    file_lbl = W.Label("No file chosen")

    code_box = W.HBox([pdb_code], layout=W.Layout(align_items="center"))
    
    # Upload box: Colab -> button, Local -> FileUpload widget
    if IS_COLAB:
        upload_box = W.HBox([btn_choose_file, upload_prog, file_lbl], layout=W.Layout(align_items="center", gap="10px"))
    else:
        upload_box = W.HBox([file_upload, upload_prog, file_lbl], layout=W.Layout(align_items="center", gap="10px"))


    btn_load = W.Button(description="Load / Detect Chains", button_style="info", icon="search", layout=W.Layout(width="260px"))

    all_chains = W.Checkbox(
        value=False,
        description="All Chains",
        indent=False,
        style={"description_width": "initial"},
        layout=W.Layout(width="120px", min_width="120px", flex="0 0 120px"),
    )

    chains_label = W.HTML("<b>Select Chains:</b>", layout=W.Layout(width="120px"))
    chains_wrap = W.Box(
        [],
        layout=W.Layout(
            flex="1 1 auto",
            width="auto",
            min_width="220px",
            display="flex",
            flex_flow="row wrap",
            align_items="center",
            gap="10px",
            border="1px solid #e5e7eb",
            border_radius="12px",
            padding="8px 10px",
        ),
    )
    chain_row = W.HBox([chains_label, chains_wrap], layout=W.Layout(align_items="center", gap="12px", width="100%"))

    gnm_row, get_gnm_cut = _list_or_custom_float("GNM cutoff (Å):", [7, 8, 9, 10, 11, 12, 13, 20], 10.0, 1.0, 100.0)
    anm_row, get_anm_cut = _list_or_custom_float("ANM cutoff (Å):", [10, 13, 15, 18, 20, 23, 36], 18.0, 1.0, 100.0)

    progress = W.IntProgress(value=0, min=0, max=4, description="Progress:", bar_style="")
    btn_run_fortran = W.Button(description="Run HingeProt", button_style="primary", icon="play", layout=W.Layout(width="320px"))
    btn_clear = W.Button(description="Clear", button_style="warning", icon="trash", layout=W.Layout(width="180px"))

    table_box = W.VBox([], layout=W.Layout(width="100%", gap="8px"))

    status_box = W.HTML('<div class="hp-pre">Load a PDB to detect chains.</div>')

    def _set_status(text: str):
        status_box.value = f'<div class="hp-pre">{_safe_html(text)}</div>'

    def _download_file(path: str):
        path = os.path.abspath(path)
        if (not os.path.exists(path)) or os.path.getsize(path) == 0:
            _set_status(f"ERROR: file not found or empty:\n{path}")
            return
        # Colab'ta doğrudan indir
        try:
            from google.colab import files  # type: ignore
            files.download(path)
            return
        except Exception:
            # Local Jupyter vb. için link göster
            display(FileLink(path))


    # ---------- viewer (boxed like form; same width & stretchable height) ----------
    CARD_W = 620
    CARD_PAD = 16  # hp-card left+right padding
    OUT_PAD = 6    # viewer_out padding
    OUT_BORDER = 1
    INNER_W = CARD_W - 2 * CARD_PAD

    # IMPORTANT: avoid horizontal scrollbar by making py3Dmol slightly smaller than inner area
    VIEW_W = INNER_W - 2 * (OUT_PAD + OUT_BORDER)
    VIEW_H = 280

    viewer_out = W.Output(
        layout=W.Layout(
            width="100%",
            height=f"{VIEW_H}px",
            border="1px solid #e5e7eb",
            border_radius="12px",
            padding=f"{OUT_PAD}px",
            overflow="hidden",
        )
    )
    viewer_title = W.HTML('<div style="font-weight:800; margin:2px 0 8px 2px;">3D Viewer</div>')

    viewer_card = W.VBox(
        [viewer_title, viewer_out],
        layout=W.Layout(
            width=f"{CARD_W}px",
            gap="10px",
            height="100%",            # allow stretch
            display="flex",
            flex_flow="column",
            align_items="stretch",
        ),
    )
    viewer_card.add_class("hp-card")

    def _viewer_placeholder(msg: str = "Load a PDB to preview it here."):
        with viewer_out:
            clear_output(wait=True)
            print(msg)

    _viewer_placeholder()

    # ---------- state ----------
    state = {
        "pdb_text": None,
        "pdb_filename": None,
        "pdb_path": None,
        "run_dir": None,
        "pdb_tag": None,
        "upload_name": None,
        "upload_bytes": None,
        "detected_chains": [],
        "chain_cbs": {},
        "chain_colors": {},     # NEW
        "manual_selection": (),
        "_syncing": False,
        "hingeprot_dir": None,
        "last_out_dir": None,
        "mode_files": {},   # {1: "/path/...mode1.pdb", 2: "/path/...mode2.pdb"}
        "mode_steps_dir": None,  # NEW: STEP_* frame'lerinin aranacağı klasör (dest_out_dir)

    }
    global LAST_UI_STATE
    LAST_UI_STATE = state

    # ---------- chain selection ----------
    def _selected_chains() -> list[str]:
        detected = state.get("detected_chains", [])
        return [ch for ch in detected if ch in state["chain_cbs"] and state["chain_cbs"][ch].value]

    def _set_selection(sel: list[str]):
        detected = state.get("detected_chains", [])
        sel = [c for c in sel if c in detected]
        state["_syncing"] = True
        try:
            for ch, cb in state["chain_cbs"].items():
                cb.value = (ch in sel)
        finally:
            state["_syncing"] = False

    def _update_all_checkbox_from_selection():
        if state["_syncing"]:
            return
        detected = state.get("detected_chains", [])
        if not detected:
            return
        sel = _selected_chains()
        all_now = (len(sel) == len(detected)) and (len(detected) > 0)

        state["_syncing"] = True
        try:
            all_chains.value = all_now
        finally:
            state["_syncing"] = False

        if not all_now:
            state["manual_selection"] = tuple(sel)

    def _refresh_viewer():
        if not state.get("pdb_text"):
            _viewer_placeholder()
            return

        py3Dmol = _ensure_py3dmol()
        pdb_text = state["pdb_text"]

        detected = state.get("detected_chains", [])
        if detected:
            if all_chains.value:
                selected = list(detected)
            else:
                selected = _selected_chains()
        else:
            selected = []

        chain_colors: dict[str, str] = state.get("chain_colors", {}) or {}

        v = py3Dmol.view(width=VIEW_W, height=VIEW_H)
        v.addModel(pdb_text, "pdb")
        v.setBackgroundColor("white")

        # default: everything grey
        v.setStyle({}, {"cartoon": {"color": "lightgray"}})

        # selected: each chain different color
        for ch in selected:
            col = chain_colors.get(ch, "red")
            v.setStyle({"chain": ch}, {"cartoon": {"color": col}})

        # ligands / hetero atoms (exclude waters)
        # show as sticks+spheres
        v.setStyle(
            {"hetflag": True, "not": {"resn": ["HOH", "WAT", "DOD"]}},
            {"stick": {}, "sphere": {"scale": 0.25}},
        )

        v.zoomTo()

        raw = _html_with_unique_divid(v._make_html())
        doc = _wrap_html_doc(raw)
        iframe_html = _iframe_from_html_doc(doc, VIEW_W, VIEW_H)

        with viewer_out:
            clear_output(wait=True)
            display(W.HTML(iframe_html))



    # ---------- input visibility ----------
    def _sync_input_visibility(*_):
        if input_mode.value == "code":
            code_box.layout.display = ""
            upload_box.layout.display = "none"
        else:
            code_box.layout.display = "none"
            upload_box.layout.display = ""

    _sync_input_visibility()
    input_mode.observe(lambda ch: _sync_input_visibility(), names="value")

    # ---------- uploader callbacks ----------
    cb_name = f"hingeprot_uploader_{uuid.uuid4().hex}"
    cb_prog = f"hingeprot_uploadprog_{uuid.uuid4().hex}"

    def _js_upload_progress_callback(payload):
        try:
            pct = int(payload.get("pct", 0))
            pct = max(0, min(100, pct))
            upload_prog.value = pct
            upload_prog.bar_style = "info" if pct < 100 else "success"
        except Exception:
            pass

    def _js_upload_callback(payload):
        try:
            name = payload.get("name", "upload.pdb")
            data_b64 = payload.get("data_b64", "")
            if not data_b64:
                _set_status("Upload callback received empty data.")
                return
            data = base64.b64decode(data_b64.encode("utf-8"))
            state["upload_name"] = name
            state["upload_bytes"] = data
            file_lbl.value = name
            upload_prog.value = 100
            upload_prog.bar_style = "success"
            _set_status(f"Uploaded file: {name} ({len(data)} bytes)\nNow click 'Load / Detect Chains'.")
        except Exception as e:
            _set_status(f"Upload callback error: {e}")



    # ---- COLAB upload wiring ----
    if IS_COLAB and colab_output is not None:
        colab_output.register_callback(cb_prog, _js_upload_progress_callback)
        colab_output.register_callback(cb_name, _js_upload_callback)

        def on_choose_file(_):
            upload_prog.value = 0
            upload_prog.bar_style = "info"

            js = f"""
            (async () => {{
              const input = document.createElement('input');
              input.type = 'file';
              input.accept = '.pdb,.ent';
              input.style.display = 'none';
              document.body.appendChild(input);

              input.onchange = async () => {{
                const file = input.files && input.files[0];
                document.body.removeChild(input);
                if (!file) return;

                const reader = new FileReader();

                reader.onloadstart = async () => {{
                  try {{
                    await google.colab.kernel.invokeFunction("{cb_prog}", [{{pct: 0}}], {{}});
                  }} catch (err) {{}}
                }};

                reader.onprogress = async (e) => {{
                  try {{
                    if (e.lengthComputable) {{
                      const pct = Math.round((e.loaded / e.total) * 100);
                      await google.colab.kernel.invokeFunction("{cb_prog}", [{{pct: pct}}], {{}});
                    }}
                  }} catch (err) {{}}
                }};

                reader.onloadend = async () => {{
                  try {{
                    await google.colab.kernel.invokeFunction("{cb_prog}", [{{pct: 100}}], {{}});
                  }} catch (err) {{}}
                }};

                reader.onload = async () => {{
                  const b64 = (reader.result || "").split(",")[1] || "";
                  await google.colab.kernel.invokeFunction(
                    "{cb_name}",
                    [{{name: file.name, data_b64: b64}}],
                    {{}}
                  );
                }};

                reader.readAsDataURL(file);
              }};

              input.click();
            }})();
            """
            colab_output.eval_js(js)

        btn_choose_file.on_click(on_choose_file)

    # ---- LOCAL upload wiring (localhost) ----
    if not IS_COLAB:
        def _on_local_upload(change):
            try:
                val = file_upload.value
                if not val:
                    return

                # ipywidgets v8: dict{filename: {content: bytes, metadata:{name,...}}}
                if isinstance(val, dict):
                    # filename key en sağlam
                    fname_key = next(iter(val.keys()))
                    item = val[fname_key] or {}
                    data = item.get("content") or item.get("data") or b""
                    meta = item.get("metadata") or {}
                    name = meta.get("name") or item.get("name") or fname_key or "upload.pdb"
                else:
                    # bazı ortamlarda list/tuple gelebilir
                    item = val[0] if isinstance(val, (list, tuple)) and val else {}
                    data = item.get("content") or item.get("data") or b""
                    name = (item.get("metadata") or {}).get("name") or item.get("name") or "upload.pdb"

                if not data:
                    _set_status("Local upload received empty file.")
                    return

                state["upload_name"] = name
                state["upload_bytes"] = bytes(data)
                file_lbl.value = name
                upload_prog.value = 100
                upload_prog.bar_style = "success"
                _set_status(f"Uploaded file: {name} ({len(state['upload_bytes'])} bytes)\nNow click 'Load / Detect Chains'.")
            except Exception as e:
                _set_status(f"Local upload error: {e}")

        file_upload.observe(_on_local_upload, names="value")

    # ---------- chain checkbox rebuild ----------
    def _on_chain_cb_change(_):
        _update_all_checkbox_from_selection()
        _refresh_viewer()

    def _rebuild_chain_checkboxes(chains: list[str], default_selected: list[str]):
        state["chain_cbs"] = {}
        items = [all_chains]
        for ch in chains:
            cb = W.Checkbox(
                value=(ch in default_selected),
                description=ch,
                indent=False,
                layout=W.Layout(width="48px", flex="0 0 48px"),
            )
            cb.observe(_on_chain_cb_change, names="value")
            state["chain_cbs"][ch] = cb
            items.append(cb)
        chains_wrap.children = items

    def _on_all_chains_toggle(ch):
        if state["_syncing"]:
            return
        detected = state.get("detected_chains", [])
        if not detected:
            return

        if ch["new"] is True:
            sel = _selected_chains()
            if len(sel) != len(detected):
                state["manual_selection"] = tuple(sel)
            _set_selection(detected)
        else:
            prev = list(state.get("manual_selection") or [])
            prev = [c for c in prev if c in detected]
            if not prev:
                prev = [detected[0]]
            _set_selection(prev)

        _update_all_checkbox_from_selection()
        _refresh_viewer()

    all_chains.observe(_on_all_chains_toggle, names="value")

    # ---------- actions ----------
    def on_load_clicked(_):
        progress.value = 0
        table_box.children = ()
        progress.bar_style = "info"
        state["last_out_dir"] = None

        try:
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

            if input_mode.value == "upload":
                if state["upload_bytes"] is None:
                    raise ValueError("Please click 'Choose file' and upload a PDB first.")
                pdb_text = state["upload_bytes"].decode("utf-8", errors="ignore")

                upname = (state.get("upload_name") or "upload.pdb").strip()
                if not re.search(r"\.(pdb|ent)$", upname, flags=re.I):
                    upname = upname + ".pdb"
                pdb_filename = os.path.basename(upname)

                base = os.path.splitext(pdb_filename)[0]
                tag = re.sub(r"[^0-9A-Za-z]+", "", base).upper() or "UPLOAD"
            else:
                code = pdb_code.value.strip()
                if not code:
                    raise ValueError("Please enter a PDB code (e.g., 3lzg).")
                code4 = code.strip().upper()
                pdb_text = _fetch_pdb_text(code4)
                pdb_filename = f"{code4.lower()}.pdb"
                tag = code4

            run_dir = os.path.join(runs_root, f"{tag}_run_{ts}")
            os.makedirs(run_dir, exist_ok=True)

            state["run_dir"] = run_dir
            state["pdb_text"] = pdb_text
            state["pdb_filename"] = pdb_filename
            state["pdb_tag"] = tag

            pdb_path = os.path.join(run_dir, pdb_filename)
            with open(pdb_path, "w", encoding="utf-8") as f:
                f.write(pdb_text)
            state["pdb_path"] = pdb_path

            chs = _detect_chains_from_text(pdb_text)
            if not chs:
                raise RuntimeError("No chains detected in the PDB.")
            state["detected_chains"] = chs

            # NEW: deterministic chain colors
            state["chain_colors"] = _assign_chain_colors(chs)

            default_sel = [chs[0]]
            state["manual_selection"] = tuple(default_sel)
            _rebuild_chain_checkboxes(chs, default_sel)

            state["_syncing"] = True
            try:
                all_chains.value = False
            finally:
                state["_syncing"] = False

            progress.value = 1
            cmap = state.get("chain_colors", {})
            preview = ", ".join([f"{c}:{cmap.get(c,'')}" for c in chs[:6]])
            if len(chs) > 6:
                preview += ", ..."

            _set_status(
                f"Loaded PDB (ID={tag})\n"
                f"Run folder: {run_dir}\n"
                f"Detected chains: {', '.join(chs)}\n"
                f"Chain colors: {preview}\n\n"
                "Viewer: unselected chains are GREY; selected chains are COLORED. "
                "Ligands (HETATM; excluding waters) are shown as sticks/spheres."
            )

            _refresh_viewer()

        except Exception as e:
            progress.bar_style = "danger"
            _set_status(f"ERROR: {e}")

    def _capture_inputs() -> dict:
        if not state.get("pdb_text"):
            raise RuntimeError("Please click 'Load / Detect Chains' first.")

        detected = state.get("detected_chains", [])
        if not detected:
            raise RuntimeError("No detected chains. Load again.")

        if all_chains.value:
            chain_list = detected
        else:
            chain_list = _selected_chains()
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
        progress.bar_style = "info"
        try:
            captured = _capture_inputs()
            global LAST_INPUTS
            LAST_INPUTS = captured
            table_box.children = (W.HTML("<div style='font-family:Arial; color:#6b7280;'>Running HingeProt…</div>"),)
            progress.value = 1
            _ensure_libg2c()

            progress.value = 2
            hp_dir = _ensure_repo(fresh=False)
            state["hingeprot_dir"] = hp_dir

            progress.value = 3
            _write_runHingeProt_pl(hp_dir, captured["gnm_cutoff_A"], captured["anm_cutoff_A"])

            pdb_filename = captured["pdb_filename"]
            pdb_abs = os.path.join(hp_dir, pdb_filename)
            with open(pdb_abs, "w", encoding="utf-8") as f:
                f.write(state["pdb_text"] or "")

            chains_str = captured["chains_str"]
            out_dir_repo = os.path.join(hp_dir, f"{pdb_filename}.{chains_str}")

            if os.path.isdir(out_dir_repo):
                shutil.rmtree(out_dir_repo, ignore_errors=True)

            cmd = f"perl ./runHingeProt.pl {pdb_filename} {chains_str}"
            r = _sh(cmd, cwd=hp_dir, timeout=1200)  # 20 dk
            if r.returncode != 0:
                raise RuntimeError(f"runHingeProt.pl failed (return code {r.returncode}).\n{r.stderr}")

            run_dir = captured["run_dir_runsroot"]
            if not run_dir or not os.path.isdir(run_dir):
                raise RuntimeError("Run folder not found. Please 'Load / Detect Chains' again.")

            dest_out_dir = os.path.join(run_dir, os.path.basename(out_dir_repo))
            if os.path.isdir(dest_out_dir):
                shutil.rmtree(dest_out_dir, ignore_errors=True)

            if os.path.isdir(out_dir_repo):
                shutil.move(out_dir_repo, dest_out_dir)
            else:
                raise RuntimeError(f"Expected output folder not found: {out_dir_repo}")

            state["last_out_dir"] = dest_out_dir

            state["mode_steps_dir"] = dest_out_dir   # NEW: STEP frameleri burada aranacak

            

            # NEW: Outputs altına bilgilendirme satırı
            fname = captured["pdb_filename"] or "PDB"
            chains_str = captured["chains_str"] or ""
            chains_pretty = ", ".join(list(chains_str)) if chains_str else "-"
            output_info.value = (
                "<div style='text-align:left; font-family:Arial, Helvetica, sans-serif;"
                "font-weight:800; color:#111827; margin:2px 0 6px 0;'>"
                f"{_safe_html(fname)} for Chain(s) : {_safe_html(chains_pretty)}"
                "</div>"
            )

            # ---------- NEW: fill SINGLE Mode Viewer ----------
            mode_files = {}
            for m in (1, 2):
                p = os.path.join(dest_out_dir, f"{pdb_filename}.mode{m}.pdb")
                if os.path.exists(p) and os.path.getsize(p) > 0:
                    mode_files[m] = p

            _set_mode_viewer_files(mode_files)



            # ---------- NEW: build table from .new.hinges report (NO calculation) ----------

            pdb_chain_path = Path(dest_out_dir) / "pdb"
            if not pdb_chain_path.exists():
                pdb_chain_path = Path(dest_out_dir) / pdb_filename
            if not pdb_chain_path.exists():
                raise RuntimeError(f"Chain PDB file not found (expected 'pdb' or '{pdb_filename}') in: {dest_out_dir}")

            # 1) REPORT: *.new.hinges
            report_file = _find_new_hinges_report(dest_out_dir, pdb_filename)
            if not report_file:
                raise RuntimeError(f"Report file not found: expected '{pdb_filename}.new.hinges' in {dest_out_dir}")

            report = parse_new_hinges_report(Path(report_file))

            if not report:
                raise RuntimeError("Parsed report is empty. Check *.new.hinges content.")

            # 2) RAW HINGES: *.hinge (for Short Flexible Fragments only)
            hinge_path = Path(dest_out_dir) / f"{pdb_filename}.hinge"
            if (not hinge_path.exists()) or (hinge_path.stat().st_size == 0):
                alt = _find_hinges_file(dest_out_dir, pdb_filename)  # sende vardı; .hinge bulur
                if alt:
                    hinge_path = Path(alt)
                else:
                    raise RuntimeError(f".hinge file not found: {hinge_path}")

            residues_by_chain = read_residue_order_from_pdb(pdb_chain_path)
            hinge_modes = parse_hinge_file(hinge_path)

            # same min_len spirit as before
            min_len = 14

            short_frags_by_mode = compute_short_flexible_fragments(
                residues_by_chain=residues_by_chain,
                hinge_modes=hinge_modes,
                report=report,
                min_len=min_len,
            )

            rigid_widget = rigidparts_report_widget_from_report(
                pdb_label=pdb_filename,
                report=report,
                short_frags_by_mode=short_frags_by_mode,
                out_dir=dest_out_dir,
                download_fn=_download_file,
            )

            table_box.children = (
                mode_viewer_card,
                rigid_widget,
            )





            progress.value = 4
            progress.bar_style = "success"

        except Exception as e:
            progress.bar_style = "danger"
            table_box.children = (W.HTML("<div style='color:#dc2626;font-weight:800;'>ERROR</div>"),)
            _set_status(f"ERROR: {e}")

    def on_clear_clicked(_):
        pdb_code.value = ""
        table_box.children = ()
        input_mode.value = "code"
        state["upload_name"] = None
        state["upload_bytes"] = None
        file_lbl.value = "No file chosen"
        upload_prog.value = 0
        upload_prog.bar_style = ""

        state["pdb_text"] = None
        state["pdb_filename"] = None
        state["pdb_path"] = None
        state["run_dir"] = None
        state["pdb_tag"] = None
        state["hingeprot_dir"] = None
        state["last_out_dir"] = None
        state["mode_steps_dir"] = None

        state["detected_chains"] = []
        state["chain_cbs"] = {}
        state["chain_colors"] = {}
        state["manual_selection"] = ()
        state["_syncing"] = False
        all_chains.value = False
        chains_wrap.children = ()
        output_info.value = ""
        
        # reset SINGLE Mode Viewer
        state["mode_files"] = {}
        mode_viewer_card.layout.display = "none"
        mode_select.disabled = True
        mode_dl.disabled = True
        _mode_placeholder()

        

        global LAST_INPUTS
        LAST_INPUTS = None

        progress.value = 0
        progress.max = 4
        progress.bar_style = ""
        _set_status("Cleared. Load a PDB to detect chains.")
        _viewer_placeholder()

    btn_load.on_click(on_load_clicked)
    btn_run_fortran.on_click(on_run_fortran_clicked)
    btn_clear.on_click(on_clear_clicked)

    # --------- CARDS (same width) ---------
    form_card = W.VBox(
        [
            input_row,          # NEW (bold Input:)
            code_box,
            upload_box,
            btn_load,
            W.HTML("<hr>"),
            chain_row,
            W.VBox([gnm_row, anm_row], layout=W.Layout(gap="8px")),
            progress,
            W.HBox([btn_run_fortran, btn_clear]),
        ],
        layout=W.Layout(width=f"{CARD_W}px", gap="10px", height="100%"),  # allow stretch
    )
    form_card.add_class("hp-card")

    output_title = W.HTML("<b>Outputs</b>")

    # NEW: run sonrası bilgi satırı
    output_info = W.HTML("", layout=W.Layout(width="100%"))

    # ---------- Outputs: SINGLE Mode Viewer (updates on Mode 1/2) ----------
    mode_title = W.HTML("<b>Mode Viewer</b>")

    mode_select = W.ToggleButtons(
        options=[("Mode 1", 1), ("Mode 2", 2)],
        value=1,
        style={"button_width": "120px"},
        layout=W.Layout(width="260px"),
    )
    mode_select.disabled = True

    mode_dl = W.Button(
        description="Download current mode.pdb",
        icon="download",
        layout=W.Layout(width="260px"),
    )
    mode_dl.disabled = True

    mode_view_container = W.Output(
        layout=W.Layout(
            width="100%",
            height="320px",
            border="1px solid #e5e7eb",
            border_radius="12px",
            padding="6px",
            overflow="hidden",
        )
    )

    # Bu kart output içinde tek bir viewer alanı olacak
    mode_viewer_card = W.VBox(
        [
            mode_title,
            W.HBox([mode_select, mode_dl], layout=W.Layout(gap="10px", align_items="center")),
            mode_view_container,
        ],
        layout=W.Layout(width="100%", gap="8px"),
    )
    mode_viewer_card.add_class("hp-card")
    mode_viewer_card.layout.display = "none"  # ilk başta gizli

    def _mode_placeholder(msg: str = "Run HingeProt to generate mode files."):
        with mode_view_container:
            clear_output(wait=True)
            print(msg)

    _mode_placeholder()

    def _render_mode(m: int):
        with mode_view_container:
            clear_output(wait=True)
            display(_make_mode_viewer(int(m)))
            
    def _set_mode_viewer_files(mode_files: dict[int, str]):
        """
        Run bitince çağrılacak:
        - hangi mode dosyaları var -> options
        - seçili mode -> viewer render
        - butonları enable et
        """
        state["mode_files"] = mode_files or {}

        if not mode_files:
            mode_viewer_card.layout.display = "none"
            mode_select.disabled = True
            mode_dl.disabled = True
            _mode_placeholder("No mode PDB produced.")
            return

        opts = [(f"Mode {m}", m) for m in sorted(mode_files.keys())]
        mode_select.options = opts
        mode_select.value = opts[0][1]          # ilk mevcut mode
        mode_select.disabled = False
        mode_dl.disabled = False
        mode_viewer_card.layout.display = ""    # görünür yap
        _render_mode(int(mode_select.value))

    def _on_mode_change(ch):
        _render_mode(int(ch["new"]))

    mode_select.observe(_on_mode_change, names="value")

    def _download_current_mode(_):
        mf = state.get("mode_files", {}) or {}
        m = int(mode_select.value)
        if m in mf:
            _download_file(mf[m])

    mode_dl.on_click(_download_current_mode)


    output_card = W.VBox(
        [output_title, output_info, table_box, status_box],
        layout=W.Layout(width="100%", gap="8px"),
    )

    output_card.add_class("hp-card")

    # wrap: dar ekranda alt alta düşsün; stretch: aynı hizada bitsin
    top_row = W.HBox(
        [form_card, viewer_card],
        layout=W.Layout(
            display="flex",
            flex_flow="row wrap",
            align_items="stretch",   # IMPORTANT: same bottom alignment
            justify_content="flex-start",
            gap="14px",
            width="100%",
        ),
    )
    web_page = W.VBox([top_row, output_card], layout=W.Layout(width="100%", gap="10px"))

    # Use RAW links so <img> works in HTML
    img_4cln = "https://raw.githubusercontent.com/enesemretas/hingeprot_fortran/main/assets/4cln.jpg"
    img_2bbm = "https://raw.githubusercontent.com/enesemretas/hingeprot_fortran/main/assets/2bbm.jpg"

    # Reference links (open in new tab)
    ref1_url = "https://pubs.acs.org/doi/10.1021/bi00188a001"
    ref2_url = "https://pubmed.ncbi.nlm.nih.gov/9218955/"
    ref3_url = "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.79.3090"
    ref4_url = "https://www.sciencedirect.com/science/article/pii/S000634950176033X"
    article_url = "https://pubmed.ncbi.nlm.nih.gov/17847101/"

    about_html = f"""
    <div class="hp-about">
      <h2>Motivation:</h2>

      <p>
        Proteins are highly flexible molecules. It is common to classify protein motions into shear and hinge motion
        <a class="hp-ref" href="{ref1_url}" target="_blank" rel="noopener noreferrer">[1]</a>.
        Shear motions are very limited and involve large number of residues. On the other hand, hinge motions are
        similar to rotations around an articulated joint and therefore can be very large. Hinge motion is characterized
        by large changes in main-chain torsional angles occurring at a localized region, which is called a hinge.
        Hinge motions usually involve a small number of residues, since even one bond can provide the required rotational
        freedom. This kind of protein motion is free of packing constraints. When a chain exhibits hinge motion at the
        region connecting two structural domains, each domain behaves as a rigid body and packing interactions can
        appear/disappear between the interfaces of those rigid bodies. Hinge motions usually occur upon binding to
        another molecule, or upon activation/deactivation of the protein.
      </p>

      <p>
        One of the most interesting examples is calmodulin. Upon binding to its ligands, there is large-scale movement
        of calmodulin involving splitting of one long helix. The total rotation of one domain relative to the other is
        upwards of 150 degrees (see images below and try the server for PDB codes 4cln and 2bbm, chainA).
      </p>

      <div class="hp-figrow">
        <div class="hp-figbox"><img src="{img_4cln}" alt="Calmodulin 4cln"></div>
        <div class="hp-figbox"><img src="{img_2bbm}" alt="Calmodulin 2bbm"></div>
      </div>

      <p>
        Therefore, hinge regions are the mechanistically informative regions of the structure and are of great importance
        in mediating cooperative motions that have functional importance.
      </p>

      <p>
        <i>HingeProt</i> is a web server for predicting rigid protein parts and the flexible hinge regions connecting them
        in the native topology of protein chains by employing elastic network (EN) models.
        <i>HingeProt</i> makes use of both Gaussian Network Model (GNM)
        <a class="hp-ref" href="{ref2_url}" target="_blank" rel="noopener noreferrer">[2]</a>,
        <a class="hp-ref" href="{ref3_url}" target="_blank" rel="noopener noreferrer">[3]</a>
        and Anisotropic Network models (ANM)
        <a class="hp-ref" href="{ref4_url}" target="_blank" rel="noopener noreferrer">[4]</a>.
      </p>

      <p>
        <i>HingeProt</i> server focuses on the prediction of the rigid parts and the hinge regions using a single static
        conformation of a protein structure. The hinge regions are the mechanistically informative regions of the
        structure and are of importance in mediating cooperative motions that have functional importance. GNM calculates
        the mean-square fluctuations and the correlation between the fluctuations of residues in the most dominant
        (slowest two) modes, which were shown to overlap with known protein motions. These suggest hinge regions and the
        cooperation between them. ANM provides the direction of the fluctuations in the corresponding modes.
      </p>

      <p>
        <i>HingeProt</i> is expected to be useful in a range of potential applications, especially in prediction
        protein-protein association by flexible docking and in refinement of the structure of the modeled complexes.
        <i>HingeProt</i> predictions are also helpful in fitting flexible hinge-bent protein structures into EM density
        maps and refining the EM structures. In addition, hinge regions can help in understanding functional mechanisms
        of macromolecular structures and assemblies.
      </p>

      <p>
        Given an input protein chain, <i>HingeProt</i> identifies the rigid parts and the hinges connecting them, and the
        direction of the fluctuation of each residue in the slowest two modes.
      </p>

      <h3>Method:</h3>

      <p>
        GNM <a class="hp-ref" href="{ref2_url}" target="_blank" rel="noopener noreferrer">[2]</a>,
        <a class="hp-ref" href="{ref3_url}" target="_blank" rel="noopener noreferrer">[3]</a> and its extension ANM
        <a class="hp-ref" href="{ref4_url}" target="_blank" rel="noopener noreferrer">[4]</a> are coarse-grained residue
        level elastic network models. GNM predicts the relative magnitudes of the fluctuations, whereas ANM predicts the
        directionalities of the collective motions in addition to their magnitudes. GNM results are more robust, and thus
        are preferentially used for evaluating square displacements in low frequency modes
        <a class="hp-ref" href="{ref4_url}" target="_blank" rel="noopener noreferrer">[4]</a>. Here GNM is used to
        calculate mean-square fluctuations and correlation between the fluctuations of residues, and ANM to generate the
        conformations that describe the fluctuations of residues from the average, X-ray, structure in the principal
        directions of motion.
      </p>

      <p>
        In <i>HingeProt</i>, first GNM decomposes the fluctuations of N residues of a structure into a series of N-1
        nonzero modes, given the Cartesian coordinates of Cα atoms. The eigenvectors corresponding to the slowest first
        and second modes are extracted. The square of these vectors describes the mean-square fluctuations (the
        autocorrelations) of residues from equilibrium positions along the principal coordinates (first and second modes
        here). Minima of mean square fluctuations at a given mode describe the flexible joints of the structure, i.e.
        the hinge regions, which connect the rigid units and mobile loops.
      </p>

      <p>
        GNM fluctuations being isotropic by definition, the directions of fluctuations are characterized by ANM. ANM
        predicts the fluctuations of N residues in the x, y and z directions from the average structure (X-ray or NMR) in
        3N-6 ANM nonzero modes. After mapping the ANM modes to GNM modes by comparing the square fluctuations, between
        the resulting modes in the two models, the directions of the fluctuations of residues in the slowest first and
        second modes of GNM are obtained by ANM analysis. As the fluctuations are symmetric with respect to the
        equilibrium positions, ANM predicted deformed structures could be obtained by adding and subtracting the
        fluctuations of each residue to/from its equilibrium position.
      </p>

      <p>
        For more detailed background of the work see the
        <a class="hp-ref" href="{article_url}" target="_blank" rel="noopener noreferrer">article</a>.
      </p>

      <div class="hp-refs">
        <h3 style="margin-top:10px;">References :</h3>
        <ol>
          <li id="ref1">
            <a class="hp-ref" href="{ref1_url}" target="_blank" rel="noopener noreferrer">[1]</a>
            Gerstein M, Lesk A. M. Lesk, Chothia C. (1994) Structural Mechanisms for Domain Movements in Proteins,
            Biochemistry 33(22), 6739-6749
          </li>
          <li id="ref2">
            <a class="hp-ref" href="{ref2_url}" target="_blank" rel="noopener noreferrer">[2]</a>
            Bahar, I., Atilgan A. R., Erman, B. (1997) Direct evaluation of thermal fluctuations in proteins using a
            single-parameter harmonic potential. Folding and Design, 2, 173-181.
          </li>
          <li id="ref3">
            <a class="hp-ref" href="{ref3_url}" target="_blank" rel="noopener noreferrer">[3]</a>
            Haliloglu, T., Bahar I, Erman B. (1997) Gaussian Dynamics of Proteins, Physical review letters, 79, 3090-3093
          </li>
          <li id="ref4">
            <a class="hp-ref" href="{ref4_url}" target="_blank" rel="noopener noreferrer">[4]</a>
            Atilgan, A. R., Durell, A. R., Jernigan, R. L., Demirel, M. C., Keskin, O., Bahar, I. (2001),
            Anisotropy of fluctuation dynamics of proteins with an elastic network model. Biophysical Journal, 80, 505-515
          </li>
        </ol>
      </div>
    </div>
    """


    about_page = W.HTML(about_html, layout=W.Layout(width="100%"))
   
    # ---------- References tab content (reads citations.txt next to ui.py) ----------
    base_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    CITATIONS_PATH = os.path.join(base_dir, "citations.txt")
    citations_text = _read_optional_file(CITATIONS_PATH)

    note_html = ""
    if not citations_text.strip():
        note_html = (
            '<div class="hp-note">'
            'Citations list not found. Create <b>citations.txt</b> next to <b>ui.py</b> and paste one citation per line.'
            "</div>"
        )

    items = []
    for line in (citations_text.splitlines() if citations_text else []):
        s = line.strip()
        if not s:
            continue
        s = re.sub(r"^[•\-\*]\s*", "", s)  # varsa madde işaretini temizle
        items.append(f"<li>{_linkify_citation_line(s)}</li>")

    ul_html = "<ul>" + "\n".join(items) + "</ul>" if items else "<ul></ul>"

    refs_html = f"""
    <div class="hp-refpage">
      <h2>References:</h2>

      <div class="sub">If you use this program, please cite the following:</div>
      <div class="cite-main">
        Emekli U, Schneidman-Duhovny D, Wolfson HJ, Nussinov R, Haliloglu T. (2008)
        HingeProt: Automated Prediction of Hinges in Protein Structures. Proteins, 70(4):1219-27.
      </div>

      <h3>Citations:</h3>
      {ul_html}
      {note_html}
    </div>
    """


    help_page = W.HTML("<div></div>")
    refs_page = W.HTML(refs_html, layout=W.Layout(width="100%"))

    main_view = W.VBox([web_page], layout=W.Layout(width="100%"))

    def _switch_page(_):
        key = nav.value
        if key == "web":
            main_view.children = [web_page]
        elif key == "about":
            main_view.children = [about_page]
        elif key == "help":
            main_view.children = [help_page]
        elif key == "refs":
            main_view.children = [refs_page]
        else:
            main_view.children = [web_page]

    nav.observe(_switch_page, names="value")
    _switch_page(None)

    display(css, header_bar, main_view)
    return None
