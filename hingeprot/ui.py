from __future__ import annotations

import os
import re
import datetime
import base64
import uuid
import json
import shutil
import subprocess
from pathlib import Path

import requests
import ipywidgets as W
from IPython.display import display


# Expose captured inputs to the notebook
LAST_UI_STATE: dict | None = None
LAST_INPUTS: dict | None = None


def get_last_inputs() -> dict | None:
    return LAST_INPUTS


# ----------------------------- shell helpers -----------------------------
def _sh(cmd: str, cwd: str | None = None, timeout: int | None = None) -> subprocess.CompletedProcess:
    """
    Run a bash command (Colab-friendly). Never raises; caller inspects returncode.
    """
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _which(tool: str) -> bool:
    return _sh(f"command -v {tool} >/dev/null 2>&1").returncode == 0


def _ldconfig_has_libg2c() -> bool:
    r = _sh(r'ldconfig -p | grep -E "libg2c\.so\.0" >/dev/null 2>&1')
    return r.returncode == 0


def _ensure_libg2c(log: callable) -> None:
    """
    Install libg2c.so.0 runtime (amd64) if missing.
    """
    if _ldconfig_has_libg2c():
        log("libg2c.so.0 already available (ldconfig).")
        return

    log("Installing libg2c.so.0 runtime (amd64) ...")
    os.makedirs("/content", exist_ok=True)
    os.chdir("/content")

    # download
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

    # install debs (allow fail) then fix deps
    r = _sh(f"dpkg -i {deb1} {deb2} || true")
    if r.stdout.strip():
        log(r.stdout.strip())
    if r.stderr.strip():
        log(r.stderr.strip())

    r = _sh("apt-get -y -qq -f install")
    if r.returncode != 0:
        raise RuntimeError(f"apt-get -f install failed:\n{r.stderr}")

    _sh("ldconfig")

    # verify
    if not _ldconfig_has_libg2c():
        raise RuntimeError("libg2c.so.0 still not found after installation.")
    log("✅ libg2c.so.0 installed and visible to ldconfig.")


def _ensure_repo(log: callable, fresh: bool = False) -> str:
    """
    Ensure /content/hingeprot_fortran exists and has hingeprot/ inside.
    Returns hingeprot directory path.
    """
    root = "/content/hingeprot_fortran"
    hp = os.path.join(root, "hingeprot")
    url = "https://github.com/enesemretas/hingeprot_fortran.git"

    # Safety: do not delete if we're currently running from inside that tree
    here = os.path.abspath(__file__)
    running_inside = here.startswith(os.path.abspath(root) + os.sep)

    if fresh:
        if running_inside:
            log("⚠️ Fresh clone requested but ui.py is running inside hingeprot_fortran; skipping rm -rf for safety.")
        else:
            log("Fresh clone: removing existing /content/hingeprot_fortran ...")
            shutil.rmtree(root, ignore_errors=True)

    if not os.path.isdir(hp):
        log("Cloning hingeprot_fortran repository ...")
        os.makedirs("/content", exist_ok=True)
        os.chdir("/content")
        r = _sh(f"git clone {url}")
        if r.returncode != 0:
            raise RuntimeError(f"git clone failed:\n{r.stderr}")
    else:
        log("Repository already present: /content/hingeprot_fortran")

    # sanity
    if not os.path.isdir(hp):
        raise RuntimeError("Repo clone seems incomplete: missing /content/hingeprot_fortran/hingeprot")

    # show pwd/ls as requested (into log)
    r = _sh("pwd", cwd=hp)
    log(f"pwd: {r.stdout.strip()}")
    r = _sh("ls -lah", cwd=hp)
    log("ls -lah:\n" + (r.stdout.strip() or ""))

    return hp


def _write_runHingeProt_pl(hingeprot_dir: str, gnm_cut: float, anm_cut: float, log: callable) -> str:
    """
    Create runHingeProt.pl with user-provided cutoffs.
    Note: part1 typically expects numeric; we pass rounded ints for stability.
    """
    gnm_i = int(round(float(gnm_cut)))
    anm_i = int(round(float(anm_cut)))
    if abs(gnm_i - float(gnm_cut)) > 1e-6 or abs(anm_i - float(anm_cut)) > 1e-6:
        log(f"NOTE: Cutoffs rounded to integers for part1: GNM {gnm_cut}→{gnm_i}, ANM {anm_cut}→{anm_i}")

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

    # perms
    _sh("chmod +x ./runHingeProt.pl || true", cwd=hingeprot_dir)
    _sh("find . -maxdepth 1 -type f -exec chmod +x {} \\; || true", cwd=hingeprot_dir)

    log(f"✅ runHingeProt.pl written at: {pl_path}")
    return pl_path


# ----------------------------- UI -----------------------------
def launch(runs_root: str = "/content/hingeprot_runs", save_json: bool = False):
    """
    Colab UI:
      - Collect PDB (code or upload), detect chains, collect cutoffs
      - Save Inputs -> LAST_INPUTS
      - Run HingeProt (Fortran) -> installs libs, ensures repo, writes runHingeProt.pl, runs perl
    """
    from google.colab import output  # colab-only
    output.enable_custom_widget_manager()

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

    # ---------- UI elements ----------
    css = W.HTML(r"""
    <style>
    .hp-card {border:1px solid #e5e7eb; border-radius:14px; padding:14px 16px; margin:10px 0; background:#fff;}
    .hp-banner{
      border:1px solid #e5e7eb;
      border-radius:16px;
      padding:14px 18px;
      margin:10px 0 12px 0;
      background:#fff;
      display:flex;
      align-items:center;
      gap:16px;
      box-shadow: 0 1px 0 rgba(0,0,0,0.03);
    }
    .hp-dot{ width:14px; height:14px; background:#ef4444; border-radius:999px; }
    .hp-title{ font-size:34px; font-weight:900; letter-spacing:0.5px; line-height: 1.0; margin:0;
               color:#111827; font-family: Arial, Helvetica, sans-serif; }
    .hp-title .prot{ color:#ef4444; }
    .hp-underline{ height:3px; width:280px; background:#111827; margin-top:6px; border-radius:999px; opacity:0.9; }
    .hp-tagline{ margin-top:6px; font-size:16px; font-weight:800; color:#dc2626; font-family: Arial, Helvetica, sans-serif; }
    .hp-pre{ white-space:pre-wrap; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
             font-size: 13px; line-height: 1.35; background:#0b1020; color:#e5e7eb; padding:12px; border-radius:12px; border:1px solid #1f2937;}
    </style>
    """)

    header = W.HTML(r"""
    <div class="hp-banner">
      <div class="hp-dot"></div>
      <div>
        <div class="hp-title">HINGE<span class="prot">prot</span></div>
        <div class="hp-underline"></div>
        <div class="hp-tagline">Colab UI + Fortran Runner</div>
      </div>
    </div>
    """)

    input_mode = W.ToggleButtons(
        options=[("Enter PDB code", "code"), ("Upload PDB file", "upload")],
        value="code",
        description="Input:",
        style={"description_width": "60px", "button_width": "170px"},
        layout=W.Layout(width="420px"),
    )

    pdb_code = W.Text(
        value="",
        description="PDB code:",
        placeholder="e.g., 3lzg",
        style={"description_width": "80px"},
        layout=W.Layout(width="420px"),
    )

    btn_choose_file = W.Button(description="Choose file", icon="upload", layout=W.Layout(width="180px"))
    file_lbl = W.Label("No file chosen")

    code_box = W.HBox([pdb_code], layout=W.Layout(align_items="center"))
    upload_box = W.HBox([btn_choose_file, file_lbl], layout=W.Layout(align_items="center", gap="10px"))

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
    btn_save = W.Button(description="Save Inputs", button_style="success", icon="save", layout=W.Layout(width="200px"))
    btn_run_fortran = W.Button(description="Run HingeProt (Fortran)", button_style="primary", icon="play", layout=W.Layout(width="280px"))
    btn_clear = W.Button(description="Clear", button_style="warning", icon="trash", layout=W.Layout(width="180px"))

    status_box = W.HTML('<div class="hp-pre">Load a PDB to detect chains.</div>')

    def _set_status(text: str):
        status_box.value = f'<div class="hp-pre">{_safe_html(text)}</div>'

    # ---------- state ----------
    state = {
        "pdb_text": None,
        "pdb_filename": None,   # e.g. 3lzg.pdb
        "pdb_path": None,
        "run_dir": None,
        "pdb_tag": None,
        "upload_name": None,
        "upload_bytes": None,
        "detected_chains": [],
        "chain_cbs": {},
        "manual_selection": (),
        "_syncing": False,
        "inputs": None,
        "hingeprot_dir": None,
        "last_out_dir": None,
    }
    global LAST_UI_STATE
    LAST_UI_STATE = state

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

    # ---------- uploader callback ----------
    cb_name = f"hingeprot_uploader_{uuid.uuid4().hex}"

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
            _set_status(f"Uploaded file: {name} ({len(data)} bytes)\nNow click 'Load / Detect Chains'.")
        except Exception as e:
            _set_status(f"Upload callback error: {e}")

    output.register_callback(cb_name, _js_upload_callback)

    def on_choose_file(_):
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
        output.eval_js(js)

    btn_choose_file.on_click(on_choose_file)

    # ---------- chain selection logic ----------
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

    def _on_chain_cb_change(_):
        _update_all_checkbox_from_selection()

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

    all_chains.observe(_on_all_chains_toggle, names="value")

    # ---------- logger ----------
    def log(msg: str):
        # append to status box
        current = re.sub(r"<[^>]*>", "", status_box.value)  # strip tags roughly
        if "Load a PDB" in current:
            current = ""
        text = (current.strip() + "\n" + msg).strip()
        _set_status(text[-8000:])  # cap

    # ---------- actions ----------
    def on_load_clicked(_):
        progress.value = 0
        progress.bar_style = "info"
        state["inputs"] = None
        state["last_out_dir"] = None

        try:
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = os.path.join(runs_root, f"run_{ts}")
            os.makedirs(run_dir, exist_ok=True)
            state["run_dir"] = run_dir

            # Determine PDB text + filename/tag
            if input_mode.value == "upload":
                if state["upload_bytes"] is None:
                    raise ValueError("Please click 'Choose file' and upload a PDB first.")
                pdb_text = state["upload_bytes"].decode("utf-8", errors="ignore")

                upname = (state.get("upload_name") or "upload.pdb").strip()
                # keep extension if present
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

            state["pdb_text"] = pdb_text
            state["pdb_filename"] = pdb_filename
            state["pdb_tag"] = tag

            # Save a copy to disk (runs_root)
            pdb_path = os.path.join(run_dir, pdb_filename)
            with open(pdb_path, "w", encoding="utf-8") as f:
                f.write(pdb_text)
            state["pdb_path"] = pdb_path

            # Detect chains
            chs = _detect_chains_from_text(pdb_text)
            if not chs:
                raise RuntimeError("No chains detected in the PDB.")
            state["detected_chains"] = chs

            default_sel = [chs[0]]
            state["manual_selection"] = tuple(default_sel)
            _rebuild_chain_checkboxes(chs, default_sel)

            state["_syncing"] = True
            try:
                all_chains.value = False
            finally:
                state["_syncing"] = False

            progress.value = 1
            _set_status(
                f"Loaded PDB (tag={tag})\n"
                f"Saved (runs_root): {pdb_path}\n"
                f"Detected chains: {', '.join(chs)}\n\n"
                "Now select chains, set cutoffs, then click 'Save Inputs' or directly 'Run HingeProt (Fortran)'."
            )

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

        captured = {
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
            "input_mode": input_mode.value,
            "pdb_code": pdb_code.value.strip() if input_mode.value == "code" else None,
            "upload_name": state.get("upload_name") if input_mode.value == "upload" else None,
            "pdb_tag": state.get("pdb_tag"),
            "pdb_filename": state.get("pdb_filename"),
            "pdb_path_runsroot": state.get("pdb_path"),
            "run_dir_runsroot": state.get("run_dir"),
            "detected_chains": list(detected),
            "all_chains": bool(all_chains.value),
            "selected_chains": list(chain_list),
            "chains_str": "".join(chain_list),   # IMPORTANT: no separators
            "gnm_cutoff_A": float(get_gnm_cut()),
            "anm_cutoff_A": float(get_anm_cut()),
        }
        return captured

    def on_save_clicked(_):
        try:
            captured = _capture_inputs()
            state["inputs"] = captured
            global LAST_INPUTS
            LAST_INPUTS = captured

            if save_json and captured.get("run_dir_runsroot"):
                outp = os.path.join(captured["run_dir_runsroot"], "inputs.json")
                with open(outp, "w", encoding="utf-8") as f:
                    json.dump(captured, f, ensure_ascii=False, indent=2)

            progress.value = max(progress.value, 2)
            progress.bar_style = "success"
            _set_status(
                "✅ Inputs captured (LAST_INPUTS)\n\n"
                f"pdb_filename   : {captured['pdb_filename']}\n"
                f"chains         : {captured['chains_str']}\n"
                f"GNM cutoff (Å) : {captured['gnm_cutoff_A']}\n"
                f"ANM cutoff (Å) : {captured['anm_cutoff_A']}\n\n"
                "Next: 'Run HingeProt (Fortran)'"
            )
        except Exception as e:
            progress.bar_style = "danger"
            _set_status(f"ERROR: {e}")

    def on_run_fortran_clicked(_):
        progress.bar_style = "info"
        try:
            captured = _capture_inputs()
            state["inputs"] = captured
            global LAST_INPUTS
            LAST_INPUTS = captured

            log("=== Step 1/4: Ensure libg2c.so.0 ===")
            progress.value = 1
            _ensure_libg2c(log)

            log("\n=== Step 2/4: Ensure repo ===")
            progress.value = 2
            hp_dir = _ensure_repo(log, fresh=False)
            state["hingeprot_dir"] = hp_dir

            log("\n=== Step 3/4: Write runHingeProt.pl (with UI cutoffs) ===")
            progress.value = 3
            _write_runHingeProt_pl(hp_dir, captured["gnm_cutoff_A"], captured["anm_cutoff_A"], log)

            # Place pdb file in hingeprot dir with expected filename
            pdb_filename = captured["pdb_filename"]
            pdb_abs = os.path.join(hp_dir, pdb_filename)
            with open(pdb_abs, "w", encoding="utf-8") as f:
                f.write(state["pdb_text"] or "")

            log(f"✅ PDB written to: {pdb_abs}")

            # Execute like:
            #   perl ./runHingeProt.pl 3lzg.pdb AB
            chains_str = captured["chains_str"]
            out_dir = os.path.join(hp_dir, f"{pdb_filename}.{chains_str}")

            # optional: remove previous output dir to avoid mixing
            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir, ignore_errors=True)

            log("\n=== Step 4/4: Run perl ./runHingeProt.pl ... ===")
            cmd = f"perl ./runHingeProt.pl {pdb_filename} {chains_str}"
            r = _sh(cmd, cwd=hp_dir)
            if r.stdout.strip():
                log("STDOUT:\n" + r.stdout.strip())
            if r.stderr.strip():
                log("STDERR:\n" + r.stderr.strip())
            if r.returncode != 0:
                raise RuntimeError(f"runHingeProt.pl failed (return code {r.returncode}).")

            state["last_out_dir"] = out_dir

            # Show directory listing (like your !ls -lah / !ls -lah outdir)
            r = _sh("ls -lah", cwd=hp_dir)
            log("\nls -lah (hingeprot dir):\n" + (r.stdout.strip() or ""))

            if os.path.isdir(out_dir):
                r2 = _sh(f"ls -lah {shlex_quote(os.path.basename(out_dir))}", cwd=hp_dir)
                log(f"\nls -lah {os.path.basename(out_dir)}:\n" + (r2.stdout.strip() or ""))
            else:
                log(f"\nWARNING: expected output dir not found: {out_dir}")

            progress.value = 4
            progress.bar_style = "success"
            log("\n✅ Done.")

        except Exception as e:
            progress.bar_style = "danger"
            log(f"\nERROR: {e}")

    def on_clear_clicked(_):
        pdb_code.value = ""
        input_mode.value = "code"
        state["upload_name"] = None
        state["upload_bytes"] = None
        file_lbl.value = "No file chosen"

        state["pdb_text"] = None
        state["pdb_filename"] = None
        state["pdb_path"] = None
        state["run_dir"] = None
        state["pdb_tag"] = None
        state["hingeprot_dir"] = None
        state["last_out_dir"] = None

        state["detected_chains"] = []
        state["chain_cbs"] = {}
        state["manual_selection"] = ()
        state["_syncing"] = False
        all_chains.value = False
        chains_wrap.children = ()

        state["inputs"] = None
        global LAST_INPUTS
        LAST_INPUTS = None

        progress.value = 0
        progress.max = 4
        progress.bar_style = ""
        _set_status("Cleared. Load a PDB to detect chains.")

    # small helper for safe ls on relative outdir
    def shlex_quote(s: str) -> str:
        return "'" + s.replace("'", "'\"'\"'") + "'"

    btn_load.on_click(on_load_clicked)
    btn_save.on_click(on_save_clicked)
    btn_run_fortran.on_click(on_run_fortran_clicked)
    btn_clear.on_click(on_clear_clicked)

    form_card = W.VBox([
        W.HTML('<div class="hp-card">'),
        input_mode,
        code_box,
        upload_box,
        btn_load,
        W.HTML("<hr>"),
        chain_row,
        W.VBox([gnm_row, anm_row], layout=W.Layout(gap="8px")),
        progress,
        W.HBox([btn_save, btn_run_fortran, btn_clear]),
        W.HTML("</div>"),
    ])

    output_card = W.VBox([
        W.HTML('<div class="hp-card"><b>Status / Logs</b></div>'),
        status_box,
    ])

    display(css, header, form_card, output_card)
    return state
