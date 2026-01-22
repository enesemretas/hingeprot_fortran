from __future__ import annotations

import os
import re
import json
import datetime
import base64
import uuid
import shutil
import subprocess

import requests
import ipywidgets as W
from IPython.display import display, clear_output, HTML


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

    here = os.path.abspath(__file__)
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


# ----------------------------- UI -----------------------------
def launch(runs_root: str = "/content/hingeprot_runs"):
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
      max-width: 1100px;
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
      padding:10px;
      border-radius:10px;
      border:1px solid #111827;
    }
    .hp-figbox img{
      display:block;
      max-width: 420px;
      width: 420px;
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
    upload_prog = W.IntProgress(value=0, min=0, max=100, description="", layout=W.Layout(width="160px"))
    upload_prog.bar_style = ""
    file_lbl = W.Label("No file chosen")

    code_box = W.HBox([pdb_code], layout=W.Layout(align_items="center"))
    upload_box = W.HBox([btn_choose_file, upload_prog, file_lbl], layout=W.Layout(align_items="center", gap="10px"))

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

    status_box = W.HTML('<div class="hp-pre">Load a PDB to detect chains.</div>')

    def _set_status(text: str):
        status_box.value = f'<div class="hp-pre">{_safe_html(text)}</div>'

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

        raw = v._make_html()
        raw = _html_with_unique_divid(raw)

        with viewer_out:
            clear_output(wait=True)
            display(HTML(raw))

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

    output.register_callback(cb_prog, _js_upload_progress_callback)
    output.register_callback(cb_name, _js_upload_callback)

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
        output.eval_js(js)

    btn_choose_file.on_click(on_choose_file)

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
            r = _sh(cmd, cwd=hp_dir)
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

            # ensure requested filename exists (alias)
            src_hinge = os.path.join(dest_out_dir, f"{pdb_filename}.hinge")
            alias_hinges = os.path.join(dest_out_dir, f"{pdb_filename}.new.hinges")
            if os.path.exists(src_hinge) and not os.path.exists(alias_hinges):
                try:
                    shutil.copyfile(src_hinge, alias_hinges)
                except Exception:
                    pass

            hinges_fp = _find_hinges_file(dest_out_dir, pdb_filename)
            if hinges_fp:
                txt = _read_text_file(hinges_fp, max_lines=900)
                _set_status(txt)
            else:
                _set_status(
                    "Run completed, but hinges file not found.\n"
                    f"Expected: {pdb_filename}.new.hinges (or fallback)\n"
                    f"Output folder: {dest_out_dir}"
                )

            progress.value = 4
            progress.bar_style = "success"

        except Exception as e:
            progress.bar_style = "danger"
            _set_status(f"ERROR: {e}")

    def on_clear_clicked(_):
        pdb_code.value = ""
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

        state["detected_chains"] = []
        state["chain_cbs"] = {}
        state["chain_colors"] = {}
        state["manual_selection"] = ()
        state["_syncing"] = False
        all_chains.value = False
        chains_wrap.children = ()

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

    output_title = W.HTML("<b>Hinges Output</b>")
    output_card = W.VBox([output_title, status_box], layout=W.Layout(width="100%", gap="8px"))
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


    about_page = W.HTML(about_html)

    help_page = W.HTML("<div></div>")
    refs_page = W.HTML("<div></div>")

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
