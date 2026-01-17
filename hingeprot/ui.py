# ui.py  (INPUT-ONLY version for Colab)
from __future__ import annotations

import os
import re
import datetime
import base64
import uuid

import requests
import ipywidgets as W
from IPython.display import display

# Expose captured inputs to the notebook
LAST_UI_STATE: dict | None = None
LAST_INPUTS: dict | None = None


def get_last_inputs() -> dict | None:
    """Return the last captured inputs (after clicking 'Save Inputs')."""
    return LAST_INPUTS


def launch(runs_root: str = "/content/hingeprot_runs"):
    """
    Launch the UI in the current notebook kernel (Colab).
    This INPUT-ONLY version:
      - collects PDB (code or upload),
      - detects chains,
      - collects user selections (chains + cutoffs),
      - stores them into LAST_INPUTS.
    """
    from google.colab import output  # colab-only
    output.enable_custom_widget_manager()

    os.makedirs(runs_root, exist_ok=True)

    # ---------- helpers ----------
    def _fetch_pdb(pdb_code: str) -> str:
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

    def _list_or_custom_float(
        label: str,
        options,
        default_value: float,
        minv: float,
        maxv: float,
        step: float = 0.1,
        label_width: str = "120px",
        toggle_width: str = "180px",
        value_width: str = "240px",
    ):
        opts = [float(x) for x in options]
        default_value = float(default_value)
        if default_value not in opts:
            opts = sorted(set(opts + [default_value]))
        else:
            opts = sorted(set(opts))

        lbl = W.Label(label, layout=W.Layout(width=label_width))
        toggle = W.ToggleButtons(
            options=[("List", "list"), ("Custom", "custom")],
            value="list",
            layout=W.Layout(width=toggle_width),
            style={"button_width": "80px"},
        )
        dropdown = W.Dropdown(options=opts, value=default_value, layout=W.Layout(width=value_width))
        fbox = W.BoundedFloatText(
            value=default_value, min=minv, max=maxv, step=step, layout=W.Layout(width=value_width)
        )
        value_box = W.Box([dropdown], layout=W.Layout(align_items="center"))

        def _on_toggle(ch):
            value_box.children = [dropdown] if ch["new"] == "list" else [fbox]

        toggle.observe(_on_toggle, names="value")

        def get_value() -> float:
            return float(dropdown.value) if toggle.value == "list" else float(fbox.value)

        row = W.HBox([lbl, toggle, value_box], layout=W.Layout(align_items="center", gap="12px"))
        return row, get_value

    def _colab_choose_file_js(cb_name: str):
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

    # ---------- UI ----------
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
    .hp-title{
      font-size:34px; font-weight:900; letter-spacing:0.5px; line-height: 1.0; margin:0;
      color:#111827; font-family: Arial, Helvetica, sans-serif;
    }
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
        <div class="hp-tagline">Input Collector (Colab UI)</div>
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
        placeholder="e.g., 4cln",
        style={"description_width": "80px"},
        layout=W.Layout(width="420px"),
    )

    btn_choose_file = W.Button(description="Choose file", icon="upload", layout=W.Layout(width="180px"))
    file_lbl = W.Label("No file chosen")

    code_box = W.HBox([pdb_code], layout=W.Layout(align_items="center"))
    upload_box = W.HBox([btn_choose_file, file_lbl], layout=W.Layout(align_items="center", gap="10px"))

    btn_load = W.Button(
        description="Load / Detect Chains",
        button_style="info",
        icon="search",
        layout=W.Layout(width="260px"),
    )

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

    gnm_row, get_gnm_cut = _list_or_custom_float(
        "GNM cutoff (Å):", options=[7, 8, 9, 10, 11, 12, 13, 20], default_value=10.0, minv=1.0, maxv=100.0
    )
    anm_row, get_anm_cut = _list_or_custom_float(
        "ANM cutoff (Å):", options=[10, 13, 15, 18, 20, 23, 36], default_value=18.0, minv=1.0, maxv=100.0
    )

    progress = W.IntProgress(value=0, min=0, max=2, description="Progress:", bar_style="")
    btn_save = W.Button(description="Save Inputs", button_style="success", icon="save", layout=W.Layout(width="320px"))
    btn_clear = W.Button(description="Clear", button_style="warning", icon="trash", layout=W.Layout(width="180px"))

    # output
    status_box = W.HTML('<div class="hp-pre">Load a PDB to detect chains.</div>')

    def _set_status(text: str):
        status_box.value = f'<div class="hp-pre">{_safe_html(text)}</div>'

    state = {
        # raw pdb
        "pdb_text": None,
        "pdb_path": None,          # optional saved path
        "run_dir": None,           # optional run dir
        "upload_name": None,
        "upload_bytes": None,

        # chain selection
        "detected_chains": [],
        "chain_cbs": {},
        "manual_selection": (),
        "_syncing": False,

        # tags
        "pdb_tag": None,

        # captured inputs
        "inputs": None,
    }

    # expose state globally for notebook access/debug
    global LAST_UI_STATE
    LAST_UI_STATE = state

    def _sync_input_visibility(*_):
        if input_mode.value == "code":
            code_box.layout.display = ""
            upload_box.layout.display = "none"
        else:
            code_box.layout.display = "none"
            upload_box.layout.display = ""

    _sync_input_visibility()
    input_mode.observe(lambda ch: _sync_input_visibility(), names="value")

    # uploader callback (unique)
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
            _set_status(f"Uploaded file received: {name} ({len(data)} bytes)\nNow click 'Load / Detect Chains'.")
        except Exception as e:
            _set_status(f"Upload callback error: {e}")

    output.register_callback(cb_name, _js_upload_callback)

    def on_choose_file(_):
        _colab_choose_file_js(cb_name)

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

    # ---------- actions ----------
    def on_load_clicked(_):
        progress.value = 0
        progress.bar_style = "info"
        state["inputs"] = None

        try:
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = os.path.join(runs_root, f"run_{ts}")
            os.makedirs(run_dir, exist_ok=True)
            state["run_dir"] = run_dir

            # Determine
