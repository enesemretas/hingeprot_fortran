cat > /content/hingeprot_fortran/hingeprot/ui.py <<'PY'
from __future__ import annotations

import os
import re
import datetime
import base64
import uuid
import json

import requests
import ipywidgets as W
from IPython.display import display

# Expose captured inputs to the notebook
LAST_UI_STATE: dict | None = None
LAST_INPUTS: dict | None = None


def get_last_inputs() -> dict | None:
    """Return the last captured inputs (after clicking 'Save Inputs')."""
    return LAST_INPUTS


def launch(runs_root: str = "/content/hingeprot_runs", save_json: bool = False):
    """
    INPUT-ONLY Colab UI:
      - Collects PDB (by code or upload)
      - Detects chains and lets user select
      - Collects GNM/ANM cutoffs
      - Saves everything into LAST_INPUTS (and optionally JSON)
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
        fbox = W.BoundedFloatText(value=default_value, min=minv, max=maxv, step=step, layout=W.Layout(width=value_width))
        value_box = W.Box([dropdown], layout=W.Layout(align_items="center"))

        def _on_toggle(ch):
            value_box.children = [dropdown] if ch["new"] == "list" else [fbox]

        toggle.observe(_on_toggle, names="value")

        def get_value() -> float:
            return float(dropdown.value) if toggle.value == "list" else float(fbox.value)

        row = W.HBox([lbl, toggle, value_box], layout=W.Layout(align_items="center", gap="12px"))
        return row, get_value

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

    gnm_row, get_gnm_cut = _list_or_custom_float(
        "GNM cutoff (Å):", options=[7, 8, 9, 10, 11, 12, 13, 20], default_value=10.0, minv=1.0, maxv=100.0
    )
    anm_row, get_anm_cut = _list_or_custom_float(
        "ANM cutoff (Å):", options=[10, 13, 15, 18, 20, 23, 36], default_value=18.0, minv=1.0, maxv=100.0
    )

    progress = W.IntProgress(value=0, min=0, max=2, description="Progress:", bar_style="")
    btn_save = W.Button(description="Save Inputs", button_style="success", icon="save", layout=W.Layout(width="320px"))
    btn_clear = W.Button(description="Clear", button_style="warning", icon="trash", layout=W.Layout(width="180px"))

    status_box = W.HTML('<div class="hp-pre">Load a PDB to detect chains.</div>')

    def _set_status(text: str):
        status_box.value = f'<div class="hp-pre">{_safe_html(text)}</div>'

    # ---------- state ----------
    state = {
        "pdb_text": None,
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

            # Determine pdb text
            if input_mode.value == "upload":
                if state["upload_bytes"] is None:
                    raise ValueError("Please click 'Choose file' and upload a PDB first.")
                pdb_text = state["upload_bytes"].decode("utf-8", errors="ignore")
                base = os.path.splitext(os.path.basename(state.get("upload_name") or "UPLOAD"))[0]
                tag = re.sub(r"[^0-9A-Za-z]+", "", base).upper() or "UPLOAD"
            else:
                code = pdb_code.value.strip()
                if not code:
                    raise ValueError("Please enter a PDB code (e.g., 4cln).")
                pdb_text = _fetch_pdb_text(code)
                tag = code.strip().upper()

            state["pdb_text"] = pdb_text
            state["pdb_tag"] = tag

            # Save a copy to disk for later steps
            pdb_path = os.path.join(run_dir, "input.pdb")
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
                f"Saved: {pdb_path}\n"
                f"Detected chains: {', '.join(chs)}\n\n"
                "Now select chains and click 'Save Inputs'."
            )

        except Exception as e:
            progress.bar_style = "danger"
            _set_status(f"ERROR: {e}")

    def on_save_clicked(_):
        progress.bar_style = "info"
        try:
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
                "input_mode": input_mode.value,  # "code" or "upload"
                "pdb_code": pdb_code.value.strip() if input_mode.value == "code" else None,
                "upload_name": state.get("upload_name") if input_mode.value == "upload" else None,
                "pdb_tag": state.get("pdb_tag"),
                "pdb_path": state.get("pdb_path"),
                "run_dir": state.get("run_dir"),
                "detected_chains": list(detected),
                "all_chains": bool(all_chains.value),
                "selected_chains": list(chain_list),
                "gnm_cutoff_A": float(get_gnm_cut()),
                "anm_cutoff_A": float(get_anm_cut()),
            }

            state["inputs"] = captured

            global LAST_INPUTS
            LAST_INPUTS = captured

            # optional JSON write
            if save_json and captured.get("run_dir"):
                outp = os.path.join(captured["run_dir"], "inputs.json")
                with open(outp, "w", encoding="utf-8") as f:
                    json.dump(captured, f, ensure_ascii=False, indent=2)

            progress.value = 2
            progress.bar_style = "success"

            _set_status(
                "✅ Inputs captured (saved into LAST_INPUTS)\n\n"
                f"pdb_tag        : {captured['pdb_tag']}\n"
                f"input_mode     : {captured['input_mode']}\n"
                f"pdb_path       : {captured['pdb_path']}\n"
                f"run_dir        : {captured['run_dir']}\n"
                f"chains         : {', '.join(captured['selected_chains'])}\n"
                f"GNM cutoff (Å) : {captured['gnm_cutoff_A']}\n"
                f"ANM cutoff (Å) : {captured['anm_cutoff_A']}\n\n"
                "Notebook'ta erişim:\n"
                "  from ui import get_last_inputs\n"
                "  get_last_inputs()"
            )

        except Exception as e:
            progress.bar_style = "danger"
            _set_status(f"ERROR: {e}")

    def on_clear_clicked(_):
        pdb_code.value = ""
        input_mode.value = "code"
        state["upload_name"] = None
        state["upload_bytes"] = None
        file_lbl.value = "No file chosen"

        state["pdb_text"] = None
        state["pdb_path"] = None
        state["run_dir"] = None
        state["pdb_tag"] = None

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
        progress.max = 2
        progress.bar_style = ""
        _set_status("Cleared. Load a PDB to detect chains.")

    btn_load.on_click(on_load_clicked)
    btn_save.on_click(on_save_clicked)
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
        W.HBox([btn_save, btn_clear]),
        W.HTML("</div>"),
    ])

    output_card = W.VBox([
        W.HTML('<div class="hp-card"><b>Status / Captured Inputs</b></div>'),
        status_box,
    ])

    display(css, header, form_card, output_card)
    return state
PY
