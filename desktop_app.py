"""Tkinter analysis workflows with shared calculations and durable result projects.

Run `python desktop_app.py`. Packaging and OS acceptance are documented in
docs/desktop_prototype.md; this entrypoint does not import Streamlit.
"""
from __future__ import annotations

import multiprocessing as mp
from pathlib import Path
import tempfile


def _analysis_worker(sender, output_path, data, settings):
    try:
        from mfrm_app.desktop_analysis import run_analysis
        from mfrm_app.project import save_project

        analysis = run_analysis(data, **settings, progress=lambda phase: sender.send(("progress", phase)))
        save_project(output_path, analysis)
        sender.send(("complete", ""))
    except Exception as exc:
        sender.send(("failed", f"{type(exc).__name__}: {exc}"))
    finally:
        sender.close()


class AnalysisWorker:
    """Single owned process; cancellation discards its unpublished result."""

    def __init__(self):
        self.process = None
        self.receiver = None
        self.directory = None
        self.terminal = None

    def start(self, data, settings):
        if self.process is not None:
            raise RuntimeError("An analysis is already active")
        self.directory = tempfile.TemporaryDirectory(prefix="mfrm-analysis-")
        self.output_path = Path(self.directory.name) / "result.mfrm"
        context = mp.get_context("spawn")
        self.receiver, sender = context.Pipe(duplex=False)
        self.process = context.Process(target=_analysis_worker, args=(sender, self.output_path, data, settings))
        try:
            self.process.start()
        except Exception:
            self.process = None
            self.close()
            raise
        finally:
            sender.close()

    def poll(self):
        if self.process is None:
            return []
        events = []
        while True:
            try:
                if not self.receiver.poll():
                    break
                status, detail = self.receiver.recv()
            except (EOFError, BrokenPipeError):
                # Windows can report pipe EOF from poll(), after the terminal
                # message was received but before the child has fully exited.
                # Keep that message and still verify exit status/result below.
                break
            if status in {"complete", "failed"}:
                self.terminal = (status, detail)
            else:
                events.append((status, detail))
        if not self.process.is_alive():
            self.process.join()
            terminal = self.terminal or ("failed", f"Worker exited without a result (exit={self.process.exitcode})")
            if terminal[0] == "complete" and self.process.exitcode == 0:
                try:
                    from mfrm_app.project import load_project
                    terminal = ("complete", load_project(self.output_path))
                except Exception as exc:
                    terminal = ("failed", str(exc))
            elif terminal[0] == "complete":
                terminal = ("failed", f"Worker exit={self.process.exitcode}")
            self.close()
            events.append(terminal)
        return events

    def cancel(self):
        """Stop and clean only this worker; never publish a partial project."""
        self.close()
        return ("cancelled", "Calculation cancelled; start a new analysis to run again.")

    def close(self):
        if self.process is not None:
            if self.process.is_alive():
                self.process.terminate()
            self.process.join()
            self.process.close()
        if self.receiver is not None:
            self.receiver.close()
        if self.directory is not None:
            self.directory.cleanup()
        self.process = self.receiver = self.directory = self.terminal = None


def demonstration_data():
    """Small artificial ratings for trying the workflow, not a validation study."""
    import pandas as pd
    rows = [(f"P{person:02}", f"R{rater}", f"T{task}", (person + 2 * rater + task) % 4)
            for person in range(1, 13) for rater in range(2) for task in range(2)]
    return pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Score"])


PALETTE = {
    "canvas_bg": "#F1F5F9",
    "card_bg": "#FFFFFF",
    "border": "#CBD5E1",
    "header_bg": "#0F172A",
    "header_fg": "#F8FAFC",
    "header_sub": "#94A3B8",
    "header_badge_bg": "#1E3A8A",
    "header_badge_fg": "#93C5FD",
    "nav_bg": "#FFFFFF",
    "primary": "#2563EB",
    "primary_hover": "#1D4ED8",
    "accent": "#0F766E",
    "accent_hover": "#115E59",
    "step_ready_bg": "#EFF6FF",
    "step_ready_fg": "#1D4ED8",
    "danger_bg": "#FEE2E2",
    "danger_fg": "#B91C1C",
    "text": "#0F172A",
    "muted": "#475569",
    "row_even": "#FFFFFF",
    "row_odd": "#F8FAFC",
    "table_head_bg": "#1E293B",
    "table_head_fg": "#F8FAFC",
    "banner": {
        "ok": ("#10B981", "#ECFDF5", "#065F46"),
        "info": ("#2563EB", "#EFF6FF", "#1E40AF"),
        "warn": ("#F59E0B", "#FFFBEB", "#92400E"),
        "error": ("#EF4444", "#FEF2F2", "#991B1B"),
    },
}


def _ui_font_family() -> str:
    import sys
    if sys.platform == "darwin":
        return "Hiragino Sans"
    if sys.platform == "win32":
        return "Yu Gothic UI"
    return "TkDefaultFont"


class DesktopApp:
    def __init__(self, root, language="en"):
        import tkinter as tk
        from tkinter import ttk
        from mfrm_app.engine import APP_VERSION

        self.language = language
        self.root, self.tk, self.ttk = root, tk, ttk
        self.font_family = _ui_font_family()
        self.worker = AnalysisWorker()
        self.analysis = self.data = None
        self.unsaved = self.busy = False
        self.widgets, self.canvases = [], []
        self.figure_cache = {}
        self.data_label = ""
        self.anchor_records, self.group_anchor_records = [], []
        self._configure_theme(root)
        root.title(self.tr(f"MFRM Desktop 試作 — {APP_VERSION}"))
        root.geometry("1120x800")
        root.minsize(860, 660)
        root.configure(bg=PALETTE["canvas_bg"])
        root.protocol("WM_DELETE_WINDOW", self.close)

        # Top brand bar
        header = tk.Frame(root, bg=PALETTE["header_bg"], padx=20, pady=10)
        header.pack(fill="x")
        brand_left = tk.Frame(header, bg=PALETTE["header_bg"])
        brand_left.pack(side="left")
        tk.Label(
            brand_left, text="MFRM-Desktop", bg=PALETTE["header_bg"], fg=PALETTE["header_fg"],
            font=(self.font_family, 15, "bold"),
        ).pack(side="left")
        self.header_sub_label = tk.Label(
            brand_left, text=self.tr("オフライン実データ分析ワークステーション"),
            bg=PALETTE["header_bg"], fg=PALETTE["header_sub"], font=(self.font_family, 11),
        )
        self.header_sub_label._i18n_static = True
        self.header_sub_label.pack(side="left", padx=(12, 0))
        tk.Label(
            header, text=f"v{APP_VERSION}", bg=PALETTE["header_badge_bg"], fg=PALETTE["header_badge_fg"],
            font=(self.font_family, 10, "bold"), padx=8, pady=2,
        ).pack(side="right", padx=(8, 0))
        ttk.Button(header, text=self.tr("ヘルプ"), style="Header.TButton", command=self.show_help).pack(side="right", padx=(6, 0))
        self.language_choice = tk.StringVar(value="English" if language == "en" else "日本語")
        language_box = ttk.Combobox(
            header, textvariable=self.language_choice, values=["English", "日本語"], state="readonly", width=10
        )
        language_box.pack(side="right", padx=6)
        language_box.bind("<<ComboboxSelected>>", self.change_language)
        self.widgets.append(language_box)

        # Step pill bar
        nav_outer = tk.Frame(root, bg=PALETTE["nav_bg"], highlightbackground=PALETTE["border"], highlightthickness=1)
        nav_outer.pack(fill="x")
        navigation = ttk.Frame(nav_outer, style="Nav.TFrame", padding=(18, 8))
        navigation.pack(fill="x")
        self.step_buttons = {}
        for key, label in (("data", self.tr("1. データ")), ("settings", self.tr("2. 設定")), ("results", self.tr("3. 結果"))):
            button = ttk.Button(navigation, text=label, command=lambda step=key: self.show_step(step))
            button.pack(side="left")
            self.step_buttons[key] = button
            if key != "results":
                ttk.Label(navigation, text="  ›  ", style="NavArrow.TLabel").pack(side="left")

        frame = ttk.Frame(root, style="Canvas.TFrame", padding=(18, 14))
        frame.pack(fill="both", expand=True)
        self.presentation = tk.StringVar(value="explained")
        self.gd_swap = tk.BooleanVar(value=False)
        self._syncing_gd_cols = False
        self.status = tk.StringVar()
        self.status_label = self.wrapped_label(frame, textvariable=self.status, style="Status.TLabel")
        self.status_label.pack(side="bottom", fill="x", pady=(6, 0))
        self.pages = {key: ttk.Frame(frame, style="Canvas.TFrame") for key in self.step_buttons}

        # Data: choose a source, convert wide sheets if needed, inspect the preview, then continue.
        page = self.pages["data"]
        toolbar_card = ttk.Frame(page, style="Card.TFrame", padding=(12, 10))
        toolbar_card.pack(fill="x", pady=(0, 10))
        for label, command, bstyle in (
            (self.tr("CSVを開く"), self.open_csv, "Primary.TButton"),
            (self.tr("サンプル集（8種）…"), self.open_sample_gallery, "Accent.TButton"),
            (self.tr("横長→縦長変換…"), self.open_wide_converter_manual, "Secondary.TButton"),
            (self.tr("サンプル"), self.use_demo, "Secondary.TButton"),
            (self.tr("G/D 2相サンプル"), self.use_gd_demo, "Secondary.TButton"),
            (self.tr("プロジェクトを開く"), self.open_project, "Secondary.TButton"),
        ):
            button = ttk.Button(toolbar_card, text=label, style=bstyle, command=command)
            button.pack(side="left", padx=(0, 6))
            self.widgets.append(button)
        ttk.Button(toolbar_card, text=self.tr("CSV形式"), style="Secondary.TButton", command=self.show_csv_guide).pack(side="right")

        self.data_banner_var = tk.StringVar(value=self.tr("「CSVを開く」または「サンプル集（8種）…」からデータを読み込んでください。"))
        self.data_banner_frame, self.input_label = self._make_banner(page, self.data_banner_var, kind="info")
        self.data_banner_frame.pack(fill="x", pady=(0, 10))

        bottom = ttk.Frame(page, style="Canvas.TFrame")
        bottom.pack(side="bottom", fill="x", pady=(10, 0))
        self.next_button = ttk.Button(
            bottom, text=self.tr("次へ"), style="Primary.TButton",
            command=lambda: self.show_step("settings"), state="disabled",
        )
        self.next_button.pack(side="right")
        self.preview = ttk.Frame(page, style="Card.TFrame", padding=10)
        self.preview.pack(fill="both", expand=True)

        # Settings: the scoring rubric owns the range; observed scores are context.
        page = self.pages["settings"]
        wf_bar = ttk.Frame(page, style="Card.TFrame", padding=(12, 8))
        wf_bar.pack(fill="x", pady=(0, 10))
        ttk.Label(wf_bar, text=self.tr("分析ワークフロー："), style="CardBold.TLabel").pack(side="left", padx=(0, 8))
        self.workflow = tk.StringVar(value="MFRM")
        choice = ttk.Combobox(wf_bar, textvariable=self.workflow, values=["MFRM", "G/D-Study"], state="readonly", width=22)
        choice.pack(side="left")
        choice.bind("<<ComboboxSelected>>", self.change_workflow)
        self.widgets.append(choice)
        ttk.Button(wf_bar, text=self.tr("設定の説明"), style="Secondary.TButton", command=self.show_settings_help).pack(side="right")

        self.settings_tabs = ttk.Notebook(page)
        self.settings_tabs.pack(fill="both", expand=True)
        self.settings_pages, self.settings_contents = {}, {}
        for key, title in (
            ("data", self.tr("列・採点基準")),
            ("model", self.tr("モデル・推定")),
            ("anchors", self.tr("アンカー・等化")),
            ("extras", self.tr("バイアス・区間")),
            ("gd", self.tr("G/D-Study計画")),
        ):
            panel = ttk.Frame(self.settings_tabs, style="Card.TFrame")
            canvas = tk.Canvas(panel, bg=PALETTE["card_bg"], highlightthickness=0, borderwidth=0)
            scrollbar = ttk.Scrollbar(panel, orient="vertical", command=canvas.yview)
            scrollbar.pack(side="right", fill="y")
            canvas.pack(fill="both", expand=True)
            canvas.configure(yscrollcommand=scrollbar.set)
            content = ttk.Frame(canvas, style="Card.TFrame", padding=14)
            window = canvas.create_window(0, 0, window=content, anchor="nw")
            content.bind("<Configure>", lambda e, c=canvas: c.configure(scrollregion=c.bbox("all")))
            canvas.bind("<Configure>", lambda e, c=canvas, w=window: c.itemconfigure(w, width=e.width))
            self.settings_tabs.add(panel, text=title)
            self.settings_pages[key], self.settings_contents[key] = panel, content

        page = self.settings_contents["data"]
        columns = ttk.LabelFrame(page, text=self.tr("列"), padding=12)
        columns.pack(fill="x")
        self.person, self.score = tk.StringVar(), tk.StringVar()
        for index, (label, variable) in enumerate(((self.tr("対象者"), self.person), (self.tr("得点"), self.score))):
            ttk.Label(columns, text=label).grid(row=0, column=index, sticky="w", padx=(0, 18))
            box = ttk.Combobox(columns, textvariable=variable, state="readonly", width=22)
            box.grid(row=1, column=index, sticky="w", padx=(0, 18), pady=(4, 14))
            self.widgets.append(box)
            if index == 0:
                self.person_box = box
            else:
                self.score_box = box
        ttk.Label(columns, text=self.tr("評価者・課題など（相）")).grid(row=2, column=0, sticky="w")
        facet_frame = ttk.Frame(columns)
        facet_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(4, 0))
        self.facets = tk.Listbox(
            facet_frame, selectmode="multiple", exportselection=False, height=5,
            bg="#FFFFFF", fg=PALETTE["text"], selectbackground="#DBEAFE", selectforeground="#1E3A8A",
            highlightthickness=1, highlightbackground=PALETTE["border"], relief="flat",
        )
        self.facets.pack(side="left", fill="x", expand=True)
        scroll = ttk.Scrollbar(facet_frame, orient="vertical", command=self.facets.yview)
        scroll.pack(side="right", fill="y")
        self.facets.configure(yscrollcommand=scroll.set)
        self.facets.bind("<<ListboxSelect>>", self.refresh_choices)
        self.widgets.append(self.facets)

        self.mfrm_preflight = tk.StringVar()
        self.mfrm_preflight_banner, _ = self._make_banner(page, self.mfrm_preflight, kind="info")
        self.mfrm_preflight_banner.pack(fill="x", pady=(12, 0))

        rubric = ttk.LabelFrame(page, text=self.tr("採点基準"), padding=14)
        rubric.pack(fill="x", pady=14)
        self.minimum, self.maximum = tk.StringVar(), tk.StringVar()
        self.range_entries = []
        for index, (label, variable) in enumerate(((self.tr("最低点"), self.minimum), (self.tr("最高点"), self.maximum))):
            ttk.Label(rubric, text=label).grid(row=0, column=index, sticky="w", padx=(0, 18))
            entry = ttk.Entry(rubric, textvariable=variable, width=10)
            entry.grid(row=1, column=index, sticky="w", padx=(0, 18), pady=(4, 0))
            self.widgets.append(entry)
            self.range_entries.append(entry)
        self.observed_scores = tk.StringVar()
        ttk.Label(rubric, textvariable=self.observed_scores, style="Muted.TLabel").grid(
            row=2, column=0, columnspan=2, sticky="w", pady=(8, 0)
        )
        self.score.trace_add("write", self.update_observed_scores)
        self.person.trace_add("write", self.update_facet_choices)
        self.score.trace_add("write", self.update_facet_choices)
        self.person.trace_add("write", self._sync_person_to_gd)
        self.score.trace_add("write", self._sync_score_to_gd)

        page = self.pages["settings"]
        self.input_error = tk.StringVar()
        self.error_banner, _ = self._make_banner(page, self.input_error, kind="error")
        self.input_error.trace_add("write", self._toggle_error_banner)
        actions = ttk.Frame(page, style="Canvas.TFrame")
        actions.pack(side="bottom", fill="x", pady=(10, 0))
        self.back_button = ttk.Button(actions, text=self.tr("戻る"), style="Secondary.TButton", command=lambda: self.show_step("data"))
        self.back_button.pack(side="left")
        self.run_button = ttk.Button(actions, text=self.tr("分析"), style="Primary.TButton", command=self.run, state="disabled")
        self.run_button.pack(side="right")
        self.cancel_button = ttk.Button(actions, text=self.tr("中止"), style="Danger.TButton", command=self.cancel, state="disabled")

        self.settings_tabs.pack_forget()
        self.settings_tabs.pack(fill="both", expand=True)
        self.build_extended_settings()

        # Results: only result actions and outputs remain on this page.
        page = self.pages["results"]
        actions = ttk.Frame(page, style="Card.TFrame", padding=(12, 8))
        actions.pack(fill="x", pady=(0, 10))
        self.result_label = ttk.Label(actions, text="", style="CardBold.TLabel")
        self.result_label.pack(side="left")
        self.project_button = ttk.Button(actions, text=self.tr("プロジェクトを保存"), style="Primary.TButton", command=self.save, state="disabled")
        self.project_button.pack(side="right", padx=(8, 0))
        self.report_button = ttk.Button(actions, text=self.tr("レポートを保存"), style="Accent.TButton", command=self.export_report, state="disabled")
        self.report_button.pack(side="right", padx=(8, 0))
        ttk.Button(actions, text=self.tr("分析情報"), style="Secondary.TButton", command=self.show_analysis_info).pack(side="right", padx=(8, 0))
        self.qualification = tk.StringVar()
        self.qualification_banner, self.qualification_label = self._make_banner(page, self.qualification, kind="warn")
        self.draft_notice = tk.StringVar()
        self.draft_banner, self.draft_label = self._make_banner(page, self.draft_notice, kind="info")
        self.notebook = ttk.Notebook(page)
        self.notebook.pack(fill="both", expand=True)
        self.show_step("data")
        if root.tk.call("tk", "windowingsystem") == "aqua":
            root.createcommand("::tk::mac::OpenDocument", self.open_documents)
        root.after(100, self.poll)

    def _configure_theme(self, root):
        style = self.ttk.Style(root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        ff = self.font_family
        style.configure(".", font=(ff, 11), background=PALETTE["card_bg"], foreground=PALETTE["text"])
        style.configure("Canvas.TFrame", background=PALETTE["canvas_bg"])
        style.configure("Card.TFrame", background=PALETTE["card_bg"], relief="solid", borderwidth=1, bordercolor=PALETTE["border"])
        style.configure("Nav.TFrame", background=PALETTE["nav_bg"])
        style.configure("NavArrow.TLabel", background=PALETTE["nav_bg"], foreground=PALETTE["header_sub"], font=(ff, 12, "bold"))
        style.configure("CardBold.TLabel", background=PALETTE["card_bg"], foreground=PALETTE["text"], font=(ff, 11, "bold"))
        style.configure("Muted.TLabel", background=PALETTE["card_bg"], foreground=PALETTE["muted"], font=(ff, 10))
        style.configure("Status.TLabel", background=PALETTE["canvas_bg"], foreground=PALETTE["muted"], font=(ff, 10))
        style.configure("TLabelframe", background=PALETTE["card_bg"], bordercolor=PALETTE["border"], relief="solid", borderwidth=1)
        style.configure("TLabelframe.Label", background=PALETTE["card_bg"], foreground="#1E3A8A", font=(ff, 11, "bold"))
        style.configure("TCheckbutton", background=PALETTE["card_bg"], foreground=PALETTE["text"])
        style.map("TCheckbutton", background=[("active", PALETTE["card_bg"])])
        style.configure("TRadiobutton", background=PALETTE["card_bg"], foreground=PALETTE["text"])
        style.map("TRadiobutton", background=[("active", PALETTE["card_bg"])])
        style.configure("TCombobox", fieldbackground="#FFFFFF", background="#F8FAFC", foreground=PALETTE["text"], arrowsize=13, padding=4)
        style.map("TCombobox", fieldbackground=[("readonly", "#FFFFFF"), ("disabled", "#F1F5F9")], foreground=[("disabled", "#94A3B8")])
        style.configure("TEntry", fieldbackground="#FFFFFF", foreground=PALETTE["text"], padding=4)

        # Button hierarchy
        style.configure("TButton", background="#FFFFFF", foreground=PALETTE["text"], bordercolor=PALETTE["border"], padding=(12, 6), font=(ff, 11))
        style.map("TButton", background=[("active", "#F1F5F9"), ("disabled", "#F8FAFC")], foreground=[("disabled", "#94A3B8")])
        style.configure("Secondary.TButton", background="#FFFFFF", foreground=PALETTE["text"], bordercolor=PALETTE["border"], padding=(12, 6), font=(ff, 11))
        style.map("Secondary.TButton", background=[("active", "#EFF6FF"), ("disabled", "#F8FAFC")], foreground=[("active", "#1D4ED8"), ("disabled", "#94A3B8")])
        style.configure("Primary.TButton", background=PALETTE["primary"], foreground="#FFFFFF", bordercolor=PALETTE["primary"], padding=(14, 6), font=(ff, 11, "bold"))
        style.map("Primary.TButton", background=[("active", PALETTE["primary_hover"]), ("disabled", "#CBD5E1")], foreground=[("disabled", "#64748B")])
        style.configure("Accent.TButton", background=PALETTE["accent"], foreground="#FFFFFF", bordercolor=PALETTE["accent"], padding=(13, 6), font=(ff, 11, "bold"))
        style.map("Accent.TButton", background=[("active", PALETTE["accent_hover"]), ("disabled", "#CBD5E1")], foreground=[("disabled", "#64748B")])
        style.configure("Danger.TButton", background=PALETTE["danger_bg"], foreground=PALETTE["danger_fg"], bordercolor="#FCA5A5", padding=(12, 6), font=(ff, 11, "bold"))
        style.map("Danger.TButton", background=[("active", "#FECACA"), ("disabled", "#F8FAFC")], foreground=[("disabled", "#94A3B8")])
        style.configure("Header.TButton", background="#1E293B", foreground="#F8FAFC", bordercolor="#334155", padding=(10, 4), font=(ff, 10))
        style.map("Header.TButton", background=[("active", "#334155")], foreground=[("active", "#FFFFFF")])
        style.configure("Current.TButton", background=PALETTE["primary"], foreground="#FFFFFF", bordercolor=PALETTE["primary"], padding=(14, 5), font=(ff, 11, "bold"))
        style.map("Current.TButton", background=[("active", PALETTE["primary_hover"])], foreground=[("active", "#FFFFFF")])
        style.configure("StepReady.TButton", background=PALETTE["step_ready_bg"], foreground=PALETTE["step_ready_fg"], bordercolor="#BFDBFE", padding=(14, 5), font=(ff, 11))
        style.map("StepReady.TButton", background=[("active", "#DBEAFE")])

        # Notebook tabs
        style.configure("TNotebook", background=PALETTE["canvas_bg"], borderwidth=0, tabmargins=(2, 4, 2, 0))
        style.configure("TNotebook.Tab", background="#E2E8F0", foreground="#475569", padding=(14, 7), font=(ff, 11))
        style.map(
            "TNotebook.Tab",
            background=[("selected", "#FFFFFF"), ("active", "#EFF6FF"), ("disabled", "#F1F5F9")],
            foreground=[("selected", "#1D4ED8"), ("active", "#1E3A8A"), ("disabled", "#94A3B8")],
        )

        # Treeview tables
        style.configure("Treeview", background="#FFFFFF", fieldbackground="#FFFFFF", foreground=PALETTE["text"], rowheight=26, bordercolor=PALETTE["border"], font=(ff, 10))
        style.configure("Treeview.Heading", background=PALETTE["table_head_bg"], foreground=PALETTE["table_head_fg"], padding=(8, 6), font=(ff, 10, "bold"), relief="flat")
        style.map("Treeview.Heading", background=[("active", "#334155")])
        style.map("Treeview", background=[("selected", "#DBEAFE")], foreground=[("selected", "#1E3A8A")])

    def _make_banner(self, parent, textvariable, kind="info"):
        bar_color, bg_color, fg_color = PALETTE["banner"][kind]
        outer = self.tk.Frame(parent, bg=bg_color, highlightbackground=bar_color, highlightthickness=1)
        stripe = self.tk.Frame(outer, bg=bar_color, width=5)
        stripe.pack(side="left", fill="y")
        lbl = self.tk.Label(
            outer, textvariable=textvariable, bg=bg_color, fg=fg_color,
            justify="left", anchor="w", padx=10, pady=7, font=(self.font_family, 10),
        )
        lbl.pack(side="left", fill="x", expand=True)
        lbl.bind("<Configure>", lambda e, w=lbl: w.configure(wraplength=max(120, e.width - 24)))
        outer._stripe, outer._label = stripe, lbl
        return outer, lbl

    def _set_banner_kind(self, banner_frame, kind):
        bar_color, bg_color, fg_color = PALETTE["banner"][kind]
        banner_frame.configure(bg=bg_color, highlightbackground=bar_color)
        banner_frame._stripe.configure(bg=bar_color)
        banner_frame._label.configure(bg=bg_color, fg=fg_color)

    def _toggle_error_banner(self, *args):
        if self.input_error.get().strip():
            self.error_banner.pack(fill="x", before=self.settings_tabs, pady=(0, 8))
        else:
            self.error_banner.pack_forget()

    def tr(self, message):
        from mfrm_app.desktop_i18n import text
        return text(message, getattr(self, 'language', 'ja'))

    def view_settings(self):
        from mfrm_app.desktop_i18n import text
        out = dict(figure=text(self.map_choice.get(), 'ja'), facet=('すべて' if self.subset_box.current() == 0 else self.map_subset.get()),
                   statistic=text(self.statistic.get(), 'ja'), presentation=self.presentation.get(), language=self.language)
        if self.gd_swap.get():
            out['gd_swap'] = '1'
        return out

    def change_language(self, event=None):
        if self.busy:
            return
        view = self.view_settings() if self.analysis is not None else None
        self.language = 'ja' if self.language_choice.get() == '日本語' else 'en'

        def visit(widget):
            if widget.winfo_class() in ('TLabel', 'TButton', 'TCheckbutton', 'TRadiobutton', 'TLabelframe') or getattr(widget, '_i18n_static', False):
                widget.configure(text=self.tr(str(widget.cget('text'))))
            if widget.winfo_class() == 'TNotebook':
                for tab in widget.tabs():
                    widget.tab(tab, text=self.tr(widget.tab(tab, 'text')))
            for child in widget.winfo_children():
                visit(child)

        visit(self.root)
        self.root.title(self.tr(self.root.title()))
        # Translate only fixed application choices, never user column/level names.
        for variable in (self.interval_mode, self.joint_mode, self.gd_fixed_mode, self.gd_plan_mode):
            variable.set(self.tr(variable.get()))
        fixed_vars = {str(self.interval_mode), str(self.joint_mode), str(self.gd_fixed_mode), str(self.gd_plan_mode)}
        for panel in self.settings_contents.values():
            def choices(widget):
                if widget.winfo_class() == 'TCombobox' and str(widget.cget('textvariable')) in fixed_vars:
                    widget.configure(values=[self.tr(v) for v in widget.cget('values')])
                for child in widget.winfo_children():
                    choices(child)
            choices(panel)
        self.status.set(self.tr(self.status.get()))
        self.input_error.set(self.tr(self.input_error.get()))
        self.data_banner_var.set(self.tr(self.data_banner_var.get()))
        self.update_observed_scores()
        self.update_gd_labels()
        self.update_mfrm_preflight()
        self.refresh_anchor_lists(clear_error=False)
        for widget in (self.bias_list, self.interval_list):
            selected = widget.curselection()
            items = [self.tr(widget.get(i)) for i in range(widget.size())]
            widget.delete(0, 'end')
            for item in items:
                widget.insert('end', item)
            for index in selected:
                widget.selection_set(index)
        if self.analysis is not None:
            page, unsaved = self.step, self.unsaved
            self.show_analysis(self.analysis, view_override=view)
            self.unsaved = unsaved
            self.show_step(page)

    def show_figure_caption(self):
        from mfrm_app.figure_presentation import caption_text
        return self.guide_window(self.tr('タイトル・図注'),
            caption_text(self.current_figure._mfrm_caption, self.analysis['analysis_id']))

    def wrapped_label(self, parent, **options):
        # Let the allocated width govern wrapping, including after a language switch.
        label = self.ttk.Label(parent, width=1, wraplength=1, **options)
        label.bind("<Configure>", lambda e: label.configure(wraplength=max(1, e.width)))
        return label

    def _styled_listbox(self, parent, *, height=4, selectmode="browse"):
        return self.tk.Listbox(
            parent, height=height, selectmode=selectmode, exportselection=False,
            bg="#FFFFFF", fg=PALETTE["text"], selectbackground="#DBEAFE", selectforeground="#1E3A8A",
            highlightthickness=1, highlightbackground=PALETTE["border"], relief="flat",
        )

    def option_field(self, frame, row, label, value="", choices=None, column=0):
        variable = self.tk.StringVar(value=value)
        self.ttk.Label(frame, text=label, wraplength=150).grid(row=row, column=column, sticky="w", padx=(0, 8), pady=5)
        widget = (self.ttk.Combobox(frame, textvariable=variable, values=choices, state="readonly", width=16)
                  if choices is not None else self.ttk.Entry(frame, textvariable=variable, width=16))
        frame.columnconfigure(column+1, weight=1)
        widget.grid(row=row, column=column+1, sticky="ew", pady=5)
        self.widgets.append(widget)
        return variable, widget

    def build_extended_settings(self):
        panel = self.settings_contents["model"]
        self.model_defaults = dict(model="RSM", method="JMLE", step_facet="", quad_points="15",
                                   population_prior_sd="1", maxit="400", reltol="1e-6")
        self.model_vars, self.model_boxes = {}, {}
        for row, (key, label, choices) in enumerate([
            ("model", self.tr("モデル"), ["RSM", "PCM"]), ("method", self.tr("推定法"), ["JMLE", "MML"]),
            ("step_facet", self.tr("PCMの段階相"), []), ("quad_points", self.tr("MMLの求積点数"), None),
            ("population_prior_sd", self.tr("MMLの固定母集団SD"), None), ("maxit", self.tr("最大反復数"), None), ("reltol", self.tr("停止許容誤差"), None)]):
            v, box = self.option_field(panel, row, label, self.model_defaults[key], choices)
            self.model_vars[key], self.model_boxes[key] = v, box
            if choices is not None:
                box.bind("<<ComboboxSelected>>", self.refresh_choices)
        self.wrapped_label(panel, style="Muted.TLabel", text=self.tr("RSM：共通の段階。PCM：選んだ相の水準別の段階。\nMMLは固定SD・EMです。指定区間は対象ごとに数値条件を確認します。")).grid(row=7, column=0, columnspan=2, sticky="ew", pady=10)

        # Anchors / Scale Equating tab
        anchor_panel = self.settings_contents["anchors"]
        self.mfrm_preflight_anchor_banner, _ = self._make_banner(anchor_panel, self.mfrm_preflight, kind="info")
        self.mfrm_preflight_anchor_banner.pack(fill="x", pady=(0, 10))

        fixed_box = self.ttk.LabelFrame(anchor_panel, text=self.tr("固定アンカー（水準のロジット値を固定：Status A）"), padding=10)
        fixed_box.pack(fill="x")
        self.anchor_facet, self.anchor_facet_box = self.option_field(fixed_box, 0, self.tr("相"), choices=[])
        self.anchor_level, self.anchor_level_box = self.option_field(fixed_box, 0, self.tr("水準"), choices=[], column=2)
        self.anchor_value, _ = self.option_field(fixed_box, 1, self.tr("アンカー値"), "0.0")
        self.anchor_facet_box.bind("<<ComboboxSelected>>", self.refresh_choices)
        self.anchor_list = self._styled_listbox(fixed_box, height=4)
        self.anchor_list.grid(row=2, column=0, columnspan=3, sticky="ew", pady=(6, 0), padx=(0, 8))
        self.anchor_list.bind("<<ListboxSelect>>", self._on_select_anchor)
        a_btns = self.ttk.Frame(fixed_box)
        a_btns.grid(row=1, column=3, rowspan=2, sticky="ns")
        for label, command in [
            (self.tr("追加"), self.add_anchor),
            (self.tr("選択を削除"), self.remove_anchor),
            (self.tr("CSVから一括読込…"), self.import_anchor_csv),
            (self.tr("全解除"), self.clear_anchors),
        ]:
            btn = self.ttk.Button(a_btns, text=label, style="Secondary.TButton", command=command)
            btn.pack(fill="x", pady=2)
            self.widgets.append(btn)

        group_box = self.ttk.LabelFrame(anchor_panel, text=self.tr("グループアンカー（指定水準群の平均ロジットを固定：Status G）"), padding=10)
        group_box.pack(fill="x", pady=10)
        self.g_anchor_facet, self.g_anchor_facet_box = self.option_field(group_box, 0, self.tr("相"), choices=[])
        self.g_anchor_group, _ = self.option_field(group_box, 0, self.tr("グループ名"), "G1", column=2)
        self.g_anchor_mean, _ = self.option_field(group_box, 1, self.tr("グループ平均"), "0.0")
        self.g_anchor_facet_box.bind("<<ComboboxSelected>>", self.refresh_choices)
        self._last_g_anchor_facet = None
        self.ttk.Label(group_box, text=self.tr("グループに含める水準（2つ以上選択）"), wraplength=150).grid(row=2, column=0, sticky="nw", pady=4)
        g_lvl_frame = self.ttk.Frame(group_box)
        g_lvl_frame.grid(row=2, column=1, columnspan=2, sticky="ew", pady=4, padx=(0, 8))
        self.g_anchor_levels_list = self._styled_listbox(g_lvl_frame, height=4, selectmode="multiple")
        self.g_anchor_levels_list.pack(side="left", fill="x", expand=True)
        g_scroll = self.ttk.Scrollbar(g_lvl_frame, orient="vertical", command=self.g_anchor_levels_list.yview)
        g_scroll.pack(side="right", fill="y")
        self.g_anchor_levels_list.configure(yscrollcommand=g_scroll.set)
        self.widgets.append(self.g_anchor_levels_list)
        self.g_anchor_list = self._styled_listbox(group_box, height=3)
        self.g_anchor_list.grid(row=3, column=0, columnspan=3, sticky="ew", pady=(6, 0), padx=(0, 8))
        g_btns = self.ttk.Frame(group_box)
        g_btns.grid(row=2, column=3, rowspan=2, sticky="ns")
        for label, command in [
            (self.tr("グループを追加・更新"), self.add_group_anchor),
            (self.tr("選択を削除"), self.remove_group_anchor),
            (self.tr("全解除"), self.clear_group_anchors),
        ]:
            btn = self.ttk.Button(g_btns, text=label, style="Secondary.TButton", command=command)
            btn.pack(fill="x", pady=2)
            self.widgets.append(btn)
        self.wrapped_label(
            anchor_panel, style="Muted.TLabel",
            text=self.tr("アンカー指定時は主効果・適合度・条件付きバイアスを等化尺度上で推定します（固定SD MMLの指定区間・共同交互作用は総和ゼロ制約専用のためアンカー未使用時のみ有効）。"),
        ).pack(fill="x", pady=(4, 0))

        # Extras: Bias, Specified Intervals, Joint Interaction
        panel = self.settings_contents["extras"]
        self.bias_requests, self.interval_requests = [], []
        bias = self.ttk.LabelFrame(panel, text=self.tr("バイアス（一覧に登録した条件で計算）"), padding=10)
        bias.pack(fill="x")
        self.bias_a, self.bias_a_box = self.option_field(bias, 0, self.tr("相A"), choices=[])
        self.bias_b, self.bias_b_box = self.option_field(bias, 0, self.tr("相B"), choices=[], column=2)
        self.bias_bound, _ = self.option_field(bias, 1, self.tr("探索範囲の絶対値"), "10")
        self.bias_omit = self.tk.BooleanVar(value=True)
        check = self.ttk.Checkbutton(bias, text=self.tr("極端な水準を除外"), variable=self.bias_omit)
        check.grid(row=1, column=2, columnspan=2, sticky="w"); self.widgets.append(check)
        self.bias_list = self._styled_listbox(bias, height=3)
        self.bias_list.grid(row=2, column=0, columnspan=3, sticky="ew", padx=(0, 8), pady=(4, 0))
        self.bias_list.bind("<<ListboxSelect>>", lambda e: self.edit_request("bias"))
        buttons = self.ttk.Frame(bias); buttons.grid(row=2, column=3, rowspan=2, sticky="ns")
        for label, command in [(self.tr("追加"), self.add_bias), (self.tr("選択を更新"), lambda: self.add_bias(update=True)), (self.tr("選択を削除"), lambda: self.remove_request("bias"))]:
            button = self.ttk.Button(buttons, text=label, style="Secondary.TButton", command=command)
            button.pack(fill="x", pady=2); self.widgets.append(button)
        self.wrapped_label(bias, text=self.tr("有意性・バイアスの有無の結論は保留。"), foreground="#92400E").grid(row=3, column=0, columnspan=3, sticky="ew", pady=(4, 0))
        box = self.ttk.LabelFrame(panel, text=self.tr("指定区間：固定SD MML（一覧に登録した条件で計算）"), padding=10)
        box.pack(fill="x", pady=10)
        self.interval_facet, self.interval_facet_box = self.option_field(box, 0, self.tr("相"), choices=[])
        self.interval_level, _ = self.option_field(box, 0, self.tr("信頼水準"), ".95", column=2)
        self.interval_mode, mode_box = self.option_field(box, 1, self.tr("対象"), self.tr("全水準"), [self.tr("全水準"), self.tr("A − B（事前指定）")])
        self.interval_a, self.interval_a_box = self.option_field(box, 2, self.tr("水準A"), choices=[])
        self.interval_b, self.interval_b_box = self.option_field(box, 2, self.tr("水準B"), choices=[], column=2)
        self.interval_facet_box.bind("<<ComboboxSelected>>", self.refresh_choices)
        self.interval_list = self._styled_listbox(box, height=3)
        self.interval_list.grid(row=3, column=0, columnspan=3, sticky="ew", padx=(0, 8), pady=(4, 0))
        self.interval_list.bind("<<ListboxSelect>>", lambda e: self.edit_request("interval"))
        buttons = self.ttk.Frame(box); buttons.grid(row=3, column=3, rowspan=2, sticky="ns")
        for label, command in [(self.tr("追加"), self.add_interval), (self.tr("選択を更新"), lambda: self.add_interval(update=True)), (self.tr("選択を削除"), lambda: self.remove_request("interval"))]:
            button = self.ttk.Button(buttons, text=label, style="Secondary.TButton", command=command)
            button.pack(fill="x", pady=2); self.widgets.append(button)
        self.wrapped_label(box, text=self.tr("保留対象は区間を描かず理由を保持します。"), foreground="#92400E").grid(row=4, column=0, columnspan=3, sticky="ew", pady=(4, 0))
        box = self.ttk.LabelFrame(panel, text=self.tr("共同交互作用：主効果・段階も同時に再推定する別モデル"), padding=10)
        box.pack(fill="x", pady=10)
        self.joint_enabled = self.tk.BooleanVar(value=False)
        toggle = self.ttk.Checkbutton(box, text=self.tr("次回の分析に追加（MML・固定SD 1のみ）"), variable=self.joint_enabled)
        toggle.grid(row=0, column=0, columnspan=4, sticky="w"); self.widgets.append(toggle)
        self.joint_a, self.joint_a_box = self.option_field(box, 1, self.tr("相A"), choices=[])
        self.joint_b, self.joint_b_box = self.option_field(box, 1, self.tr("相B"), choices=[], column=2)
        self.joint_level, _ = self.option_field(box, 2, self.tr("信頼水準"), ".95")
        self.joint_mode, _ = self.option_field(box, 2, self.tr("対象"), self.tr("全セル"), [self.tr("全セル"), self.tr("全セル＋指定差")], column=2)
        self.joint_cell_vars, self.joint_cell_boxes = {}, {}
        for i, key in enumerate(("a0", "b0", "a1", "b1")):
            label = (self.tr("差の左側") if i < 2 else self.tr("差の右側")) + (self.tr("：A水準") if key[0] == 'a' else self.tr("：B水準"))
            var, widget = self.option_field(box, 3+i//2, label, choices=[], column=2*(i%2))
            self.joint_cell_vars[key], self.joint_cell_boxes[key] = var, widget
        for widget in (self.joint_a_box, self.joint_b_box):
            widget.bind("<<ComboboxSelected>>", self.refresh_choices)
        self.wrapped_label(box, style="Muted.TLabel", text=self.tr("数値条件を満たすセル・指定差に点ごとの区間と近似p値を表示します。\nHolm補正の対象は、この分析の全セルと指定差です。区間は同時区間ではありません。\n元の主効果モデルと固定較正の診断結果は保持します。公平性の判定にはなりません。")).grid(row=5, column=0, columnspan=4, sticky="ew", pady=8)

        # G/D-Study tab
        panel = self.settings_contents["gd"]
        self.gd_enabled = self.tk.BooleanVar(value=False)
        check = self.ttk.Checkbutton(panel, text=self.tr("MFRMと併せてG/D-Studyも計算（G/D-Study入口では常に実行）"), variable=self.gd_enabled)
        check.grid(row=0, column=0, columnspan=4, sticky="w"); self.widgets.append(check)
        self.gd_object, self.gd_object_box = self.option_field(panel, 1, self.tr("測定対象の列"), choices=[])
        self.gd_score, self.gd_score_box = self.option_field(panel, 1, self.tr("観測得点の列"), choices=[], column=2)
        self.gd_object_box.bind("<<ComboboxSelected>>", lambda e: self.update_gd_facet_choices())
        self.gd_score_box.bind("<<ComboboxSelected>>", lambda e: self.update_gd_facet_choices())
        self.ttk.Label(panel, text=self.tr("ランダム相（1〜2列）"), wraplength=150).grid(row=2, column=0, sticky="nw", pady=4)
        gd_list_frame = self.ttk.Frame(panel)
        gd_list_frame.grid(row=2, column=1, columnspan=3, sticky="ew", pady=4)
        self.gd_facets = self._styled_listbox(gd_list_frame, selectmode="multiple", height=4)
        self.gd_facets.pack(side="left", fill="x", expand=True)
        self.widgets.append(self.gd_facets)
        self.gd_facets.bind("<<ListboxSelect>>", self.update_gd_labels)
        self.gd_swap_order_button = self.ttk.Button(gd_list_frame, text=self.tr("順序入替"), style="Secondary.TButton", command=self.swap_gd_random_facets, state="disabled")
        self.gd_swap_order_button.pack(side="left", padx=(8, 0))
        self.widgets.append(self.gd_swap_order_button)
        self.gd_facet_count_var = self.tk.StringVar()
        self.ttk.Label(panel, textvariable=self.gd_facet_count_var, style="Muted.TLabel").grid(row=3, column=1, columnspan=3, sticky="w", pady=(0, 4))

        gd_fixed_choices = [
            self.tr("指定しない（選択した相のみで均衡交差）"),
            self.tr("1水準を抽出"),
            self.tr("全水準の合計得点（完全均衡を確認）"),
            self.tr("全水準の平均得点（完全均衡を確認）"),
        ]
        self.gd_fixed_mode, self.gd_fixed_mode_box = self.option_field(
            panel, 4, self.tr("固定相の扱い（3相以上のデータ）"), gd_fixed_choices[0], gd_fixed_choices
        )
        self.gd_fixed_mode_box.grid(row=4, column=1, columnspan=3, sticky="ew", pady=4)
        self.gd_fixed_mode_box.bind("<<ComboboxSelected>>", lambda e: self.update_gd_facet_choices())
        self.gd_fixed_facet, self.gd_fixed_facet_box = self.option_field(panel, 5, self.tr("固定相の列"), choices=[])
        self.gd_fixed_level, self.gd_fixed_level_box = self.option_field(panel, 5, self.tr("抽出する水準"), choices=[], column=2)
        self.gd_fixed_facet_box.bind("<<ComboboxSelected>>", lambda e: self.update_gd_facet_choices())
        self.gd_fixed_level_box.bind("<<ComboboxSelected>>", self.update_gd_labels)

        self.gd_preflight = self.tk.StringVar()
        self.gd_preflight_banner, self.gd_preflight_label = self._make_banner(panel, self.gd_preflight, kind="ok")
        self.gd_preflight_banner.grid(row=6, column=0, columnspan=4, sticky="ew", pady=(6, 10))

        plan_choices = [self.tr("観測水準数を使う"), self.tr("計画数を指定")]
        self.gd_plan_mode, self.gd_plan_mode_box = self.option_field(
            panel, 7, self.tr("D-Study計画の指定"), plan_choices[0], plan_choices
        )
        self.gd_plan_mode_box.bind("<<ComboboxSelected>>", self._on_gd_plan_mode_change)
        self.gd_limit, _ = self.option_field(panel, 7, self.tr("計画表の上限"), "10", column=2)
        self.gd_count_labels, self.gd_counts, self.gd_count_entries = [], [], []
        for i in range(2):
            label = self.ttk.Label(panel, text=self.tr(f"計画の相{i+1}の数"), wraplength=150)
            label.grid(row=8+i, column=0, sticky="w")
            v = self.tk.StringVar()
            entry = self.ttk.Entry(panel, textvariable=v, width=16)
            entry.grid(row=8+i, column=1, sticky="ew", pady=4)
            v.trace_add("write", lambda *a: self._sync_gd_plan_mode_from_entries())
            self.gd_count_labels.append(label); self.gd_counts.append(v); self.gd_count_entries.append(entry); self.widgets.append(entry)
        self.gd_population, _ = self.option_field(panel, 10, self.tr("対象母集団（宣言）"))
        self.wrapped_label(panel, style="Muted.TLabel", text=self.tr("元の観測得点を使用します。MFRM推定・採点範囲の指定は不要です。\n完全均衡交差、各セル1観測、単位重み。3相以上のデータは「固定相の扱い」で1水準抽出または完全均衡を確認した合計／平均集約を選べます。\nG／Phiは点推定です。区間は保留です。")).grid(row=11, column=0, columnspan=4, sticky="ew", pady=10)
        self.refresh_choices()

    def _sync_person_to_gd(self, *args):
        if self._syncing_gd_cols or not hasattr(self, "gd_object"): return
        self._syncing_gd_cols = True
        try:
            if self.person.get() and (not self.gd_object.get() or not self.gd_enabled.get()):
                self.gd_object.set(self.person.get())
                self.update_gd_facet_choices()
        finally:
            self._syncing_gd_cols = False

    def _sync_score_to_gd(self, *args):
        if self._syncing_gd_cols or not hasattr(self, "gd_score"): return
        self._syncing_gd_cols = True
        try:
            if self.score.get() and (not self.gd_score.get() or not self.gd_enabled.get()):
                self.gd_score.set(self.score.get())
                self.update_gd_facet_choices()
        finally:
            self._syncing_gd_cols = False

    def _gd_fixed_mode_code(self):
        from mfrm_app.desktop_i18n import text
        ja_val = text(self.gd_fixed_mode.get(), 'ja')
        return {
            "1水準を抽出": "level",
            "全水準の合計得点（完全均衡を確認）": "sum",
            "全水準の平均得点（完全均衡を確認）": "mean",
        }.get(ja_val, "none")

    def _set_gd_fixed_mode_code(self, code):
        ja_label = {
            "level": "1水準を抽出",
            "sum": "全水準の合計得点（完全均衡を確認）",
            "mean": "全水準の平均得点（完全均衡を確認）",
        }.get(code, "指定しない（選択した相のみで均衡交差）")
        self.gd_fixed_mode.set(self.tr(ja_label))

    def _gd_fixed_condition_dict(self):
        mode = self._gd_fixed_mode_code()
        if mode == "none":
            return None
        facet = self.gd_fixed_facet.get().strip()
        level = self.gd_fixed_level.get().strip() if mode == "level" else None
        if not facet or (mode == "level" and not level):
            raise ValueError(self.tr("固定相の列と水準を選択してください。"))
        return dict(facet=facet, mode=mode, level=level)

    def _on_gd_plan_mode_change(self, event=None):
        from mfrm_app.desktop_i18n import text
        if text(self.gd_plan_mode.get(), 'ja') == "観測水準数を使う":
            for v in self.gd_counts:
                v.set("")
        self.update_gd_labels()

    def _sync_gd_plan_mode_from_entries(self):
        from mfrm_app.desktop_i18n import text
        if any(v.get().strip() for v in self.gd_counts):
            if text(self.gd_plan_mode.get(), 'ja') != "計画数を指定":
                self.gd_plan_mode.set(self.tr("計画数を指定"))
                self._refresh_gd_entry_states()

    def _refresh_gd_entry_states(self):
        from mfrm_app.desktop_i18n import text
        if self.busy:
            return
        selected = [self.gd_facets.get(i) for i in self.gd_facets.curselection()]
        custom = text(self.gd_plan_mode.get(), 'ja') == "計画数を指定"
        for i, entry in enumerate(self.gd_count_entries):
            entry.configure(state="normal" if (custom and i < len(selected)) else "disabled")
        mode_code = self._gd_fixed_mode_code()
        self.gd_fixed_facet_box.configure(state="readonly" if mode_code != "none" else "disabled")
        self.gd_fixed_level_box.configure(state="readonly" if mode_code == "level" else "disabled")
        self.gd_swap_order_button.configure(state="normal" if len(selected) == 2 else "disabled")

    def swap_gd_random_facets(self):
        selected = list(self.gd_facets.curselection())
        if len(selected) != 2:
            return
        items = [self.gd_facets.get(i) for i in range(self.gd_facets.size())]
        i0, i1 = selected[0], selected[1]
        items[i0], items[i1] = items[i1], items[i0]
        self.gd_facets.delete(0, "end")
        for idx, col in enumerate(items):
            self.gd_facets.insert("end", col)
            if idx in (i0, i1):
                self.gd_facets.selection_set(idx)
        c0, c1 = self.gd_counts[0].get(), self.gd_counts[1].get()
        self.gd_counts[0].set(c1); self.gd_counts[1].set(c0)
        self.update_gd_labels()

    def update_gd_facet_choices(self, *, reset=False):
        if self.data is None or not hasattr(self, "gd_facets"):
            return
        obj = self.gd_object.get()
        scr = self.gd_score.get()
        fixed_candidates = [c for c in self.data.columns if c not in {obj, scr, "Weight"}]
        self.gd_fixed_facet_box.configure(values=fixed_candidates)
        mode_code = self._gd_fixed_mode_code()
        prev_selected = [] if reset else [self.gd_facets.get(i) for i in self.gd_facets.curselection()]

        if reset:
            if 1 <= len(fixed_candidates) <= 2:
                tentative_chosen = list(fixed_candidates)
            elif len(fixed_candidates) >= 3:
                preferred = [c for c in ("Rater", "Task", "Item", "評価者", "課題") if c in fixed_candidates]
                tentative_chosen = (preferred + [c for c in fixed_candidates if c not in preferred])[:2]
            else:
                tentative_chosen = []
        else:
            tentative_chosen = [c for c in prev_selected if c in fixed_candidates]

        if reset or self.gd_fixed_facet.get() not in fixed_candidates:
            unselected = [c for c in fixed_candidates if c not in tentative_chosen]
            pool = unselected or fixed_candidates
            default_fixed = next(
                (c for c in ("Criterion", "観点") if c in pool),
                pool[0] if pool else "",
            )
            self.gd_fixed_facet.set(default_fixed)

        ff = self.gd_fixed_facet.get()
        levels = list(dict.fromkeys(self.data[ff].dropna().astype(str))) if ff in self.data.columns else []
        levels = [lv for lv in levels if lv.strip()]
        self.gd_fixed_level_box.configure(values=levels)
        if reset or self.gd_fixed_level.get() not in levels:
            self.gd_fixed_level.set(levels[0] if levels else "")

        excluded = {obj, scr, "Weight"}
        if mode_code != "none" and ff:
            excluded.add(ff)
        candidates = [c for c in self.data.columns if c not in excluded]
        chosen = [c for c in tentative_chosen if c in candidates]
        ordered = [*chosen, *[c for c in candidates if c not in chosen]]
        self.gd_facets.delete(0, "end")
        for idx, col in enumerate(ordered):
            self.gd_facets.insert("end", col)
            if col in chosen:
                self.gd_facets.selection_set(idx)
        self.update_gd_labels()

    def change_workflow(self, event=None):
        gd = self.workflow.get() in {"G/D-Study", "G/D Study"}
        if gd and self.data is not None:
            if not self.gd_object.get() and self.person.get():
                self.gd_object.set(self.person.get())
            if not self.gd_score.get() and self.score.get():
                self.gd_score.set(self.score.get())
            self.update_gd_facet_choices()
        for key in ("data", "model", "anchors", "extras"):
            self.settings_tabs.tab(self.settings_pages[key], state="disabled" if gd else "normal")
        self.settings_tabs.select(self.settings_pages["gd" if gd else "data"])

    def _facet_levels_for_name(self, facet_name):
        if self.data is None or not facet_name:
            return []
        col = self.person.get() if facet_name == "Person" else facet_name
        if col not in self.data.columns:
            return []
        vals = [str(v).strip() for v in self.data[col].dropna().tolist()]
        return list(dict.fromkeys(v for v in vals if v))

    def refresh_choices(self, event=None):
        if not hasattr(self, "model_boxes"):
            return
        facets = [self.facets.get(i) for i in self.facets.curselection()]
        all_facets = ["Person", *facets]
        self.model_boxes["step_facet"].configure(values=facets)
        self.model_boxes["step_facet"].configure(state="readonly" if self.model_vars["model"].get() == "PCM" else "disabled")
        for key in ("quad_points", "population_prior_sd"):
            self.model_boxes[key].configure(state="normal" if self.model_vars["method"].get() == "MML" else "disabled")
        for box in (self.bias_a_box, self.bias_b_box):
            box.configure(values=all_facets)
        self.interval_facet_box.configure(values=facets)
        facet = self.interval_facet.get()
        levels = self._facet_levels_for_name(facet)
        for box in (self.interval_a_box, self.interval_b_box):
            box.configure(values=levels)
        for box in (self.joint_a_box, self.joint_b_box):
            box.configure(values=facets)
        for key, box in self.joint_cell_boxes.items():
            name = self.joint_a.get() if key[0] == 'a' else self.joint_b.get()
            box.configure(values=self._facet_levels_for_name(name))

        # Anchor controls
        self.anchor_facet_box.configure(values=all_facets)
        if self.anchor_facet.get() not in all_facets:
            self.anchor_facet.set(facets[0] if facets else "Person")
        a_levels = self._facet_levels_for_name(self.anchor_facet.get())
        self.anchor_level_box.configure(values=a_levels)
        if self.anchor_level.get() not in a_levels:
            self.anchor_level.set(a_levels[0] if a_levels else "")

        self.g_anchor_facet_box.configure(values=all_facets)
        if self.g_anchor_facet.get() not in all_facets:
            self.g_anchor_facet.set(facets[0] if facets else "Person")
        gf = self.g_anchor_facet.get()
        g_levels = self._facet_levels_for_name(gf)
        prev_g_sel = (
            {self.g_anchor_levels_list.get(i) for i in self.g_anchor_levels_list.curselection()}
            if self._last_g_anchor_facet == gf else set()
        )
        self._last_g_anchor_facet = gf
        self.g_anchor_levels_list.delete(0, "end")
        for idx, lv in enumerate(g_levels):
            self.g_anchor_levels_list.insert("end", lv)
            if lv in prev_g_sel:
                self.g_anchor_levels_list.selection_set(idx)

        self.update_mfrm_preflight()
        self.update_gd_labels()

    def update_mfrm_preflight(self):
        if not hasattr(self, "mfrm_preflight"):
            return
        if self.data is None:
            self.mfrm_preflight.set("")
            return
        from mfrm_app.desktop_analysis import preview_mfrm_connectivity
        facets = [self.facets.get(i) for i in self.facets.curselection()]
        info = preview_mfrm_connectivity(
            self.data, person_col=self.person.get(), score_col=self.score.get(),
            facet_cols=facets, anchors=self.anchor_records, group_anchors=self.group_anchor_records,
        )
        msg = info["summary_ja"] if self.language == "ja" else info["summary_en"]
        self.mfrm_preflight.set(msg)
        kind_map = {
            "ok": "ok",
            "anchored_ok": "info",
            "disconnected": "warn",
            "incomplete": "warn",
            "anchor_error": "error",
        }
        bkind = kind_map.get(info["status"], "info")
        self._set_banner_kind(self.mfrm_preflight_banner, bkind)
        if hasattr(self, "mfrm_preflight_anchor_banner"):
            self._set_banner_kind(self.mfrm_preflight_anchor_banner, bkind)

    def update_gd_labels(self, event=None):
        if not hasattr(self, "gd_facets"):
            return
        facets = [self.gd_facets.get(i) for i in self.gd_facets.curselection()]
        for i, label in enumerate(self.gd_count_labels):
            label.configure(text=self.tr(f"計画の {facets[i]} の数") if i < len(facets) else self.tr(f"計画の相{i+1}の数"))
        if len(facets) == 1:
            self.gd_facet_count_var.set(self.tr(f"選択中：1列（横軸: {facets[0]}）"))
        elif len(facets) == 2:
            self.gd_facet_count_var.set(self.tr(f"選択中：2列（横軸: {facets[0]} / 系列: {facets[1]}）"))
        else:
            self.gd_facet_count_var.set(self.tr(f"選択中：{len(facets)}列（1〜2列を選択）"))
        self._refresh_gd_entry_states()
        if self.data is None:
            self.gd_preflight.set("")
            return
        from mfrm_app.desktop_analysis import preview_gd_design
        fc = None
        if self._gd_fixed_mode_code() != "none":
            mode = self._gd_fixed_mode_code()
            ff = self.gd_fixed_facet.get().strip()
            fl = self.gd_fixed_level.get().strip() if mode == "level" else None
            fc = dict(facet=ff, mode=mode, level=fl)
        info = preview_gd_design(
            self.data, object_facet=self.gd_object.get(), score_col=self.gd_score.get(),
            random_facets=facets, fixed_condition=fc,
        )
        msg = info["summary_ja"] if self.language == "ja" else info["summary_en"]
        self.gd_preflight.set(msg)
        self._set_banner_kind(self.gd_preflight_banner, "ok" if info["ok"] else "warn")

    def _on_select_anchor(self, event=None):
        sel = self.anchor_list.curselection()
        if not sel or sel[0] >= len(self.anchor_records):
            return
        rec = self.anchor_records[sel[0]]
        self.anchor_facet.set(rec["Facet"])
        self.refresh_choices()
        self.anchor_level.set(rec["Level"])
        self.anchor_value.set(f"{rec['Anchor']:g}")

    def add_anchor(self):
        from mfrm_app.basic_analysis import normalize_anchor_records, normalize_group_anchor_records
        try:
            facet = self.anchor_facet.get().strip()
            level = self.anchor_level.get().strip()
            if not facet or not level:
                raise ValueError(self.tr("入力データにある水準を選択してください。"))
            val = self.number(self.anchor_value, self.tr("アンカー値"))
            facets = [self.facets.get(i) for i in self.facets.curselection()]
            updated = [r for r in self.anchor_records if not (r["Facet"] == facet and r["Level"] == level)]
            updated.append({"Facet": facet, "Level": level, "Anchor": val})
            norm_fixed = normalize_anchor_records(updated, person_col=self.person.get(), facet_cols=facets)
            if self.group_anchor_records:
                normalize_group_anchor_records(
                    self.group_anchor_records, person_col=self.person.get(), facet_cols=facets, fixed_records=norm_fixed,
                )
            self.anchor_records = norm_fixed
            self.refresh_anchor_lists()
        except Exception as exc:
            self.input_error.set(self.tr(str(exc)))

    def remove_anchor(self):
        for i in reversed(self.anchor_list.curselection()):
            if i < len(self.anchor_records):
                self.anchor_records.pop(i)
        self.refresh_anchor_lists()

    def clear_anchors(self):
        self.anchor_records = []
        self.refresh_anchor_lists()

    def import_anchor_csv(self, path=None):
        from tkinter import filedialog
        import pandas as pd
        from mfrm_app.desktop_guidance import read_desktop_csv
        from mfrm_app.basic_analysis import normalize_anchor_records, normalize_group_anchor_records
        if path is None:
            path = filedialog.askopenfilename(
                parent=self.root, title=self.tr("アンカーCSV（Facet, Level, Anchor または Estimate 列）を選択"),
                filetypes=[("CSV", "*.csv")],
            )
        if not path:
            return
        try:
            df = read_desktop_csv(path)
            col_map = {c.strip().lower(): c for c in df.columns}
            f_col = col_map.get("facet") or col_map.get("相")
            l_col = col_map.get("level") or col_map.get("element") or col_map.get("水準")
            a_col = col_map.get("anchor") or col_map.get("estimate") or col_map.get("measure") or col_map.get("アンカー値")
            s_col = col_map.get("status")
            if not f_col or not l_col or not a_col:
                raise ValueError(self.tr("アンカーCSVには Facet, Level, Anchor（または Estimate）列が必要です。"))
            imported = []
            for _, row in df.iterrows():
                if s_col:
                    st = str(row[s_col]).strip().upper()
                    if st and st not in {"A", "ANCHOR", "FIXED"}:
                        continue
                raw_val = pd.to_numeric(str(row[a_col]).strip(), errors="coerce")
                if pd.isna(raw_val):
                    continue
                f_val = str(row[f_col]).strip()
                l_val = str(row[l_col]).strip()
                if f_val and l_val:
                    imported.append({"Facet": f_val, "Level": l_val, "Anchor": float(raw_val)})
            if not imported:
                raise ValueError(self.tr("アンカーCSVには Facet, Level, Anchor（または Estimate）列が必要です。"))
            facets = [self.facets.get(i) for i in self.facets.curselection()]
            norm_imp = normalize_anchor_records(imported, person_col=self.person.get(), facet_cols=facets)
            imp_keys = {(r["Facet"], r["Level"]) for r in norm_imp}
            merged = [r for r in self.anchor_records if (r["Facet"], r["Level"]) not in imp_keys] + norm_imp
            norm_merged = normalize_anchor_records(merged, person_col=self.person.get(), facet_cols=facets)
            if self.group_anchor_records:
                normalize_group_anchor_records(
                    self.group_anchor_records, person_col=self.person.get(), facet_cols=facets, fixed_records=norm_merged,
                )
            self.anchor_records = norm_merged
            self.refresh_anchor_lists()
            self.status.set(self.tr(f"アンカーを{len(norm_imp)}件読み込みました。"))
        except Exception as exc:
            self.input_error.set(self.tr(str(exc)))

    def _grouped_anchor_Display_items(self):
        groups = []
        seen = {}
        for r in self.group_anchor_records:
            key = (r["Facet"], r["Group"])
            if key not in seen:
                seen[key] = {"Facet": r["Facet"], "Group": r["Group"], "GroupValue": r["GroupValue"], "Levels": []}
                groups.append(seen[key])
            seen[key]["Levels"].append(r["Level"])
        return groups

    def add_group_anchor(self):
        from mfrm_app.basic_analysis import normalize_group_anchor_records
        try:
            facet = self.g_anchor_facet.get().strip()
            group = self.g_anchor_group.get().strip() or "G1"
            mean_val = self.number(self.g_anchor_mean, self.tr("グループ平均"))
            levels = [self.g_anchor_levels_list.get(i) for i in self.g_anchor_levels_list.curselection()]
            if len(levels) < 2:
                raise ValueError(self.tr("グループアンカーには2つ以上の水準を選択してください。"))
            facets = [self.facets.get(i) for i in self.facets.curselection()]
            kept = [r for r in self.group_anchor_records if not (r["Facet"] == facet and r["Group"] == group)]
            added = [{"Facet": facet, "Level": lv, "Group": group, "GroupValue": mean_val} for lv in levels]
            norm_grp = normalize_group_anchor_records(
                kept + added, person_col=self.person.get(), facet_cols=facets, fixed_records=self.anchor_records,
            )
            self.group_anchor_records = norm_grp
            self.refresh_anchor_lists()
        except Exception as exc:
            self.input_error.set(self.tr(str(exc)))

    def remove_group_anchor(self):
        groups = self._grouped_anchor_Display_items()
        sel = self.g_anchor_list.curselection()
        if not sel or sel[0] >= len(groups):
            return
        target = (groups[sel[0]]["Facet"], groups[sel[0]]["Group"])
        self.group_anchor_records = [
            r for r in self.group_anchor_records if (r["Facet"], r["Group"]) != target
        ]
        self.refresh_anchor_lists()

    def clear_group_anchors(self):
        self.group_anchor_records = []
        self.refresh_anchor_lists()

    def refresh_anchor_lists(self, *, clear_error=True):
        if not hasattr(self, "anchor_list"):
            return
        self.anchor_list.delete(0, "end")
        for r in self.anchor_records:
            self.anchor_list.insert("end", self.tr(f"{r['Facet']} : {r['Level']}  =  {r['Anchor']:+.4g} logits"))
        self.g_anchor_list.delete(0, "end")
        for g in self._grouped_anchor_Display_items():
            lvl_str = ", ".join(g["Levels"])
            self.g_anchor_list.insert(
                "end", self.tr(f"{g['Facet']} [Group={g['Group']}, 平均={g['GroupValue']:+.4g}] : {lvl_str}")
            )
        if clear_error:
            self.input_error.set("")
        self.update_mfrm_preflight()

    def add_bias(self, update=False):
        from mfrm_app.desktop_analysis import _bias_requests
        try:
            request = dict(facet_a=self.bias_a.get(), facet_b=self.bias_b.get(), max_abs=float(self.bias_bound.get()), omit_extreme=self.bias_omit.get())
            facets = ["Person", *[self.facets.get(i) for i in self.facets.curselection()]]
            self.bias_requests = _bias_requests(self.replace_request("bias", request, update), facets)
            self.refresh_request_lists()
        except Exception as exc:
            self.input_error.set(self.tr(str(exc)))

    def add_interval(self, update=False):
        from mfrm_app.desktop_analysis import _interval_requests
        try:
            request = dict(facet=self.interval_facet.get(), level=float(self.interval_level.get()),
                difference=[self.interval_a.get(), self.interval_b.get()] if self.interval_mode.get() != self.tr("全水準") else None)
            requests = _interval_requests(self.replace_request("interval", request, update), [self.facets.get(i) for i in self.facets.curselection()], self.model_vars["method"].get())
            if request['difference'] and any(x not in self.interval_a_box['values'] for x in request['difference']):
                raise ValueError(self.tr("入力データにある水準を選択してください。"))
            self.interval_requests = requests
            self.refresh_request_lists()
        except Exception as exc:
            self.input_error.set(self.tr(str(exc)))

    def replace_request(self, kind, request, update):
        requests = list(getattr(self, kind + "_requests"))
        selected = getattr(self, kind + "_list").curselection()
        if update:
            if not selected:
                raise ValueError(self.tr("更新する条件を一覧から選択してください。"))
            requests[selected[0]] = request
        else:
            requests.append(request)
        return requests

    def edit_request(self, kind):
        selected = getattr(self, kind + "_list").curselection()
        if not selected:
            return
        request = getattr(self, kind + "_requests")[selected[0]]
        if kind == "bias":
            self.bias_a.set(request["facet_a"]); self.bias_b.set(request["facet_b"])
            self.bias_bound.set(str(request["max_abs"])); self.bias_omit.set(request["omit_extreme"])
        else:
            self.interval_facet.set(request["facet"]); self.interval_level.set(str(request["level"]))
            self.interval_mode.set(self.tr("A − B（事前指定）") if request["difference"] else self.tr("全水準"))
            self.interval_a.set(request["difference"][0] if request["difference"] else "")
            self.interval_b.set(request["difference"][1] if request["difference"] else "")
        self.refresh_choices()

    def remove_request(self, kind):
        widget, requests = getattr(self, kind + "_list"), getattr(self, kind + "_requests")
        for i in reversed(widget.curselection()):
            requests.pop(i)
        self.refresh_request_lists()

    def refresh_request_lists(self):
        for kind in ("bias", "interval"):
            widget = getattr(self, kind + "_list"); widget.delete(0, "end")
            for r in getattr(self, kind + "_requests"):
                label = (self.tr(f"{r['facet_a']} × {r['facet_b']} · 範囲 ±{r['max_abs']:g} · 極端水準{'除外' if r['omit_extreme'] else '保持'}")
                         if kind == "bias" else self.tr(f"{r['facet']} · {' − '.join(r['difference']) if r['difference'] else self.tr('全水準')} · {100*r['level']:g}%"))
                widget.insert("end", label)
            if widget.size():
                widget.selection_set(0); self.edit_request(kind)
        self.input_error.set("")

    def number(self, variable, label, integer=False):
        try:
            return int(variable.get()) if integer else float(variable.get())
        except (ValueError, TypeError):
            raise ValueError(self.tr(f"{label}を{'整数' if integer else '数値'}で入力してください。")) from None

    def settings(self):
        from copy import deepcopy
        import numpy as np
        from mfrm_app.desktop_i18n import text
        from mfrm_app.desktop_analysis import _bias_requests, _interval_requests, _gd_settings
        from mfrm_app.basic_analysis import normalize_anchor_records, normalize_group_anchor_records
        standalone = self.workflow.get() in {"G/D-Study", "G/D Study"}
        gd = None
        if standalone or self.gd_enabled.get():
            facets = [self.gd_facets.get(i) for i in self.gd_facets.curselection()]
            columns = [self.gd_object.get(), self.gd_score.get(), *facets]
            if not 1 <= len(facets) <= 2 or len(set(columns)) != len(columns) or any(c not in self.data for c in columns):
                raise ValueError(self.tr("G/D-Studyの測定対象・得点・ランダム相1〜2列を、重複なく選択してください。"))
            fixed_cond = self._gd_fixed_condition_dict()
            values = [v.get().strip() for v in self.gd_counts[:len(facets)]]
            custom_mode = text(self.gd_plan_mode.get(), 'ja') == "計画数を指定"
            selected = (
                {f: self.number(v, self.tr(f"計画の {f} の数"), integer=True) for f, v in zip(facets, self.gd_counts)}
                if (custom_mode or any(values)) else None
            )
            gd_raw = dict(
                object_facet=columns[0], score_col=columns[1], random_facets=facets,
                n_max_per_facet=self.number(self.gd_limit, self.tr("計画表の上限"), integer=True),
                target_population=self.gd_population.get(), selected_plan=selected,
            )
            if fixed_cond is not None:
                gd_raw["fixed_condition"] = fixed_cond
            gd = _gd_settings(gd_raw)
        if standalone:
            return dict(workflow="gd", gd_settings=gd)
        s = dict(person_col=self.person.get(), score_col=self.score.get(), facet_cols=[self.facets.get(i) for i in self.facets.curselection()],
            rating_min=self.number(self.minimum, self.tr("採点基準の最低点"), integer=True), rating_max=self.number(self.maximum, self.tr("採点基準の最高点"), integer=True), keep_original=True)
        s.update({k: v.get() for k, v in self.model_vars.items() if k in ("model", "method")})
        s.update(maxit=self.number(self.model_vars['maxit'], self.tr('最大反復数'), integer=True), reltol=self.number(self.model_vars['reltol'], self.tr('停止許容誤差')))
        if s['maxit'] < 1 or not np.isfinite(s['reltol']) or s['reltol'] <= 0:
            raise ValueError(self.tr("反復数・許容誤差には正の値を入力してください。"))
        if s['model'] == 'PCM':
            s['step_facet'] = self.model_vars['step_facet'].get()
            if s['step_facet'] not in s['facet_cols']:
                raise ValueError(self.tr("PCMの段階相を選択してください。"))
        if s['method'] == 'MML':
            s.update(quad_points=self.number(self.model_vars['quad_points'], self.tr('求積点数'), integer=True), population_prior_sd=self.number(self.model_vars['population_prior_sd'], self.tr('固定母集団SD')),
                     estimate_population_sd=False, mml_engine='EM')
            if s['quad_points'] < 3 or not np.isfinite(s['population_prior_sd']) or s['population_prior_sd'] <= 0:
                raise ValueError(self.tr("MMLの求積点数は3以上、固定SDは正の有限値にしてください。"))
        norm_anchors = normalize_anchor_records(
            deepcopy(self.anchor_records), person_col=s['person_col'], facet_cols=s['facet_cols']
        )
        norm_group_anchors = normalize_group_anchor_records(
            deepcopy(self.group_anchor_records), person_col=s['person_col'], facet_cols=s['facet_cols'],
            fixed_records=norm_anchors,
        )
        if norm_anchors:
            s['anchors'] = norm_anchors
        if norm_group_anchors:
            s['group_anchors'] = norm_group_anchors
        request = dict(workflow='mfrm', **s, bias_requests=_bias_requests(deepcopy(self.bias_requests), ['Person', *s['facet_cols']]),
            interval_requests=_interval_requests(deepcopy(self.interval_requests), s['facet_cols'], s['method']), gd_settings=gd)
        if self.joint_enabled.get():
            from mfrm_app import joint_interaction
            if s['method'] != 'MML' or s.get('population_prior_sd') != 1.:
                raise ValueError(self.tr('共同交互作用にはMML・固定母集団SD 1を指定してください。'))
            difference = [[self.joint_cell_vars[k].get() for k in pair] for pair in (('a0','b0'),('a1','b1'))] if self.joint_mode.get() != self.tr('全セル') else None
            if difference and any(self.joint_cell_vars[k].get() not in self.joint_cell_boxes[k]['values'] for k in self.joint_cell_vars):
                raise ValueError(self.tr('指定差には相A・Bの観測された水準を選択してください。'))
            request['joint_settings'] = joint_interaction.settings(dict(facet_a=self.joint_a.get(), facet_b=self.joint_b.get(),
                level=self.number(self.joint_level, self.tr('共同交互作用の信頼水準')), difference=difference), s['facet_cols'], s['method'])
        return request

    def restore_settings(self, analysis):
        from copy import deepcopy
        from mfrm_app.desktop_outputs import request_settings
        request = request_settings(analysis)
        joint = request.get('joint_settings'); self.joint_enabled.set(bool(joint))
        self.joint_a.set(joint['facet_a'] if joint else '')
        self.joint_b.set(joint['facet_b'] if joint else '')
        self.joint_level.set(str(joint['level'] if joint else .95))
        self.joint_mode.set(self.tr('全セル＋指定差') if joint and joint['difference'] else self.tr('全セル'))
        cells = joint['difference'] if joint and joint['difference'] else [['',''],['','']]
        for key, value in zip(('a0','b0','a1','b1'), [v for cell in cells for v in cell]):
            self.joint_cell_vars[key].set(value)
        self.workflow.set("G/D-Study" if request['workflow'] == 'gd' else "MFRM")
        if request['workflow'] == 'mfrm':
            self.person.set(request['person_col']); self.score.set(request['score_col'])
            self.minimum.set(str(request['rating_min'])); self.maximum.set(str(request['rating_max']))
            candidates = [self.facets.get(i) for i in range(self.facets.size())]
            self.facets.delete(0, 'end')
            for i, column in enumerate([*request['facet_cols'], *[c for c in candidates if c not in request['facet_cols']]]):
                self.facets.insert('end', column)
                if column in request['facet_cols']:
                    self.facets.selection_set(i)
            for k, v in self.model_vars.items():
                v.set(str(request.get(k, self.model_defaults[k])))
            self.anchor_records = deepcopy(request.get('anchors', []))
            self.group_anchor_records = deepcopy(request.get('group_anchors', []))
            self.refresh_anchor_lists()
            self.bias_requests, self.interval_requests = deepcopy(request['bias_requests']), deepcopy(request['interval_requests'])
            self.refresh_request_lists()
        else:
            self.anchor_records, self.group_anchor_records = [], []
            self.refresh_anchor_lists()
        gd = request.get('gd_settings'); self.gd_enabled.set(bool(gd))
        if gd:
            self.gd_object.set(gd['object_facet']); self.gd_score.set(gd['score_col'])
            self.gd_population.set(gd['target_population']); self.gd_limit.set(str(gd['n_max_per_facet']))
            fc = gd.get('fixed_condition')
            if fc:
                self._set_gd_fixed_mode_code(fc['mode'])
                self.gd_fixed_facet.set(fc['facet'])
                self.update_gd_facet_choices()
                if fc.get('level'):
                    self.gd_fixed_level.set(fc['level'])
            else:
                self._set_gd_fixed_mode_code('none')
                self.update_gd_facet_choices()
            # Preserve saved facet order; it also defines the D-Study axes.
            excluded = {gd['object_facet'], gd['score_col'], 'Weight'}
            if fc:
                excluded.add(fc['facet'])
            columns = [*gd['random_facets'], *[c for c in self.data if c not in gd['random_facets'] and c not in excluded]]
            self.gd_facets.delete(0, 'end')
            for i, c in enumerate(columns):
                self.gd_facets.insert('end', c)
                if c in gd['random_facets']:
                    self.gd_facets.selection_set(i)
            self.gd_plan_mode.set(self.tr("計画数を指定") if gd['selected_plan'] else self.tr("観測水準数を使う"))
            for i in range(2):
                f = gd['random_facets'][i] if i < len(gd['random_facets']) else None
                self.gd_counts[i].set(str(gd['selected_plan'][f]) if (gd['selected_plan'] and f in gd['selected_plan']) else '')
            self.update_gd_labels()
        self.refresh_choices(); self.change_workflow()

    def show_step(self, step):
        if self.busy or (step == "settings" and self.data is None) or (step == "results" and self.analysis is None):
            return
        for page in self.pages.values():
            page.pack_forget()
        self.pages[step].pack(fill="both", expand=True)
        self.step = step
        for key, button in self.step_buttons.items():
            allowed = key == "data" or (key == "settings" and self.data is not None) or (key == "results" and self.analysis is not None)
            btn_style = "Current.TButton" if key == step else ("StepReady.TButton" if allowed else "TButton")
            button.configure(state="normal" if allowed else "disabled", style=btn_style)
        self.status.set("")
        if step == "results":
            from mfrm_app.desktop_outputs import request_settings
            from mfrm_app.engine import stable_json_fingerprint
            try:
                same = stable_json_fingerprint(self.settings()) == stable_json_fingerprint(request_settings(self.analysis))
            except (ValueError, TypeError):
                same = False
            self.draft_notice.set(self.tr("設定は変更されています。表示中の結果は前回の分析です。") if not same else "")
            if same:
                self.draft_banner.pack_forget()
            else:
                self.draft_banner.pack(fill="x", before=self.notebook, pady=(0, 8))

    def update_observed_scores(self, *args):
        import pandas as pd
        if self.data is None or self.score.get() not in self.data:
            self.observed_scores.set("")
            return
        scores = pd.to_numeric(self.data[self.score.get()], errors="coerce").dropna()
        values = sorted(scores.unique())
        sample = ", ".join(f"{v:g}" for v in values[:12])
        self.observed_scores.set(self.tr("データ内の得点：") + (sample + (" …" if len(values) > 12 else "") if values else self.tr("数値なし")))

    def update_facet_choices(self, *args, reset=False):
        if self.data is None:
            return
        selected = ({"Rater", "Task", "Criterion", "評価者", "課題", "観点"} if reset else
                    {self.facets.get(i) for i in self.facets.curselection()})
        candidates = [c for c in self.data.columns if c not in {self.person.get(), self.score.get()}]
        self.facets.delete(0, "end")
        for index, column in enumerate(candidates):
            self.facets.insert("end", column)
            if column in selected:
                self.facets.selection_set(index)
        self.refresh_choices()

    def show_settings_help(self):
        from mfrm_app.desktop_guidance import SETTINGS_HELP
        from mfrm_app.desktop_i18n import GUIDES_EN
        if self.language == "en": SETTINGS_HELP = GUIDES_EN["SETTINGS_HELP"]
        return self.guide_window(self.tr("設定"), SETTINGS_HELP)

    def show_analysis_info(self):
        import json
        from mfrm_app.desktop_outputs import request_settings
        a = self.analysis
        text = json.dumps(request_settings(a), ensure_ascii=False, indent=2)
        text += self.tr("\n\n分析ID: ") + a["analysis_id"]
        text += self.tr("\n\nプロジェクトファイルには入力データと識別子を含みます。")
        text += self.tr("\n\n計算版・出典: ") + json.dumps(a["source"], ensure_ascii=False, indent=2)
        return self.guide_window(self.tr("分析情報"), text)

    def discard(self):
        from tkinter import messagebox
        return not self.unsaved or messagebox.askyesno(self.tr("未保存の結果"), self.tr("現在の結果を保存せずに進みますか？"), parent=self.root)

    def error(self, exc):
        from tkinter import messagebox
        self.status.set(self.tr("処理できませんでした。入力や保存先を確認してください。"))
        message = self.tr(str(exc))
        if isinstance(exc, PermissionError):
            message = self.tr("指定した場所へのアクセスが拒否されました。\n読込・保存先の権限を確認し、利用できる場所を選んでください。\n\n詳細: ") + message
        messagebox.showerror(self.tr("処理できませんでした"), message, parent=self.root)

    def guide_window(self, title, text, action_label=None, action=None):
        from tkinter.scrolledtext import ScrolledText
        window = self.tk.Toplevel(self.root)
        window.title(title)
        window.geometry("770x620")
        content = ScrolledText(window, wrap="word", padx=16, pady=12, font=("", 13))
        content.pack(fill="both", expand=True)
        content.insert("1.0", text)
        content.configure(state="disabled")
        buttons = self.ttk.Frame(window, padding=10)
        buttons.pack(fill="x")
        if action is not None:
            self.ttk.Button(buttons, text=action_label, command=action).pack(side="left")
        self.ttk.Button(buttons, text=self.tr("閉じる"), command=window.destroy).pack(side="right")
        return window

    def show_csv_guide(self):
        from mfrm_app.desktop_guidance import CSV_GUIDE
        from mfrm_app.desktop_i18n import GUIDES_EN
        if self.language == "en": CSV_GUIDE = GUIDES_EN["CSV_GUIDE"]
        return self.guide_window(self.tr("CSVの形式と見本"), CSV_GUIDE, self.tr("見本CSVを保存"), self.save_example)

    def save_example(self):
        from tkinter import filedialog
        from mfrm_app.desktop_guidance import EXAMPLE_PATH
        from mfrm_app.project import write_bytes_atomic
        name = filedialog.asksaveasfilename(parent=self.root, initialfile=self.tr("作文_見本.csv"), defaultextension=".csv", filetypes=[("CSV UTF-8", "*.csv")])
        if name:
            try:
                if Path(name).suffix.lower() != ".csv":
                    raise ValueError(self.tr("見本は .csv 拡張子で保存してください。"))
                write_bytes_atomic(name, EXAMPLE_PATH.read_bytes())
                self.status.set(self.tr(f"作文の人工データを保存しました：{name}"))
            except Exception as exc:
                self.error(exc)

    def show_help(self):
        from mfrm_app.desktop_guidance import QUICK_HELP, WEB_APP_URL
        from mfrm_app.desktop_i18n import GUIDES_EN
        if self.language == "en": QUICK_HELP = GUIDES_EN["QUICK_HELP"]
        return self.guide_window(self.tr("使い方・理論"), QUICK_HELP + "\n" + WEB_APP_URL, self.tr("Web版の理論・ヘルプ"), self.open_web_help)

    def open_web_help(self):
        import webbrowser
        from mfrm_app.desktop_guidance import WEB_APP_URL
        try:
            if not webbrowser.open_new_tab(WEB_APP_URL):
                raise RuntimeError(self.tr("ブラウザーを開けませんでした。ヘルプ欄のURLをコピーして開いてください。"))
        except Exception as exc:
            self.error(exc)

    def table_view(self, frame, table):
        import pandas as pd
        tree = self.ttk.Treeview(frame, columns=list(range(len(table.columns))), show="headings")
        tree.tag_configure("even", background=PALETTE["row_even"])
        tree.tag_configure("odd", background=PALETTE["row_odd"])
        vertical = self.ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        horizontal = self.ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        horizontal.pack(side="bottom", fill="x")
        vertical.pack(side="right", fill="y")
        tree.pack(fill="both", expand=True)
        sample_df = table.head(30)
        for index, column in enumerate(table.columns):
            col_name = str(column)
            tree.heading(index, text=col_name)
            max_chars = len(col_name)
            if len(sample_df) > 0:
                for val in sample_df.iloc[:, index]:
                    s_val = "" if pd.isna(val) else (f"{val:.4g}" if isinstance(val, float) else str(val))
                    if len(s_val) > max_chars:
                        max_chars = len(s_val)
            width = min(340, max(95, max_chars * 8 + 26))
            tree.column(index, width=width, stretch=len(table.columns) <= 6)
        for r_idx, row in enumerate(table.itertuples(index=False, name=None)):
            vals = ["" if pd.isna(x) else (f"{x:.4g}" if isinstance(x, float) else str(x)) for x in row]
            tree.insert("", "end", values=vals, tags=("even" if r_idx % 2 == 0 else "odd",))
        return tree

    def set_data(self, data, label):
        self.data = data
        self.data_label = label
        columns = list(data.columns)
        self.person_box.configure(values=columns)
        self.score_box.configure(values=columns)
        self.person.set(next((c for c in ("Person", "Examinee", "対象者", "受験者") if c in columns), ""))
        self.score.set(next((c for c in ("Score", "Rating", "得点", "評定") if c in columns), ""))
        self.minimum.set("")
        self.maximum.set("")
        self.workflow.set("MFRM")
        for name, value in self.model_defaults.items():
            self.model_vars[name].set(value)
        self.anchor_records, self.group_anchor_records = [], []
        self.bias_requests, self.interval_requests = [], []
        self.bias_list.delete(0, "end"); self.interval_list.delete(0, "end")
        self.joint_enabled.set(False); self.joint_a.set(''); self.joint_b.set('')
        self.joint_level.set('.95'); self.joint_mode.set(self.tr('全セル'))
        for var in self.joint_cell_vars.values():
            var.set('')
        self.gd_enabled.set(False)
        self._set_gd_fixed_mode_code("none")
        self.gd_plan_mode.set(self.tr("観測水準数を使う"))
        self.gd_swap.set(False)
        self.gd_object.set(self.person.get()); self.gd_score.set(self.score.get())
        self.gd_object_box.configure(values=columns); self.gd_score_box.configure(values=columns)
        self.gd_population.set(""); self.gd_limit.set("10")
        for v in self.gd_counts:
            v.set("")
        self.update_gd_facet_choices(reset=True)
        self.change_workflow()
        self.update_facet_choices(reset=True)
        self.refresh_anchor_lists()
        self.data_banner_var.set(self.tr(f"{label} · {len(data):,}行 × {len(columns)}列"))
        self._set_banner_kind(self.data_banner_frame, "ok")
        self.run_button.configure(state="normal")
        self.next_button.configure(state="normal")
        self.input_error.set("")
        self.clear_result()
        for widget in self.preview.winfo_children():
            widget.destroy()
        self.ttk.Label(self.preview, text=self.tr("先頭10行"), style="Muted.TLabel").pack(anchor="w", pady=(0, 6))
        self.table_view(self.preview, data.head(10))
        self.show_step("data")

    def clear_result(self):
        self.analysis = None
        self.unsaved = False
        self.qualification.set("")
        self.qualification_banner.pack_forget()
        self.draft_banner.pack_forget()
        for tab in self.notebook.tabs():
            self.notebook.nametowidget(tab).destroy()
        self.canvases.clear()
        self.figure_cache = {}
        self.current_figure = None
        self.project_button.configure(state="disabled")
        self.report_button.configure(state="disabled")
        self.step_buttons["results"].configure(state="disabled", style="TButton")

    def use_demo(self):
        if self.discard():
            from mfrm_app.desktop_guidance import writing_example
            try:
                self.set_data(writing_example(), self.tr("作文サンプル（人工データ）"))
                self.minimum.set("0")
                self.maximum.set("5")
            except Exception as exc:
                self.error(exc)

    def use_gd_demo(self):
        if self.discard():
            from mfrm_app.desktop_guidance import gd_two_facet_example
            try:
                self.set_data(gd_two_facet_example(), self.tr("G/D 2相サンプル（人工データ）"))
                self.minimum.set("0")
                self.maximum.set("5")
                self.workflow.set("G/D-Study")
                self.change_workflow()
            except Exception as exc:
                self.error(exc)

    def load_example_scenario(self, key: str):
        from copy import deepcopy
        from mfrm_app.desktop_examples import get_scenario, load_scenario_dataframe
        if not self.discard():
            return
        try:
            meta = get_scenario(key)
            df = load_scenario_dataframe(key)
            label = meta["title_ja"] if self.language == "ja" else meta["title_en"]
            if meta.get("is_wide"):
                self.open_wide_converter(
                    df, label, defaults=meta.get("wide_defaults"),
                    preset_rubric=(meta["rating_min"], meta["rating_max"]),
                    preset_facets=meta.get("facet_cols"),
                )
                return
            self.set_data(df, label)
            self.person.set(meta["person_col"])
            self.score.set(meta["score_col"])
            candidates = [c for c in df.columns if c not in {meta["person_col"], meta["score_col"]}]
            self.facets.delete(0, "end")
            for idx, col in enumerate(candidates):
                self.facets.insert("end", col)
                if col in meta["facet_cols"]:
                    self.facets.selection_set(idx)
            self.minimum.set(str(meta["rating_min"]))
            self.maximum.set(str(meta["rating_max"]))
            self.anchor_records = deepcopy(meta.get("anchors", []))
            self.group_anchor_records = deepcopy(meta.get("group_anchors", []))
            self.refresh_choices()
            self.refresh_anchor_lists()
            self.workflow.set(meta["workflow"])
            self.change_workflow()
        except Exception as exc:
            self.error(exc)

    def save_scenario_csv(self, key: str, parent=None):
        from tkinter import filedialog
        from mfrm_app.desktop_examples import EXAMPLES_DIR, get_scenario, load_scenario_dataframe
        from mfrm_app.project import write_bytes_atomic
        meta = get_scenario(key)
        name = filedialog.asksaveasfilename(
            parent=parent or self.root, initialfile=meta["filename"],
            defaultextension=".csv", filetypes=[("CSV UTF-8", "*.csv")],
        )
        if not name:
            return
        try:
            if Path(name).suffix.lower() != ".csv":
                raise ValueError(self.tr("見本は .csv 拡張子で保存してください。"))
            src = EXAMPLES_DIR / meta["filename"]
            payload = src.read_bytes() if src.exists() else load_scenario_dataframe(key).to_csv(index=False).encode("utf-8")
            write_bytes_atomic(name, payload)
            self.status.set(self.tr(f"サンプルデータを保存しました：{name}"))
        except Exception as exc:
            self.error(exc)

    def open_sample_gallery(self):
        from mfrm_app.desktop_examples import SCENARIO_CATALOG
        window = self.tk.Toplevel(self.root)
        window.title(self.tr("サンプルデータ集（人工データ・8シナリオ）"))
        window.geometry("920x660")
        window.minsize(760, 520)
        window.configure(bg=PALETTE["canvas_bg"])

        top = self.tk.Frame(window, bg=PALETTE["header_bg"], padx=18, pady=12)
        top.pack(fill="x")
        self.tk.Label(
            top, text=self.tr("サンプルデータ集（人工データ・8シナリオ）"),
            bg=PALETTE["header_bg"], fg=PALETTE["header_fg"], font=(self.font_family, 13, "bold"),
        ).pack(anchor="w")
        self.tk.Label(
            top, text=self.tr("試したい研究・実務シナリオを選択すると、列・採点基準（アンカー例ではアンカー値も）を自動設定して読み込みます。"),
            bg=PALETTE["header_bg"], fg=PALETTE["header_sub"], font=(self.font_family, 10),
        ).pack(anchor="w", pady=(3, 0))

        body = self.ttk.Frame(window, style="Canvas.TFrame", padding=(14, 10))
        body.pack(fill="both", expand=True)
        canvas = self.tk.Canvas(body, bg=PALETTE["canvas_bg"], highlightthickness=0, borderwidth=0)
        scroll = self.ttk.Scrollbar(body, orient="vertical", command=canvas.yview)
        scroll.pack(side="right", fill="y")
        canvas.pack(fill="both", expand=True)
        canvas.configure(yscrollcommand=scroll.set)
        cards_frame = self.ttk.Frame(canvas, style="Canvas.TFrame")
        win_id = canvas.create_window(0, 0, window=cards_frame, anchor="nw")
        cards_frame.bind("<Configure>", lambda e, c=canvas: c.configure(scrollregion=c.bbox("all")))
        canvas.bind("<Configure>", lambda e, c=canvas, w=win_id: c.itemconfigure(w, width=e.width))

        badge_palette = {
            "writing_3facet": ("#ECFDF5", "#065F46", "#10B981"),
            "clinical_4facet": ("#EFF6FF", "#1E40AF", "#2563EB"),
            "anchor_equating": ("#EEF2FF", "#3730A3", "#4F46E5"),
            "missing_ratings": ("#FFFBEB", "#92400E", "#F59E0B"),
            "sparse_spiral": ("#F0FDFA", "#115E59", "#0D9488"),
            "sparse_disconnected": ("#FEF2F2", "#991B1B", "#EF4444"),
            "wide_spreadsheet": ("#F5F3FF", "#5B21B6", "#7C3AED"),
            "gd_2facet": ("#ECFEFF", "#155E75", "#0891B2"),
        }
        window._scenario_buttons = {}
        for item in SCENARIO_CATALOG:
            key = item["key"]
            b_bg, b_fg, stripe_col = badge_palette.get(key, ("#EFF6FF", "#1E40AF", "#2563EB"))
            card = self.tk.Frame(cards_frame, bg=PALETTE["card_bg"], highlightbackground=PALETTE["border"], highlightthickness=1)
            card.pack(fill="x", pady=5)
            self.tk.Frame(card, bg=stripe_col, width=5).pack(side="left", fill="y")
            inner = self.tk.Frame(card, bg=PALETTE["card_bg"], padx=14, pady=10)
            inner.pack(side="left", fill="both", expand=True)
            head_row = self.tk.Frame(inner, bg=PALETTE["card_bg"])
            head_row.pack(fill="x")
            badge_txt = item["badge_ja"] if self.language == "ja" else item["badge_en"]
            title_txt = item["title_ja"] if self.language == "ja" else item["title_en"]
            desc_txt = item["summary_ja"] if self.language == "ja" else item["summary_en"]
            self.tk.Label(
                head_row, text=badge_txt, bg=b_bg, fg=b_fg,
                font=(self.font_family, 9, "bold"), padx=8, pady=2,
            ).pack(side="left", padx=(0, 10))
            self.tk.Label(
                head_row, text=title_txt, bg=PALETTE["card_bg"], fg=PALETTE["text"],
                font=(self.font_family, 11, "bold"), anchor="w",
            ).pack(side="left", fill="x")
            desc_lbl = self.tk.Label(
                inner, text=desc_txt, bg=PALETTE["card_bg"], fg=PALETTE["muted"],
                font=(self.font_family, 10), justify="left", anchor="w", wraplength=640,
            )
            desc_lbl.pack(fill="x", pady=(5, 0))
            inner.bind("<Configure>", lambda e, w=desc_lbl: w.configure(wraplength=max(220, e.width - 28)))

            right = self.tk.Frame(card, bg=PALETTE["card_bg"], padx=12, pady=10)
            right.pack(side="right", fill="y")
            load_btn = self.ttk.Button(
                right, text=self.tr("このサンプルを読み込む"), style="Primary.TButton",
                command=lambda k=key, w=window: (w.destroy(), self.load_example_scenario(k)),
            )
            load_btn.pack(fill="x", pady=(0, 4))
            self.ttk.Button(
                right, text=self.tr("選択サンプルをCSV保存…"), style="Secondary.TButton",
                command=lambda k=key, w=window: self.save_scenario_csv(k, parent=w),
            ).pack(fill="x")
            window._scenario_buttons[key] = load_btn

        footer = self.ttk.Frame(window, style="Canvas.TFrame", padding=(14, 8))
        footer.pack(fill="x")
        self.ttk.Button(footer, text=self.tr("閉じる"), style="Secondary.TButton", command=window.destroy).pack(side="right")
        return window

    @staticmethod
    def _convert_wide_dataframe(raw_df, *, id_cols, score_cols, new_facet_name, score_col_name):
        import numpy as np
        import pandas as pd
        from mfrm_app.io_tables import apply_wide_to_long_pivot
        total_cells = len(raw_df) * len(score_cols)
        melted = apply_wide_to_long_pivot(
            raw_df, id_cols=id_cols, score_cols=score_cols,
            new_facet_name=new_facet_name, score_col_name=score_col_name,
        )
        dropped = max(0, total_cells - len(melted))

        def _fmt_score(v):
            if pd.isna(v):
                return ""
            fv = float(v)
            if np.isfinite(fv) and fv.is_integer():
                return str(int(fv))
            return f"{fv:g}"

        for c in id_cols + [new_facet_name]:
            melted[c] = melted[c].astype(str)
        melted[score_col_name] = melted[score_col_name].map(_fmt_score)
        return melted, dropped

    def open_wide_converter_manual(self):
        from tkinter import filedialog
        from mfrm_app.desktop_guidance import read_desktop_csv
        if self.data is not None:
            return self.open_wide_converter(self.data, self.data_label, allow_raw_load=False)
        if not self.discard():
            return None
        name = filedialog.askopenfilename(parent=self.root, filetypes=[("CSV", "*.csv")])
        if not name:
            return None
        try:
            raw_df = read_desktop_csv(name)
            return self.open_wide_converter(raw_df, Path(name).name, allow_raw_load=True)
        except Exception as exc:
            self.error(exc)
            return None

    def open_wide_converter(self, raw_df, label, *, defaults=None, preset_rubric=None, preset_facets=None, allow_raw_load=True):
        from mfrm_app.io_tables import detect_wide_format_columns
        detected = detect_wide_format_columns(raw_df)
        cols = list(raw_df.columns)
        def_ids = list((defaults or {}).get("id_cols") or detected["probable_id_cols"] or cols[:1])
        def_scores = list((defaults or {}).get("score_cols") or detected["probable_score_cols"] or [c for c in cols if c not in def_ids])
        default_facet = (defaults or {}).get("new_facet_name") or (
            "Rater" if def_scores and all(str(c).upper().startswith("R") for c in def_scores) else "Item"
        )
        default_score_col = (defaults or {}).get("score_col_name") or "Score"

        window = self.tk.Toplevel(self.root)
        window.title(self.tr("横長（ワイド形式）CSVの縦長変換"))
        window.geometry("900x640")
        window.minsize(740, 520)
        window.configure(bg=PALETTE["canvas_bg"])

        top = self.tk.Frame(window, bg=PALETTE["header_bg"], padx=18, pady=12)
        top.pack(fill="x")
        self.tk.Label(
            top, text=self.tr("横長（ワイド形式）CSVの縦長変換"),
            bg=PALETTE["header_bg"], fg=PALETTE["header_fg"], font=(self.font_family, 13, "bold"),
        ).pack(anchor="w")
        self.tk.Label(
            top, text=self.tr("評価者や観点が列方向に並んだ表を、1行＝1評定の縦長形式へ変換します（空欄セルは自動除外）。"),
            bg=PALETTE["header_bg"], fg=PALETTE["header_sub"], font=(self.font_family, 10),
        ).pack(anchor="w", pady=(3, 0))

        body = self.ttk.Frame(window, style="Canvas.TFrame", padding=(14, 12))
        body.pack(fill="both", expand=True)

        cfg_card = self.ttk.Frame(body, style="Card.TFrame", padding=12)
        cfg_card.pack(fill="x", pady=(0, 10))
        for col_i in range(3):
            cfg_card.columnconfigure(col_i, weight=1)

        self.ttk.Label(cfg_card, text=self.tr("行に残す列（対象者・課題など）"), style="CardBold.TLabel").grid(row=0, column=0, sticky="w", padx=(0, 10))
        self.ttk.Label(cfg_card, text=self.tr("縦に展開する得点列（複数選択）"), style="CardBold.TLabel").grid(row=0, column=1, sticky="w", padx=(0, 10))
        id_list = self._styled_listbox(cfg_card, height=5, selectmode="multiple")
        id_list.grid(row=1, column=0, sticky="ew", padx=(0, 10), pady=(4, 0))
        score_list = self._styled_listbox(cfg_card, height=5, selectmode="multiple")
        score_list.grid(row=1, column=1, sticky="ew", padx=(0, 10), pady=(4, 0))
        for idx, c in enumerate(cols):
            id_list.insert("end", c)
            score_list.insert("end", c)
            if c in def_ids:
                id_list.selection_set(idx)
            if c in def_scores and c not in def_ids:
                score_list.selection_set(idx)

        names_frame = self.ttk.Frame(cfg_card)
        names_frame.grid(row=0, column=2, rowspan=2, sticky="nsew")
        facet_name_var = self.tk.StringVar(value=default_facet)
        score_name_var = self.tk.StringVar(value=default_score_col)
        self.ttk.Label(names_frame, text=self.tr("新しい相の列名")).pack(anchor="w")
        self.ttk.Entry(names_frame, textvariable=facet_name_var, width=20).pack(fill="x", pady=(3, 8))
        self.ttk.Label(names_frame, text=self.tr("新しい得点の列名")).pack(anchor="w")
        self.ttk.Entry(names_frame, textvariable=score_name_var, width=20).pack(fill="x", pady=(3, 0))

        preview_msg = self.tk.StringVar()
        preview_banner, _ = self._make_banner(body, preview_msg, kind="info")
        preview_banner.pack(fill="x", pady=(0, 8))

        preview_card = self.ttk.Frame(body, style="Card.TFrame", padding=10)
        preview_card.pack(fill="both", expand=True)

        def _current_melted():
            ids = [id_list.get(i) for i in id_list.curselection()]
            scs = [score_list.get(i) for i in score_list.curselection()]
            fn = facet_name_var.get().strip()
            sn = score_name_var.get().strip()
            if not ids or not scs or (set(ids) & set(scs)):
                raise ValueError(self.tr("行に残す列と縦に展開する得点列を重複なく1列以上ずつ選択してください。"))
            if not fn or not sn or fn == sn or fn in ids or sn in ids:
                raise ValueError(self.tr("新しい相の列名と得点の列名を入力してください（既存ID列と重複不可）。"))
            return self._convert_wide_dataframe(raw_df, id_cols=ids, score_cols=scs, new_facet_name=fn, score_col_name=sn)

        def _refresh_preview(*args):
            for ch in preview_card.winfo_children():
                ch.destroy()
            try:
                melted, dropped = _current_melted()
                preview_msg.set(self.tr(f"プレビュー：{len(melted):,}行 × {len(melted.columns)}列（空欄除外 {dropped:,}セル）"))
                self._set_banner_kind(preview_banner, "ok")
                self.table_view(preview_card, melted.head(8))
            except Exception as exc:
                preview_msg.set(str(exc))
                self._set_banner_kind(preview_banner, "warn")

        id_list.bind("<<ListboxSelect>>", _refresh_preview)
        score_list.bind("<<ListboxSelect>>", _refresh_preview)
        facet_name_var.trace_add("write", _refresh_preview)
        score_name_var.trace_add("write", _refresh_preview)
        _refresh_preview()

        def _apply_convert():
            try:
                melted, dropped = _current_melted()
                window.destroy()
                self.set_data(melted, label)
                if preset_facets:
                    cands = [self.facets.get(i) for i in range(self.facets.size())]
                    for idx, c in enumerate(cands):
                        if c in preset_facets:
                            self.facets.selection_set(idx)
                    self.refresh_choices()
                if preset_rubric is not None:
                    self.minimum.set(str(preset_rubric[0]))
                    self.maximum.set(str(preset_rubric[1]))
                self.status.set(self.tr(f"横長表を縦長へ変換しました：{len(melted):,}行 × {len(melted.columns)}列（空欄除外 {dropped:,}セル）"))
            except Exception as exc:
                self.error(exc)

        def _load_raw():
            window.destroy()
            self.set_data(raw_df, label)

        footer = self.ttk.Frame(window, style="Canvas.TFrame", padding=(14, 10))
        footer.pack(fill="x")
        convert_btn = self.ttk.Button(footer, text=self.tr("変換して読み込む"), style="Primary.TButton", command=_apply_convert)
        convert_btn.pack(side="right", padx=(8, 0))
        if allow_raw_load:
            raw_btn = self.ttk.Button(footer, text=self.tr("変換せずそのまま読み込む"), style="Secondary.TButton", command=_load_raw)
            raw_btn.pack(side="right", padx=(8, 0))
        else:
            raw_btn = None
        self.ttk.Button(footer, text=self.tr("閉じる"), style="Secondary.TButton", command=window.destroy).pack(side="left")
        window._convert_button = convert_btn
        window._raw_button = raw_btn
        window._id_list, window._score_list = id_list, score_list
        window._facet_name_var, window._score_name_var = facet_name_var, score_name_var
        return window

    def open_csv(self):
        from tkinter import filedialog
        from mfrm_app.desktop_guidance import read_desktop_csv
        from mfrm_app.io_tables import detect_wide_format_columns
        if not self.discard():
            return
        name = filedialog.askopenfilename(parent=self.root, filetypes=[("CSV", "*.csv")])
        if name:
            try:
                data = read_desktop_csv(name)
                label = Path(name).name
                lower_cols = {c.strip().lower() for c in data.columns}
                has_explicit_score = bool(lower_cols & {"score", "rating", "得点", "評定"})
                if not has_explicit_score and detect_wide_format_columns(data)["looks_wide"]:
                    self.open_wide_converter(data, label, allow_raw_load=True)
                else:
                    self.set_data(data, label)
            except Exception as exc:
                self.error(exc)

    def set_busy(self, busy):
        self.busy = busy
        for widget in self.widgets:
            state = "disabled" if busy else ("readonly" if isinstance(widget, self.ttk.Combobox) else "normal")
            widget.configure(state=state)
        self.run_button.configure(state="disabled" if busy or self.data is None else "normal")
        self.cancel_button.configure(state="normal" if busy else "disabled")
        self.back_button.configure(state="disabled" if busy else "normal")
        for key, button in self.step_buttons.items():
            available = key == "data" or (key == "settings" and self.data is not None) or (key == "results" and self.analysis is not None)
            button.configure(state="normal" if not busy and available else "disabled")
        if busy:
            self.cancel_button.pack(side="right", padx=8)
        else:
            self.cancel_button.pack_forget()
            self.refresh_choices()

    def run(self):
        if self.busy:
            return
        try:
            from tkinter import messagebox
            from mfrm_app.desktop_guidance import input_review
            self.input_error.set("")
            settings = self.settings()
            if settings["workflow"] == "mfrm":
                review = input_review(self.data, settings)
                if review["messages"] and not messagebox.askyesno(self.tr("入力データの確認"), "\n\n".join(self.tr(message) for message in review["messages"]) +
                        self.tr(f"\n\n計算対象は{review['included']:,}行です。この内容で分析しますか？"), parent=self.root):
                    return
            else:
                review = {"included": len(self.data)}
            if not self.discard():
                return
            self.worker.start(self.data, settings)
            self.set_busy(True)
            self.status.set(self.tr(f"計算中 · {review['included']:,}行"))
        except Exception as exc:
            self.input_error.set(self.tr(str(exc)))

    def cancel(self):
        self.worker.cancel()
        self.set_busy(False)
        self.status.set(self.tr("中止しました。"))

    def poll(self):
        for status, detail in self.worker.poll():
            if status == "progress":
                self.status.set({"estimating": self.tr("推定しています…"), "diagnostics": self.tr("診断と利用条件を確認しています…"),
                                 "tables": self.tr("表を準備しています…"), "bias": self.tr("バイアスを計算しています…"),
                                 "intervals": self.tr("指定区間の数値条件を確認しています…"), "generalizability": self.tr("G/D-Studyを計算しています…"),
                                 "joint_interaction": self.tr("主効果・段階・交互作用を共同推定し、区間を確認しています…")}.get(detail, detail))
            elif status == "complete":
                self.set_busy(False)
                try:
                    self.show_analysis(detail)
                    self.unsaved = True
                except Exception as exc:
                    self.error(exc)
            elif status == "failed":
                self.set_busy(False)
                self.status.set(self.tr("計算に失敗しました。入力と設定を確認してやり直してください。"))
                self.error(detail)
        self.root.after(100, self.poll)

    def show_analysis(self, analysis, view_override=None):
        from mfrm_app.desktop_outputs import current_tables, figure_choices, analysis_figure
        from mfrm_app.basic_analysis import qualification_snapshot
        from mfrm_app.visual_maps import map_caution
        tables = current_tables(analysis)
        choices = [(self.tr(title), options) for title, options in figure_choices(analysis)]
        first_figure = analysis_figure(analysis, **choices[0][1])
        self.clear_result()
        self.analysis, self.output_tables, self.plot_choices = analysis, tables, choices
        graph = self.ttk.Frame(self.notebook, style="Card.TFrame", padding=10)
        self.notebook.add(graph, text=self.tr("図"))
        controls = self.ttk.Frame(graph); controls.pack(fill="x")
        self.map_choice = self.tk.StringVar(value=choices[0][0])
        selector = self.ttk.Combobox(controls, textvariable=self.map_choice, values=[v[0] for v in choices], state="readonly", width=34)
        selector.pack(side="left"); selector.bind("<<ComboboxSelected>>", self.show_map)
        facets = list(analysis.get("diagnostics", {}).get("measures", {}).get("Facet", []))
        facets = list(dict.fromkeys(facets))
        self.map_subset = self.tk.StringVar(value=next((f for f in facets if f != "Person"), self.tr("すべて")))
        self.subset_box = self.ttk.Combobox(controls, textvariable=self.map_subset, values=[self.tr("すべて"), *facets], state="readonly", width=12)
        self.statistic = self.tk.StringVar(value="Infit MnSq")
        self.statistic_box = self.ttk.Combobox(controls, textvariable=self.statistic, values=["Infit MnSq", self.tr("近似Infit ZSTD")], state="readonly", width=16)
        self.statistic_box.bind("<<ComboboxSelected>>", self.show_map)
        self.subset_box.bind("<<ComboboxSelected>>", self.show_map)
        self.gd_swap_check = self.ttk.Checkbutton(controls, text=self.tr("横軸を入れ替え"), variable=self.gd_swap, command=self.show_map)
        self.ttk.Button(controls, text=self.tr("図を保存…"), style="Primary.TButton", command=self.save_figure).pack(side="right")
        self.ttk.Button(controls, text=self.tr("読み方"), style="Secondary.TButton", command=self.show_map_help).pack(side="right", padx=6)
        format_bar = self.ttk.Frame(graph); format_bar.pack(fill='x', pady=6)
        self.ttk.Label(format_bar, text=self.tr('図の形式'), style="CardBold.TLabel").pack(side='left', padx=(0, 6))
        for value, label in [('explained', '説明付き'), ('manuscript', '論文掲載用')]:
            self.ttk.Radiobutton(format_bar, text=self.tr(label), variable=self.presentation, value=value,
                command=self.show_map).pack(side='left', padx=5)
        self.ttk.Button(format_bar, text=self.tr('タイトル・図注'), style="Secondary.TButton", command=self.show_figure_caption).pack(side='right')
        self.hover_text = self.tk.StringVar()
        self.ttk.Label(graph, textvariable=self.hover_text, style="Muted.TLabel", wraplength=950).pack(side="bottom", fill="x", pady=(4, 0))
        self.plot_frame = self.ttk.Frame(graph); self.plot_frame.pack(fill="both", expand=True)
        table_page = self.ttk.Frame(self.notebook, style="Card.TFrame", padding=10); self.notebook.add(table_page, text=self.tr("表"))
        bar = self.ttk.Frame(table_page); bar.pack(fill="x")
        self.table_choice = self.tk.StringVar(value=f"1. {self.tr(tables[0][0])}")
        selector = self.ttk.Combobox(bar, textvariable=self.table_choice,
            values=[f"{i+1}. {self.tr(title)}" for i, (title, _, _) in enumerate(tables)], state="readonly", width=58)
        selector.current(0); selector.pack(side="left"); self.table_selector = selector
        selector.bind("<<ComboboxSelected>>", self.show_table)
        self.ttk.Button(bar, text=self.tr("CSVを保存…"), style="Primary.TButton", command=self.save_table).pack(side="right")
        self.table_note = self.tk.StringVar()
        self.ttk.Label(table_page, textvariable=self.table_note, style="Muted.TLabel", wraplength=920).pack(fill="x", pady=8)
        self.table_frame = self.ttk.Frame(table_page); self.table_frame.pack(fill="both", expand=True)
        self.show_table()
        self.figure_cache = {tuple(sorted(choices[0][1].items())): first_figure}
        view = view_override if view_override is not None else analysis.get("view_settings", {})
        self.presentation.set(view.get("presentation", self.presentation.get()))
        self.gd_swap.set(view.get("gd_swap") == "1")
        raw_fig = view.get("figure", "")
        if raw_fig == "D Study · G / Phi":
            raw_fig = "D-Study · G / Phi"
        view = dict(view, figure=self.tr(raw_fig), statistic=self.tr(view.get("statistic", "")))
        if view.get("facet") == "すべて" and view["facet"] not in facets:
            view["facet"] = self.tr("すべて")
        if view.get("figure") in [v[0] for v in choices]:
            self.map_choice.set(view["figure"])
        if view.get("facet") in [self.tr("すべて"), *facets]:
            self.map_subset.set(view["facet"])
        if view.get("statistic") in ["Infit MnSq", self.tr("近似Infit ZSTD")]:
            self.statistic.set(view["statistic"])
        self.show_map()
        message = ""
        if "result" in analysis:
            current = qualification_snapshot(analysis["result"])
            message = (self.tr("一般表の推論は保留です：") + map_caution(analysis["result"], lang=self.language) if current["calibration_hold_reason"] else "")
            if current != analysis["qualification_at_save"]:
                message += self.tr(" 保存時から利用条件が変わったため、現在の判定を適用しています。")
        if analysis.get("bias_analyses"):
            message += self.tr(" 固定較正のバイアス診断では有無の推論は保留です。")
        if 'joint_analysis' in analysis:
            message += self.tr(" 共同交互作用は別モデルの点ごとの近似推論です。利用可否は専用の表・図で確認してください。")
        if gd_entry := analysis.get("gd_analysis"):
            gd_res = gd_entry.get("result", {})
            if not gd_res.get("available"):
                gd_reason = self.tr(gd_res.get("reason", ""))
                message += self.tr(f" G/D-Studyは計算できませんでした：{gd_reason}")
        self.qualification.set(message)
        if message:
            self.qualification_banner.pack(fill="x", before=self.notebook, pady=(0, 8))
        self.result_label.configure(text=self.tr(self.data_label))
        self.project_button.configure(state="normal"); self.report_button.configure(state="normal")
        self.show_step("results")

    def show_table(self, event=None):
        _, table, note = self.output_tables[self.table_selector.current()]
        for widget in self.table_frame.winfo_children(): widget.destroy()
        self.table_note.set(self.tr(note))
        self.table_view(self.table_frame, table)

    def save_table(self):
        from tkinter import filedialog
        from mfrm_app.project import write_bytes_atomic
        from mfrm_app.desktop_outputs import current_tables
        index = self.table_selector.current()
        name = filedialog.asksaveasfilename(parent=self.root, initialfile=f"table_{index+1}.csv", defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if name:
            try:
                if Path(name).suffix.lower() != ".csv": raise ValueError(self.tr(".csv 拡張子で保存してください。"))
                title, table, note = current_tables(self.analysis)[index]
                output = table.copy()
                output["AnalysisID"] = self.analysis["analysis_id"]
                output["Table"] = title; output["Note"] = note
                write_bytes_atomic(name, output.to_csv(index=False).encode("utf-8-sig"))
                self.status.set(self.tr(f"表と利用条件を保存しました：{name}"))
            except Exception as exc: self.error(exc)

    def figure_options(self):
        options = dict(next(o for title, o in self.plot_choices if title == self.map_choice.get()))
        if options["kind"] == "pathway":
            options.update(facet=self.map_subset.get() if self.subset_box.current() != 0 else None,
                           statistic="Infit" if self.statistic.get() == "Infit MnSq" else "InfitZSTD")
        elif options["kind"] == "gd" and self.gd_swap.get():
            options["gd_swap"] = True
        options["presentation"] = self.presentation.get()
        return options

    def show_map_help(self):
        from mfrm_app.visual_maps import MAP_GUIDES, MAP_GUIDES_EN
        if self.language == "en": MAP_GUIDES = MAP_GUIDES_EN
        kind = self.figure_options()["kind"]
        notes = {
            "bias": self.tr("校正を固定した条件付きのセル推定です。主効果の厳しさとは別です。空白セルの状態と評定数は表で確認してください。色の大小からバイアスの有無・公平性を結論にしません。"),
            "joint": self.tr("主効果・段階・交互作用を同時に推定した別モデルです。区間は点ごとの近似区間です。表にはWald近似p値と、全セル・指定差を対象としたHolm補正を示します。空白や保留は「交互作用なし」を意味しません。"),
            "interval": self.tr("指定した水準・差の点ごとの近似区間です。数値条件を満たした対象だけを描きます。保留の理由は表の注記を参照してください。有限標本被覆率・個人得点の区間は保証していません。"),
            "gd": self.tr("元の観測得点のG（相対評価）とPhi（絶対評価）を、評定回数の計画別に示します。星は保存した計画です。完全均衡交差・1／2ランダム相・各セル1観測・単位重みの範囲です。負の分散成分は0として係数を算出し、符号付きの元の値も表に残します。区間は保留です。")}
        note = MAP_GUIDES[kind][1] if kind in MAP_GUIDES else notes[kind]
        if kind == "wright": note += self.tr("\n点線は選択した尺度の段階オフセットです。対象者・相の位置とは区別して読みます。PCMでは表示中の水準を確認してください。")
        return self.guide_window(self.map_choice.get(), note)

    def show_map(self, event=None):
        from mfrm_app.desktop_outputs import analysis_figure
        from mfrm_app.visual_maps import attach_map_hover
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        options = self.figure_options()
        for widget in (self.subset_box, self.statistic_box):
            if options["kind"] == "pathway": widget.pack(side="left", padx=3)
            else: widget.pack_forget()
        gd_rf = self.analysis.get("gd_analysis", {}).get("settings", {}).get("random_facets", [])
        if options["kind"] == "gd" and len(gd_rf) == 2:
            self.gd_swap_check.pack(side="left", padx=6)
        else:
            self.gd_swap_check.pack_forget()
        key = tuple(sorted(options.items()))
        if key not in self.figure_cache: self.figure_cache[key] = analysis_figure(self.analysis, **options)
        for canvas in self.canvases:
            if hasattr(canvas, "_mfrm_hover_id"): canvas.mpl_disconnect(canvas._mfrm_hover_id)
        for widget in self.plot_frame.winfo_children(): widget.destroy()
        self.canvases.clear()
        self.current_figure = self.figure_cache[key]
        canvas = FigureCanvasTkAgg(self.current_figure, master=self.plot_frame)
        canvas.draw(); canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas._mfrm_hover_id = attach_map_hover(canvas, self.current_figure, self.hover_text.set)
        self.canvases.append(canvas); self.hover_text.set("")

    def save_figure(self):
        from tkinter import filedialog
        from mfrm_app.basic_analysis import figure_bytes
        from mfrm_app.desktop_outputs import analysis_figure
        from mfrm_app.project import write_bytes_atomic
        options = self.figure_options()
        name = filedialog.asksaveasfilename(parent=self.root, initialfile=options["kind"]+".png", defaultextension=".png", filetypes=[("PNG", "*.png"), ("SVG", "*.svg")])
        if name:
            try:
                figure = analysis_figure(self.analysis, **options)
                write_bytes_atomic(name, figure_bytes(self.analysis, options["kind"], figure=figure, format=Path(name).suffix.lower().lstrip(".")))
                if self.presentation.get() == 'manuscript':
                    from mfrm_app.figure_presentation import caption_text
                    write_bytes_atomic(Path(name).with_suffix('.caption.txt'),
                        caption_text(figure._mfrm_caption, self.analysis['analysis_id']).encode('utf-8'))
                self.status.set(self.tr(f"表示中の対象・条件の図を保存しました：{name}"))
            except Exception as exc: self.error(exc)

    def open_documents(self, *names):
        if len(names) != 1:
            self.error(ValueError(self.tr("プロジェクトは一度に1つ開いてください。")))
            return
        # Use Tk's event queue for both startup and already-running Finder opens.
        self.root.after_idle(self.open_project, names[0])

    def open_project(self, name=None):
        from tkinter import filedialog
        from mfrm_app.project import load_project
        from mfrm_app.desktop_outputs import current_tables, figure_choices, analysis_figure
        if self.busy:
            self.error(ValueError(self.tr("分析の完了または中止後にプロジェクトを開いてください。")))
            return
        if not self.discard(): return
        if name is None:
            name = filedialog.askopenfilename(parent=self.root, filetypes=[(self.tr("MFRMプロジェクト"), "*.mfrm")])
        if name:
            try:
                analysis = load_project(name)
                current_tables(analysis)
                analysis_figure(analysis, **figure_choices(analysis)[0][1])
                self.set_data(analysis["input_data"], Path(name).name)
                self.restore_settings(analysis)
                self.show_analysis(analysis)
            except Exception as exc: self.error(exc)

    def save(self):
        from tkinter import filedialog
        from mfrm_app.project import save_project
        name = filedialog.asksaveasfilename(parent=self.root, defaultextension=".mfrm", filetypes=[(self.tr("MFRMプロジェクト"), "*.mfrm")])
        if name:
            try:
                if Path(name).suffix.lower() != ".mfrm":
                    raise ValueError(self.tr("プロジェクトは .mfrm 拡張子で保存してください。"))
                from mfrm_app.project import content_fingerprint
                saved = dict(self.analysis, view_settings=self.view_settings())
                saved["content_fingerprint"] = content_fingerprint(saved)
                save_project(name, saved)
                self.analysis = saved
                self.unsaved = False
                self.status.set(self.tr(f"プロジェクトを保存しました：{name}（入力データを含む）"))
            except Exception as exc:
                self.error(exc)

    def export_report(self):
        from tkinter import filedialog
        from mfrm_app.desktop_outputs import report_html
        from mfrm_app.project import write_bytes_atomic
        name = filedialog.asksaveasfilename(parent=self.root, defaultextension=".html", filetypes=[("HTML", "*.html")])
        if name:
            try:
                if Path(name).suffix.lower() != ".html":
                    raise ValueError(self.tr("レポートは .html 拡張子で保存してください。"))
                statistic = "Infit" if self.statistic.get() == "Infit MnSq" else "InfitZSTD"
                facet = self.map_subset.get() if self.subset_box.current() != 0 else None
                write_bytes_atomic(name, report_html(self.analysis, pathway_statistic=statistic, pathway_facet=facet, presentation=self.presentation.get(), language=self.language, gd_swap=self.gd_swap.get()).encode("utf-8"))
                self.status.set(self.tr(f"レポートを保存しました：{name}"))
            except Exception as exc:
                self.error(exc)

    def close(self):
        if self.discard():
            self.worker.close()
            self.root.destroy()


def main():
    import argparse
    parser = argparse.ArgumentParser(description="MFRM Desktop prototype")
    parser.add_argument("project", nargs="?", type=Path, help="Open a saved .mfrm project")
    parser.add_argument("--smoke-test", type=Path, help="Run the worker/save smoke check in this output directory")
    parser.add_argument("--check-project", type=Path, help="Restore without any estimation (requires --output-dir)")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.smoke_test or args.check_project:
        if args.check_project and not args.output_dir:
            parser.error("--check-project requires --output-dir")
        output_dir = args.smoke_test or args.output_dir
        try:
            if args.smoke_test:
                smoke_check(args.smoke_test)
            else:
                smoke_reopen(args.check_project, args.output_dir)
        except Exception:
            # A windowed Windows executable has no stderr; avoid a blocking error dialog in CI.
            import traceback
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
            raise SystemExit(1) from None
        return
    import tkinter as tk
    root = tk.Tk()
    app = DesktopApp(root)
    if args.project:
        root.after_idle(app.open_project, args.project)
    root.mainloop()


def smoke_check(output_dir):
    """Exercise the real spawned worker, cancellation, failure and restart."""
    import json
    import time
    from mfrm_app.project import save_project
    from mfrm_app.basic_analysis import figure_bytes
    from mfrm_app.desktop_outputs import report_html
    from mfrm_app.desktop_guidance import writing_example
    from mfrm_app.desktop_examples import EXAMPLES_DIR, SCENARIO_CATALOG, load_scenario_dataframe
    output_dir.mkdir(parents=True, exist_ok=True)
    worker = AnalysisWorker()
    example = writing_example()
    assert example.shape == (960, 5)
    assert len(SCENARIO_CATALOG) == 8
    for item in SCENARIO_CATALOG:
        assert (EXAMPLES_DIR / item["filename"]).is_file(), item["filename"]
        assert len(load_scenario_dataframe(item["key"])) > 0
    data = demonstration_data()
    data.loc[1, "Score"] = float("nan")  # Exercise filtered-row alignment in the frozen worker.
    settings = dict(person_col="Person", facet_cols=["Rater", "Task"], score_col="Score", rating_min=0, rating_max=3)
    events = []
    def drain():
        terminal = None
        while worker.process is not None:
            for status, detail in worker.poll():
                events.append(status if status != "progress" else detail)
                if status in {"complete", "failed"}:
                    terminal = (status, detail)
            if worker.process is not None:
                time.sleep(.05)
        return terminal
    try:
        worker.start(data, settings)
        abandoned_path = worker.output_path
        assert worker.cancel()[0] == "cancelled" and not abandoned_path.exists()
        events.append("cancelled")
        worker.start(data, dict(settings, rating_max=0))
        status, detail = drain()
        assert status == "failed", detail
        worker.start(data, settings)
        status, analysis = drain()
        assert status == "complete", analysis
        assert len(analysis["scorefile"]) == len(data) - 1 and 1 not in analysis["scorefile"].index
        save_project(output_dir / "作文.mfrm", analysis)
        report = report_html(analysis)
        assert report.count("<svg") == 3
        (output_dir / "report.html").write_text(report, encoding="utf-8")
        (output_dir / "wright.png").write_bytes(figure_bytes(analysis, "wright", format="png"))
        record = dict(analysis_id=analysis["analysis_id"], events=events, worker_result="complete",
                      optimizer_success=bool(analysis["result"]["opt"].success),
                      bundled_example_rows=len(example), report_svg_count=3)
        record['extended'] = smoke_extended(worker, drain, output_dir)
        (output_dir / "smoke.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    finally:
        worker.close()


def smoke_extended(worker, drain, output_dir):
    """Check the shipped workflows in the frozen worker, not only the basic demo."""
    import pandas as pd
    from mfrm_app import project, desktop_outputs, structural_intervals, joint_interaction
    from mfrm_app.basic_analysis import figure_bytes
    data = pd.DataFrame([dict(Person=f'{p:03d}', Rater=f'R{r}', Task=f'T{t}', Score=(p+2*r+t)%3)
                         for p in range(30) for r in range(3) for t in range(2)])
    records = []
    for model, method in (('RSM','JMLE'), ('PCM','JMLE'), ('RSM','MML'), ('PCM','MML'), ('G/D','ANOVA')):
        gd = dict(object_facet='Person', random_facets=['Rater','Task'], selected_plan={'Rater':4,'Task':3})
        if model == 'G/D':
            settings = dict(workflow='gd', gd_settings=gd)
        else:
            settings = dict(person_col='Person', facet_cols=['Rater','Task'], score_col='Score',
                rating_min=0, rating_max=2, model=model, method=method, maxit=100,
                step_facet='Task' if model=='PCM' else None, quad_points=61 if method=='MML' else 15,
                bias_requests=[dict(facet_a='Rater',facet_b='Task')],gd_settings=gd)
            if method == 'MML':
                settings.update(interval_requests=[dict(facet='Rater',difference=['R0','R2'],level=.9)],
                                joint_settings=dict(facet_a='Rater',facet_b='Task',difference=[['R0','T0'],['R2','T1']]))
        worker.start(data, settings)
        status, analysis = drain()
        assert status == 'complete', analysis
        checks = dict(model=model, method=method, gd_available=analysis['gd_analysis']['result']['available'])
        assert checks['gd_available']
        if method == 'MML':
            ordinary = structural_intervals.read_structural_intervals(analysis['result'], analysis['interval_analyses'][0]['record'])
            joint = joint_interaction.read_joint_interaction(analysis['result'], analysis['joint_analysis']['record'])
            assert ordinary.CIEligible.all() and joint.CIEligible.all()
            checks.update(structural_targets=len(ordinary), joint_targets=len(joint),
                          structural_eligible=int(ordinary.CIEligible.sum()), joint_eligible=int(joint.CIEligible.sum()))
        destination = output_dir / (model.replace('/','')+'_'+method)
        destination.mkdir(exist_ok=True)
        project.save_project(destination/'analysis.mfrm', analysis)
        report = desktop_outputs.report_html(analysis)
        (destination/'report.html').write_text(report, encoding='utf-8')
        choices = desktop_outputs.figure_choices(analysis)
        assert report.count('<svg') == len(choices) and '<script' not in report
        for i, (title, options) in enumerate(choices):
            figure = desktop_outputs.analysis_figure(analysis, **options)
            for format in ('png','svg'):
                (destination/f'figure_{i}.{format}').write_bytes(figure_bytes(analysis,options['kind'],figure=figure,format=format))
        for i, (title, table, note) in enumerate(desktop_outputs.current_tables(analysis)):
            table.assign(Table=title, Note=note).to_csv(destination/f'table_{i}.csv',index=False)
        checks['figures'] = len(choices)
        records.append(checks)
    return records


def smoke_reopen(path, output_dir):
    import json
    from mfrm_app import engine, project, desktop_outputs, bias_analysis, generalizability, structural_intervals, joint_interaction, mml_joint
    def forbidden(*args, **kwargs):
        raise AssertionError("Reopening must not estimate or recompute diagnostics")
    engine.mfrm_estimate = engine.mfrm_diagnostics = forbidden
    bias_analysis.estimate_bias_interaction = generalizability.compute_generalizability_from_data = forbidden
    structural_intervals.compute_structural_intervals = joint_interaction.fit_joint_interaction = forbidden
    mml_joint.GroupedMmlLikelihood.evaluate = forbidden
    analysis = project.load_project(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "report.html").write_text(desktop_outputs.report_html(analysis), encoding="utf-8")
    (output_dir / "reopen.json").write_text(json.dumps({"analysis_id": analysis["analysis_id"], "estimation_calls": 0}), encoding="utf-8")


if __name__ == "__main__":
    mp.freeze_support()
    main()
