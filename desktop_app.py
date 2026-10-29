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


class DesktopApp:
    def __init__(self, root, language="en"):
        import tkinter as tk
        from tkinter import ttk
        from mfrm_app.engine import APP_VERSION

        self.language = language
        self.root, self.tk, self.ttk = root, tk, ttk
        self.worker = AnalysisWorker()
        self.analysis = self.data = None
        self.unsaved = self.busy = False
        self.widgets, self.canvases = [], []
        self.figure_cache = {}
        self.data_label = ""
        root.title(self.tr(f"MFRM Desktop 試作 — {APP_VERSION}"))
        root.geometry("1060x760")
        root.minsize(800, 640)
        root.protocol("WM_DELETE_WINDOW", self.close)
        frame = ttk.Frame(root, padding=16)
        frame.pack(fill="both", expand=True)
        navigation = ttk.Frame(frame)
        navigation.pack(fill="x", pady=(0, 16))
        ttk.Style(root).configure("Current.TButton", font=("", 13, "bold"))
        self.step_buttons = {}
        for key, label in (("data", self.tr("データ")), ("settings", self.tr("設定")), ("results", self.tr("結果"))):
            button = ttk.Button(navigation, text=label, command=lambda step=key: self.show_step(step))
            button.pack(side="left")
            self.step_buttons[key] = button
            if key != "results":
                ttk.Label(navigation, text="  ›  ").pack(side="left")
        ttk.Button(navigation, text=self.tr("ヘルプ"), command=self.show_help).pack(side="right")
        self.language_choice = tk.StringVar(value='English' if language == 'en' else '日本語')
        language_box = ttk.Combobox(navigation, textvariable=self.language_choice,
            values=['English', '日本語'], state='readonly', width=10)
        language_box.pack(side='right', padx=8)
        language_box.bind('<<ComboboxSelected>>', self.change_language)
        self.widgets.append(language_box)
        self.presentation = tk.StringVar(value='explained')
        self.status = tk.StringVar()
        self.status_label = self.wrapped_label(frame, textvariable=self.status)
        self.status_label.pack(side="bottom", fill="x", pady=(8, 0))
        self.pages = {key: ttk.Frame(frame) for key in self.step_buttons}

        # Data: choose a source, inspect the preview, then continue.
        page = self.pages["data"]
        toolbar = ttk.Frame(page)
        toolbar.pack(fill="x", pady=(0, 12))
        for label, command in ((self.tr("CSVを開く"), self.open_csv), (self.tr("サンプル"), self.use_demo), (self.tr("プロジェクトを開く"), self.open_project)):
            button = ttk.Button(toolbar, text=label, command=command)
            button.pack(side="left", padx=(0, 8))
            self.widgets.append(button)
        ttk.Button(toolbar, text=self.tr("CSV形式"), command=self.show_csv_guide).pack(side="right")
        self.input_label = ttk.Label(page, text="")
        self.input_label.pack(anchor="w", pady=(0, 8))
        bottom = ttk.Frame(page)
        bottom.pack(side="bottom", fill="x", pady=(12, 0))
        self.next_button = ttk.Button(bottom, text=self.tr("次へ"), command=lambda: self.show_step("settings"), state="disabled")
        self.next_button.pack(side="right")
        self.preview = ttk.Frame(page)
        self.preview.pack(fill="both", expand=True)

        # Settings: the scoring rubric owns the range; observed scores are context.
        page = self.pages["settings"]
        ttk.Button(page, text=self.tr("設定の説明"), command=self.show_settings_help).pack(anchor="e", pady=(0, 8))
        self.workflow = tk.StringVar(value="MFRM")
        choice = ttk.Combobox(page, textvariable=self.workflow, values=["MFRM", "G/D Study"], state="readonly", width=22)
        choice.pack(anchor="w", pady=(0, 8))
        choice.bind("<<ComboboxSelected>>", self.change_workflow)
        self.widgets.append(choice)
        self.settings_tabs = ttk.Notebook(page)
        self.settings_tabs.pack(fill="both", expand=True)
        self.settings_pages, self.settings_contents = {}, {}
        for key, title in (("data", self.tr("列・採点基準")), ("model", self.tr("モデル・推定")), ("extras", self.tr("バイアス・区間")), ("gd", self.tr("G/D計画"))):
            panel = ttk.Frame(self.settings_tabs)
            canvas = tk.Canvas(panel, highlightthickness=0, borderwidth=0)
            scrollbar = ttk.Scrollbar(panel, orient="vertical", command=canvas.yview)
            scrollbar.pack(side="right", fill="y")
            canvas.pack(fill="both", expand=True)
            canvas.configure(yscrollcommand=scrollbar.set)
            content = ttk.Frame(canvas, padding=12)
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
            box = ttk.Combobox(columns, textvariable=variable, state="readonly", width=21)
            box.grid(row=1, column=index, sticky="w", padx=(0, 18), pady=(4, 16))
            self.widgets.append(box)
            if index == 0:
                self.person_box = box
            else:
                self.score_box = box
        ttk.Label(columns, text=self.tr("評価者・課題など")).grid(row=2, column=0, sticky="w")
        facet_frame = ttk.Frame(columns)
        facet_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(4, 0))
        self.facets = tk.Listbox(facet_frame, selectmode="multiple", exportselection=False, height=5)
        self.facets.pack(side="left", fill="x", expand=True)
        scroll = ttk.Scrollbar(facet_frame, orient="vertical", command=self.facets.yview)
        scroll.pack(side="right", fill="y")
        self.facets.configure(yscrollcommand=scroll.set)
        self.facets.bind("<<ListboxSelect>>", self.refresh_choices)
        self.widgets.append(self.facets)
        rubric = ttk.LabelFrame(page, text=self.tr("採点基準"), padding=16)
        rubric.pack(fill="x", pady=16)
        self.minimum, self.maximum = tk.StringVar(), tk.StringVar()
        self.range_entries = []
        for index, (label, variable) in enumerate(((self.tr("最低点"), self.minimum), (self.tr("最高点"), self.maximum))):
            ttk.Label(rubric, text=label).grid(row=0, column=index, sticky="w", padx=(0, 18))
            entry = ttk.Entry(rubric, textvariable=variable, width=10)
            entry.grid(row=1, column=index, sticky="w", padx=(0, 18), pady=(4, 0))
            self.widgets.append(entry)
            self.range_entries.append(entry)
        self.observed_scores = tk.StringVar()
        ttk.Label(rubric, textvariable=self.observed_scores).grid(row=2, column=0, columnspan=2, sticky="w", pady=(8, 0))
        self.score.trace_add("write", self.update_observed_scores)
        self.person.trace_add("write", self.update_facet_choices)
        self.score.trace_add("write", self.update_facet_choices)
        page = self.pages["settings"]
        self.input_error = tk.StringVar()
        self.wrapped_label(page, textvariable=self.input_error, foreground="#a12b2b").pack(fill="x")
        actions = ttk.Frame(page)
        actions.pack(side="bottom", fill="x", pady=(12, 0))
        self.back_button = ttk.Button(actions, text=self.tr("戻る"), command=lambda: self.show_step("data"))
        self.back_button.pack(side="left")
        self.run_button = ttk.Button(actions, text=self.tr("分析"), command=self.run, state="disabled")
        self.run_button.pack(side="right")
        self.cancel_button = ttk.Button(actions, text=self.tr("中止"), command=self.cancel, state="disabled")

        self.settings_tabs.pack_forget()
        self.settings_tabs.pack(fill="both", expand=True)
        self.build_extended_settings()

        # Results: only result actions and outputs remain on this page.
        page = self.pages["results"]
        actions = ttk.Frame(page)
        actions.pack(fill="x", pady=(0, 12))
        self.result_label = ttk.Label(actions, text="")
        self.result_label.pack(side="left")
        self.project_button = ttk.Button(actions, text=self.tr("プロジェクトを保存"), command=self.save, state="disabled")
        self.project_button.pack(side="right", padx=(8, 0))
        self.report_button = ttk.Button(actions, text=self.tr("レポートを保存"), command=self.export_report, state="disabled")
        self.report_button.pack(side="right", padx=(8, 0))
        ttk.Button(actions, text=self.tr("分析情報"), command=self.show_analysis_info).pack(side="right", padx=(8, 0))
        self.qualification = tk.StringVar()
        self.qualification_label = self.wrapped_label(page, textvariable=self.qualification, foreground="#8f3f16")
        self.draft_notice = tk.StringVar()
        self.draft_label = self.wrapped_label(page, textvariable=self.draft_notice, foreground="#8f3f16")
        self.notebook = ttk.Notebook(page)
        self.notebook.pack(fill="both", expand=True)
        self.show_step("data")
        if root.tk.call("tk", "windowingsystem") == "aqua":
            root.createcommand("::tk::mac::OpenDocument", self.open_documents)
        root.after(100, self.poll)

    def tr(self, message):
        from mfrm_app.desktop_i18n import text
        return text(message, getattr(self, 'language', 'ja'))

    def view_settings(self):
        from mfrm_app.desktop_i18n import text
        return dict(figure=text(self.map_choice.get(), 'ja'), facet=('すべて' if self.subset_box.current() == 0 else self.map_subset.get()),
                    statistic=text(self.statistic.get(), 'ja'), presentation=self.presentation.get(), language=self.language)

    def change_language(self, event=None):
        if self.busy: return
        view = self.view_settings() if self.analysis is not None else None
        self.language = 'ja' if self.language_choice.get() == '日本語' else 'en'
        def visit(widget):
            if widget.winfo_class() in ('TLabel', 'TButton', 'TCheckbutton', 'TRadiobutton', 'TLabelframe'):
                widget.configure(text=self.tr(str(widget.cget('text'))))
            if widget.winfo_class() == 'TNotebook':
                for tab in widget.tabs(): widget.tab(tab, text=self.tr(widget.tab(tab, 'text')))
            for child in widget.winfo_children(): visit(child)
        visit(self.root)
        self.root.title(self.tr(self.root.title()))
        # Translate only fixed application choices, never user column/level names.
        for variable in (self.interval_mode, self.joint_mode): variable.set(self.tr(variable.get()))
        for panel in self.settings_contents.values():
            def choices(widget):
                if widget.winfo_class() == 'TCombobox' and str(widget.cget('textvariable')) in {str(self.interval_mode), str(self.joint_mode)}:
                    widget.configure(values=[self.tr(v) for v in widget.cget('values')])
                for child in widget.winfo_children(): choices(child)
            choices(panel)
        self.status.set(self.tr(self.status.get()))
        self.input_error.set(self.tr(self.input_error.get()))
        self.update_observed_scores(); self.update_gd_labels()
        for widget in (self.bias_list, self.interval_list):
            selected = widget.curselection()
            items = [self.tr(widget.get(i)) for i in range(widget.size())]
            widget.delete(0, 'end')
            for item in items: widget.insert('end', item)
            for index in selected: widget.selection_set(index)
        self.input_label.configure(text=self.tr(str(self.input_label.cget('text'))))
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
            ("step_facet", self.tr("PCMの段階ファセット"), []), ("quad_points", self.tr("MMLの求積点数"), None),
            ("population_prior_sd", self.tr("MMLの固定母集団SD"), None), ("maxit", self.tr("最大反復数"), None), ("reltol", self.tr("停止許容誤差"), None)]):
            v, box = self.option_field(panel, row, label, self.model_defaults[key], choices)
            self.model_vars[key], self.model_boxes[key] = v, box
            if choices is not None: box.bind("<<ComboboxSelected>>", self.refresh_choices)
        self.wrapped_label(panel, text=self.tr("RSM：共通の段階。PCM：選んだ要因の水準別の段階。\nMMLは固定SD・EMです。指定区間は対象ごとに数値条件を確認します。")).grid(row=7, column=0, columnspan=2, sticky="ew", pady=10)
        panel = self.settings_contents["extras"]
        self.bias_requests, self.interval_requests = [], []
        bias = self.ttk.LabelFrame(panel, text=self.tr("バイアス（一覧に登録した条件で計算）"), padding=8)
        bias.pack(fill="x")
        self.bias_a, self.bias_a_box = self.option_field(bias, 0, self.tr("要因A"), choices=[])
        self.bias_b, self.bias_b_box = self.option_field(bias, 0, self.tr("要因B"), choices=[], column=2)
        self.bias_bound, _ = self.option_field(bias, 1, self.tr("探索範囲の絶対値"), "10")
        self.bias_omit = self.tk.BooleanVar(value=True)
        check = self.ttk.Checkbutton(bias, text=self.tr("極端な水準を除外"), variable=self.bias_omit)
        check.grid(row=1, column=2, columnspan=2, sticky="w"); self.widgets.append(check)
        self.bias_list = self.tk.Listbox(bias, height=3, exportselection=False)
        self.bias_list.grid(row=2, column=0, columnspan=3, sticky="ew")
        self.bias_list.bind("<<ListboxSelect>>", lambda e: self.edit_request("bias"))
        buttons = self.ttk.Frame(bias); buttons.grid(row=2, column=3, rowspan=2)
        for label, command in [(self.tr("追加"), self.add_bias), (self.tr("選択を更新"), lambda: self.add_bias(update=True)), (self.tr("選択を削除"), lambda: self.remove_request("bias"))]:
            button = self.ttk.Button(buttons, text=label, command=command); button.pack(); self.widgets.append(button)
        self.wrapped_label(bias, text=self.tr("有意性・バイアスの有無の結論は保留。"), foreground="#8f3f16").grid(row=3, column=0, columnspan=3, sticky="ew")
        box = self.ttk.LabelFrame(panel, text=self.tr("指定区間：固定SD MML（一覧に登録した条件で計算）"), padding=8)
        box.pack(fill="x", pady=10)
        self.interval_facet, self.interval_facet_box = self.option_field(box, 0, self.tr("要因"), choices=[])
        self.interval_level, _ = self.option_field(box, 0, self.tr("信頼水準"), ".95", column=2)
        self.interval_mode, mode_box = self.option_field(box, 1, self.tr("対象"), self.tr("全水準"), [self.tr("全水準"), self.tr("A − B（事前指定）")])
        self.interval_a, self.interval_a_box = self.option_field(box, 2, self.tr("水準A"), choices=[])
        self.interval_b, self.interval_b_box = self.option_field(box, 2, self.tr("水準B"), choices=[], column=2)
        self.interval_facet_box.bind("<<ComboboxSelected>>", self.refresh_choices)
        self.interval_list = self.tk.Listbox(box, height=3, exportselection=False)
        self.interval_list.grid(row=3, column=0, columnspan=3, sticky="ew")
        self.interval_list.bind("<<ListboxSelect>>", lambda e: self.edit_request("interval"))
        buttons = self.ttk.Frame(box); buttons.grid(row=3, column=3, rowspan=2)
        for label, command in [(self.tr("追加"), self.add_interval), (self.tr("選択を更新"), lambda: self.add_interval(update=True)), (self.tr("選択を削除"), lambda: self.remove_request("interval"))]:
            button = self.ttk.Button(buttons, text=label, command=command); button.pack(); self.widgets.append(button)
        self.wrapped_label(box, text=self.tr("保留対象は区間を描かず理由を保持します。"), foreground="#8f3f16").grid(row=4, column=0, columnspan=3, sticky="ew")
        box = self.ttk.LabelFrame(panel, text=self.tr("共同交互作用：主効果・段階も同時に再推定する別モデル"), padding=8)
        box.pack(fill="x", pady=10)
        self.joint_enabled = self.tk.BooleanVar(value=False)
        toggle = self.ttk.Checkbutton(box, text=self.tr("次回の分析に追加（MML・固定SD 1のみ）"), variable=self.joint_enabled)
        toggle.grid(row=0, column=0, columnspan=4, sticky="w"); self.widgets.append(toggle)
        self.joint_a, self.joint_a_box = self.option_field(box, 1, self.tr("要因A"), choices=[])
        self.joint_b, self.joint_b_box = self.option_field(box, 1, self.tr("要因B"), choices=[], column=2)
        self.joint_level, _ = self.option_field(box, 2, self.tr("信頼水準"), ".95")
        self.joint_mode, _ = self.option_field(box, 2, self.tr("対象"), self.tr("全セル"), [self.tr("全セル"), self.tr("全セル＋指定差")], column=2)
        self.joint_cell_vars, self.joint_cell_boxes = {}, {}
        for i, key in enumerate(("a0", "b0", "a1", "b1")):
            label = (self.tr("差の左側") if i < 2 else self.tr("差の右側")) + (self.tr("：A水準") if key[0] == 'a' else self.tr("：B水準"))
            var, widget = self.option_field(box, 3+i//2, label, choices=[], column=2*(i%2))
            self.joint_cell_vars[key], self.joint_cell_boxes[key] = var, widget
        for widget in (self.joint_a_box, self.joint_b_box):
            widget.bind("<<ComboboxSelected>>", self.refresh_choices)
        self.wrapped_label(box, text=self.tr("数値条件を満たすセル・指定差に点ごとの区間と近似p値を表示します。\nHolm補正の対象は、この分析の全セルと指定差です。区間は同時区間ではありません。\n元の主効果モデルと固定較正の診断結果は保持します。公平性の判定にはなりません。")).grid(row=5, column=0, columnspan=4, sticky="ew", pady=8)
        panel = self.settings_contents["gd"]
        self.gd_enabled = self.tk.BooleanVar(value=False)
        check = self.ttk.Checkbutton(panel, text=self.tr("MFRMと併せてG/Dも計算（G/D Study入口では常に実行）"), variable=self.gd_enabled)
        check.grid(row=0, column=0, columnspan=4, sticky="w"); self.widgets.append(check)
        self.gd_object, self.gd_object_box = self.option_field(panel, 1, self.tr("測定対象の列"), choices=[])
        self.gd_score, self.gd_score_box = self.option_field(panel, 1, self.tr("観測得点の列"), choices=[], column=2)
        self.ttk.Label(panel, text=self.tr("ランダム要因（1〜2列）"), wraplength=150).grid(row=2, column=0, sticky="w")
        self.gd_facets = self.tk.Listbox(panel, selectmode="multiple", exportselection=False, height=5)
        self.gd_facets.grid(row=2, column=1, sticky="ew"); self.widgets.append(self.gd_facets)
        self.gd_facets.bind("<<ListboxSelect>>", self.update_gd_labels)
        self.gd_count_labels, self.gd_counts = [], []
        for i in range(2):
            label = self.ttk.Label(panel, text=self.tr(f"計画の要因{i+1}の数"), wraplength=150)
            label.grid(row=3+i, column=0, sticky="w")
            v = self.tk.StringVar(); entry = self.ttk.Entry(panel, textvariable=v, width=16)
            entry.grid(row=3+i, column=1, sticky="ew", pady=4)
            self.gd_count_labels.append(label); self.gd_counts.append(v); self.widgets.append(entry)
        self.gd_limit, _ = self.option_field(panel, 3, self.tr("計画表の上限"), "10", column=2)
        self.gd_population, _ = self.option_field(panel, 5, self.tr("対象母集団（宣言）"))
        self.wrapped_label(panel, text=self.tr("元の観測得点を使用します。MFRM推定・採点範囲の指定は不要です。\n完全均衡交差、各セル1観測、単位重み。欠測等の設計では理由を表示します。\n計画の数が全て空欄なら観測計画。候補は1〜6・8・10（上限と観測数で限定）。\nG／Phiは点推定です。区間は保留です。")).grid(row=6, column=0, columnspan=4, sticky="ew", pady=10)
        self.refresh_choices()

    def change_workflow(self, event=None):
        gd = self.workflow.get() == "G/D Study"
        for key in ("data", "model", "extras"):
            self.settings_tabs.tab(self.settings_pages[key], state="disabled" if gd else "normal")
        self.settings_tabs.select(self.settings_pages["gd" if gd else "data"])

    def refresh_choices(self, event=None):
        if not hasattr(self, "model_boxes"): return
        facets = [self.facets.get(i) for i in self.facets.curselection()]
        self.model_boxes["step_facet"].configure(values=facets)
        self.model_boxes["step_facet"].configure(state="readonly" if self.model_vars["model"].get() == "PCM" else "disabled")
        for key in ("quad_points", "population_prior_sd"):
            self.model_boxes[key].configure(state="normal" if self.model_vars["method"].get() == "MML" else "disabled")
        for box in (self.bias_a_box, self.bias_b_box): box.configure(values=["Person", *facets])
        self.interval_facet_box.configure(values=facets)
        facet = self.interval_facet.get()
        levels = list(dict.fromkeys(self.data[facet].dropna().astype(str))) if self.data is not None and facet in self.data else []
        for box in (self.interval_a_box, self.interval_b_box): box.configure(values=levels)
        for box in (self.joint_a_box, self.joint_b_box): box.configure(values=facets)
        for key, box in self.joint_cell_boxes.items():
            name = self.joint_a.get() if key[0] == 'a' else self.joint_b.get()
            levels = list(dict.fromkeys(self.data[name].dropna().astype(str))) if self.data is not None and name in self.data else []
            box.configure(values=levels)

    def update_gd_labels(self, event=None):
        facets = [self.gd_facets.get(i) for i in self.gd_facets.curselection()]
        for i, label in enumerate(self.gd_count_labels):
            label.configure(text=self.tr(f"計画の {facets[i]} の数") if i < len(facets) else self.tr(f"計画の要因{i+1}の数"))

    def add_bias(self, update=False):
        from mfrm_app.desktop_analysis import _bias_requests
        try:
            request = dict(facet_a=self.bias_a.get(), facet_b=self.bias_b.get(), max_abs=float(self.bias_bound.get()), omit_extreme=self.bias_omit.get())
            facets = ["Person", *[self.facets.get(i) for i in self.facets.curselection()]]
            self.bias_requests = _bias_requests(self.replace_request("bias", request, update), facets)
            self.refresh_request_lists()
        except Exception as exc: self.input_error.set(self.tr(str(exc)))

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
        except Exception as exc: self.input_error.set(self.tr(str(exc)))

    def replace_request(self, kind, request, update):
        requests = list(getattr(self, kind + "_requests"))
        selected = getattr(self, kind + "_list").curselection()
        if update:
            if not selected: raise ValueError(self.tr("更新する条件を一覧から選択してください。"))
            requests[selected[0]] = request
        else:
            requests.append(request)
        return requests

    def edit_request(self, kind):
        selected = getattr(self, kind + "_list").curselection()
        if not selected: return
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
        for i in reversed(widget.curselection()): requests.pop(i)
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
        from mfrm_app.desktop_analysis import _bias_requests, _interval_requests, _gd_settings
        standalone = self.workflow.get() == "G/D Study"
        gd = None
        if standalone or self.gd_enabled.get():
            facets = [self.gd_facets.get(i) for i in self.gd_facets.curselection()]
            columns = [self.gd_object.get(), self.gd_score.get(), *facets]
            if not 1 <= len(facets) <= 2 or len(set(columns)) != len(columns) or any(c not in self.data for c in columns):
                raise ValueError(self.tr("G/Dの測定対象・得点・ランダム要因1〜2列を、重複なく選択してください。"))
            values = [v.get().strip() for v in self.gd_counts[:len(facets)]]
            selected = {f: self.number(v, self.tr(f"計画の {f} の数"), integer=True) for f, v in zip(facets, self.gd_counts)} if any(values) else None
            gd = _gd_settings(dict(object_facet=columns[0], score_col=columns[1], random_facets=facets,
                n_max_per_facet=self.number(self.gd_limit, self.tr("計画表の上限"), integer=True), target_population=self.gd_population.get(), selected_plan=selected))
        if standalone: return dict(workflow="gd", gd_settings=gd)
        s = dict(person_col=self.person.get(), score_col=self.score.get(), facet_cols=[self.facets.get(i) for i in self.facets.curselection()],
            rating_min=self.number(self.minimum, self.tr("採点基準の最低点"), integer=True), rating_max=self.number(self.maximum, self.tr("採点基準の最高点"), integer=True), keep_original=True)
        s.update({k: v.get() for k, v in self.model_vars.items() if k in ("model", "method")})
        s.update(maxit=self.number(self.model_vars['maxit'], self.tr('最大反復数'), integer=True), reltol=self.number(self.model_vars['reltol'], self.tr('停止許容誤差')))
        if s['maxit'] < 1 or not np.isfinite(s['reltol']) or s['reltol'] <= 0: raise ValueError(self.tr("反復数・許容誤差には正の値を入力してください。"))
        if s['model'] == 'PCM':
            s['step_facet'] = self.model_vars['step_facet'].get()
            if s['step_facet'] not in s['facet_cols']: raise ValueError(self.tr("PCMの段階ファセットを選択してください。"))
        if s['method'] == 'MML':
            s.update(quad_points=self.number(self.model_vars['quad_points'], self.tr('求積点数'), integer=True), population_prior_sd=self.number(self.model_vars['population_prior_sd'], self.tr('固定母集団SD')),
                     estimate_population_sd=False, mml_engine='EM')
            if s['quad_points'] < 3 or not np.isfinite(s['population_prior_sd']) or s['population_prior_sd'] <= 0: raise ValueError(self.tr("MMLの求積点数は3以上、固定SDは正の有限値にしてください。"))
        request = dict(workflow='mfrm', **s, bias_requests=_bias_requests(deepcopy(self.bias_requests), ['Person', *s['facet_cols']]),
            interval_requests=_interval_requests(deepcopy(self.interval_requests), s['facet_cols'], s['method']), gd_settings=gd)
        if self.joint_enabled.get():
            from mfrm_app import joint_interaction
            if s['method'] != 'MML' or s.get('population_prior_sd') != 1.:
                raise ValueError(self.tr('共同交互作用にはMML・固定母集団SD 1を指定してください。'))
            difference = [[self.joint_cell_vars[k].get() for k in pair] for pair in (('a0','b0'),('a1','b1'))] if self.joint_mode.get() != self.tr('全セル') else None
            if difference and any(self.joint_cell_vars[k].get() not in self.joint_cell_boxes[k]['values'] for k in self.joint_cell_vars):
                raise ValueError(self.tr('指定差には要因A・Bの観測された水準を選択してください。'))
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
        self.workflow.set("G/D Study" if request['workflow'] == 'gd' else "MFRM")
        if request['workflow'] == 'mfrm':
            self.person.set(request['person_col']); self.score.set(request['score_col'])
            self.minimum.set(str(request['rating_min'])); self.maximum.set(str(request['rating_max']))
            candidates = [self.facets.get(i) for i in range(self.facets.size())]
            self.facets.delete(0, 'end')
            for i, column in enumerate([*request['facet_cols'], *[c for c in candidates if c not in request['facet_cols']]]):
                self.facets.insert('end', column)
                if column in request['facet_cols']: self.facets.selection_set(i)
            for k, v in self.model_vars.items(): v.set(str(request.get(k, self.model_defaults[k])))
            self.bias_requests, self.interval_requests = deepcopy(request['bias_requests']), deepcopy(request['interval_requests'])
            self.refresh_request_lists()
        gd = request.get('gd_settings'); self.gd_enabled.set(bool(gd))
        if gd:
            self.gd_object.set(gd['object_facet']); self.gd_score.set(gd['score_col'])
            self.gd_population.set(gd['target_population']); self.gd_limit.set(str(gd['n_max_per_facet']))
            # Preserve saved facet order; it also defines the D Study axes.
            columns = [*gd['random_facets'], *[c for c in self.data if c not in gd['random_facets']]]
            self.gd_facets.delete(0, 'end')
            for i, c in enumerate(columns):
                self.gd_facets.insert('end', c)
                if c in gd['random_facets']: self.gd_facets.selection_set(i)
            for i, f in enumerate(gd['random_facets'][:2]): self.gd_counts[i].set(str(gd['selected_plan'][f]) if gd['selected_plan'] else '')
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
            button.configure(state="normal" if allowed else "disabled", style="Current.TButton" if key == step else "TButton")
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
                self.draft_label.pack_forget()
            else:
                self.draft_label.pack(fill="x", before=self.notebook, pady=(0, 8))

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
        tree = self.ttk.Treeview(frame, columns=list(range(len(table.columns))), show="headings")
        vertical = self.ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        horizontal = self.ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        horizontal.pack(side="bottom", fill="x")
        vertical.pack(side="right", fill="y")
        tree.pack(fill="both", expand=True)
        for index, column in enumerate(table.columns):
            tree.heading(index, text=str(column))
            tree.column(index, width=135, stretch=False)
        for row in table.itertuples(index=False, name=None):
            tree.insert("", "end", values=[f"{x:.4g}" if isinstance(x, float) else str(x) for x in row])
        return tree

    def set_data(self, data, label):
        self.data = data
        self.data_label = label
        columns = list(data.columns)
        self.person_box.configure(values=columns)
        self.score_box.configure(values=columns)
        self.person.set(next((c for c in ("Person", "対象者", "受験者") if c in columns), ""))
        self.score.set(next((c for c in ("Score", "得点") if c in columns), ""))
        self.minimum.set("")
        self.maximum.set("")
        self.workflow.set("MFRM")
        for name, value in self.model_defaults.items():
            self.model_vars[name].set(value)
        self.bias_requests, self.interval_requests = [], []
        self.bias_list.delete(0, "end"); self.interval_list.delete(0, "end")
        self.joint_enabled.set(False); self.joint_a.set(''); self.joint_b.set('')
        self.joint_level.set('.95'); self.joint_mode.set(self.tr('全セル'))
        for var in self.joint_cell_vars.values(): var.set('')
        self.gd_enabled.set(False)
        self.gd_object.set(self.person.get()); self.gd_score.set(self.score.get())
        self.gd_object_box.configure(values=columns); self.gd_score_box.configure(values=columns)
        self.gd_facets.delete(0, "end")
        for c in columns:
            self.gd_facets.insert("end", c)
        self.gd_population.set(""); self.gd_limit.set("10")
        for v in self.gd_counts: v.set("")
        self.change_workflow()
        self.update_facet_choices(reset=True)
        self.input_label.configure(text=self.tr(f"{label} · {len(data):,}行 × {len(columns)}列"))
        self.run_button.configure(state="normal")
        self.next_button.configure(state="normal")
        self.input_error.set("")
        self.clear_result()
        for widget in self.preview.winfo_children():
            widget.destroy()
        self.table_view(self.preview, data.head(10))
        self.ttk.Label(self.preview, text=self.tr("先頭10行")).pack(anchor="w")
        self.show_step("data")

    def clear_result(self):
        self.analysis = None
        self.unsaved = False
        self.qualification.set("")
        self.qualification_label.pack_forget()
        self.draft_label.pack_forget()
        for tab in self.notebook.tabs():
            self.notebook.nametowidget(tab).destroy()
        self.canvases.clear()
        self.figure_cache = {}
        self.current_figure = None
        self.project_button.configure(state="disabled")
        self.report_button.configure(state="disabled")
        self.step_buttons["results"].configure(state="disabled")

    def use_demo(self):
        if self.discard():
            from mfrm_app.desktop_guidance import writing_example
            try:
                self.set_data(writing_example(), self.tr("作文サンプル（人工データ）"))
                self.minimum.set("0")
                self.maximum.set("5")
            except Exception as exc:
                self.error(exc)

    def open_csv(self):
        from tkinter import filedialog
        from mfrm_app.desktop_guidance import read_desktop_csv
        if not self.discard():
            return
        name = filedialog.askopenfilename(parent=self.root, filetypes=[("CSV", "*.csv")])
        if name:
            try:
                data = read_desktop_csv(name)
                self.set_data(data, Path(name).name)
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
                                 "intervals": self.tr("指定区間の数値条件を確認しています…"), "generalizability": self.tr("G/Dを計算しています…"),
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
        graph = self.ttk.Frame(self.notebook, padding=6)
        self.notebook.add(graph, text=self.tr("図"))
        controls = self.ttk.Frame(graph); controls.pack(fill="x")
        self.map_choice = self.tk.StringVar(value=choices[0][0])
        selector = self.ttk.Combobox(controls, textvariable=self.map_choice, values=[v[0] for v in choices], state="readonly", width=32)
        selector.pack(side="left"); selector.bind("<<ComboboxSelected>>", self.show_map)
        facets = list(analysis.get("diagnostics", {}).get("measures", {}).get("Facet", []))
        facets = list(dict.fromkeys(facets))
        self.map_subset = self.tk.StringVar(value=next((f for f in facets if f != "Person"), self.tr("すべて")))
        self.subset_box = self.ttk.Combobox(controls, textvariable=self.map_subset, values=[self.tr("すべて"), *facets], state="readonly", width=12)
        self.statistic = self.tk.StringVar(value="Infit MnSq")
        self.statistic_box = self.ttk.Combobox(controls, textvariable=self.statistic, values=["Infit MnSq", self.tr("近似Infit ZSTD")], state="readonly", width=16)
        self.statistic_box.bind("<<ComboboxSelected>>", self.show_map)
        self.subset_box.bind("<<ComboboxSelected>>", self.show_map)
        self.ttk.Button(controls, text=self.tr("図を保存…"), command=self.save_figure).pack(side="right")
        self.ttk.Button(controls, text=self.tr("読み方"), command=self.show_map_help).pack(side="right", padx=4)
        format_bar = self.ttk.Frame(graph); format_bar.pack(fill='x', pady=4)
        self.ttk.Label(format_bar, text=self.tr('図の形式')).pack(side='left')
        for value, label in [('explained', '説明付き'), ('manuscript', '論文掲載用')]:
            self.ttk.Radiobutton(format_bar, text=self.tr(label), variable=self.presentation, value=value,
                command=self.show_map).pack(side='left', padx=5)
        self.ttk.Button(format_bar, text=self.tr('タイトル・図注'), command=self.show_figure_caption).pack(side='right')
        self.hover_text = self.tk.StringVar()
        self.ttk.Label(graph, textvariable=self.hover_text, wraplength=950).pack(side="bottom", fill="x", pady=(4, 0))
        self.plot_frame = self.ttk.Frame(graph); self.plot_frame.pack(fill="both", expand=True)
        table_page = self.ttk.Frame(self.notebook, padding=6); self.notebook.add(table_page, text=self.tr("表"))
        bar = self.ttk.Frame(table_page); bar.pack(fill="x")
        self.table_choice = self.tk.StringVar(value=tables[0][0])
        selector = self.ttk.Combobox(bar, textvariable=self.table_choice,
            values=[f"{i+1}. {title}" for i, (title, _, _) in enumerate(tables)], state="readonly", width=58)
        selector.current(0); selector.pack(side="left"); self.table_selector = selector
        selector.bind("<<ComboboxSelected>>", self.show_table)
        self.ttk.Button(bar, text=self.tr("CSVを保存…"), command=self.save_table).pack(side="right")
        self.table_note = self.tk.StringVar()
        self.ttk.Label(table_page, textvariable=self.table_note, wraplength=920).pack(fill="x", pady=8)
        self.table_frame = self.ttk.Frame(table_page); self.table_frame.pack(fill="both", expand=True)
        self.show_table()
        self.figure_cache = {tuple(sorted(choices[0][1].items())): first_figure}
        view = view_override if view_override is not None else analysis.get("view_settings", {})
        self.presentation.set(view.get("presentation", self.presentation.get()))
        view = dict(view, figure=self.tr(view.get("figure", "")), statistic=self.tr(view.get("statistic", "")))
        if view.get("facet") == "すべて" and view["facet"] not in facets:
            view["facet"] = self.tr("すべて")
        if view.get("figure") in [v[0] for v in choices]: self.map_choice.set(view["figure"])
        if view.get("facet") in [self.tr("すべて"), *facets]: self.map_subset.set(view["facet"])
        if view.get("statistic") in ["Infit MnSq", self.tr("近似Infit ZSTD")]: self.statistic.set(view["statistic"])
        self.show_map()
        message = ""
        if "result" in analysis:
            current = qualification_snapshot(analysis["result"])
            message = (self.tr("一般表の推論は保留です：") + map_caution(analysis["result"], lang=self.language) if current["calibration_hold_reason"] else "")
            if current != analysis["qualification_at_save"]:
                message += self.tr(" 保存時から利用条件が変わったため、現在の判定を適用しています。")
        if analysis.get("bias_analyses"): message += self.tr(" 固定較正のバイアス診断では有無の推論は保留です。")
        if 'joint_analysis' in analysis: message += self.tr(" 共同交互作用は別モデルの点ごとの近似推論です。利用可否は専用の表・図で確認してください。")
        if "gd_analysis" in analysis: message += self.tr(" G/Dの区間推定は保留です。")
        self.qualification.set(message)
        if message: self.qualification_label.pack(fill="x", before=self.notebook, pady=(0, 8))
        self.result_label.configure(text=self.data_label)
        self.project_button.configure(state="normal"); self.report_button.configure(state="normal")
        self.show_step("results")

    def show_table(self, event=None):
        _, table, note = self.output_tables[self.table_selector.current()]
        for widget in self.table_frame.winfo_children(): widget.destroy()
        self.table_note.set(note)
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
            "gd": self.tr("元の観測得点のG（相対評価）とPhi（絶対評価）を、評定回数の計画別に示します。星は保存した計画です。完全均衡交差・1／2ランダムファセット・各セル1観測・単位重みの範囲です。負の分散成分は0として係数を算出し、符号付きの元の値も表に残します。区間は保留です。")}
        note = MAP_GUIDES[kind][1] if kind in MAP_GUIDES else notes[kind]
        if kind == "wright": note += self.tr("\n点線は選択した尺度の段階オフセットです。対象者・ファセットの位置とは区別して読みます。PCMでは表示中の水準を確認してください。")
        return self.guide_window(self.map_choice.get(), note)

    def show_map(self, event=None):
        from mfrm_app.desktop_outputs import analysis_figure
        from mfrm_app.visual_maps import attach_map_hover
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        options = self.figure_options()
        for widget in (self.subset_box, self.statistic_box):
            if options["kind"] == "pathway": widget.pack(side="left", padx=3)
            else: widget.pack_forget()
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
                write_bytes_atomic(name, report_html(self.analysis, pathway_statistic=statistic, pathway_facet=facet, presentation=self.presentation.get(), language=self.language).encode("utf-8"))
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
    output_dir.mkdir(parents=True, exist_ok=True)
    worker = AnalysisWorker()
    example = writing_example()
    assert example.shape == (960, 5)
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
