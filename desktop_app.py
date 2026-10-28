"""Tkinter prototype: one basic RSM/JMLE analysis and durable result projects.

Run `python desktop_app.py`. Packaging and OS acceptance are documented in
docs/desktop_prototype.md; this entrypoint does not import Streamlit.
"""
from __future__ import annotations

import multiprocessing as mp
from pathlib import Path
import tempfile


def _analysis_worker(sender, output_path, data, settings):
    try:
        from mfrm_app.basic_analysis import run_basic_analysis
        from mfrm_app.project import save_project

        analysis = run_basic_analysis(data, **settings, progress=lambda phase: sender.send(("progress", phase)))
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
        while self.receiver.poll():
            try:
                status, detail = self.receiver.recv()
            except EOFError:
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
    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk
        from mfrm_app.engine import APP_VERSION

        self.root, self.tk, self.ttk = root, tk, ttk
        self.worker = AnalysisWorker()
        self.analysis = self.data = None
        self.unsaved = self.busy = False
        self.widgets, self.canvases = [], []
        self.figure_cache = {}
        self.data_label = ""
        root.title(f"MFRM Desktop 試作 — {APP_VERSION}")
        root.geometry("1060x760")
        root.minsize(800, 640)
        root.protocol("WM_DELETE_WINDOW", self.close)
        frame = ttk.Frame(root, padding=16)
        frame.pack(fill="both", expand=True)
        navigation = ttk.Frame(frame)
        navigation.pack(fill="x", pady=(0, 16))
        ttk.Style(root).configure("Current.TButton", font=("", 13, "bold"))
        self.step_buttons = {}
        for key, label in (("data", "データ"), ("settings", "設定"), ("results", "結果")):
            button = ttk.Button(navigation, text=label, command=lambda step=key: self.show_step(step))
            button.pack(side="left")
            self.step_buttons[key] = button
            if key != "results":
                ttk.Label(navigation, text="  ›  ").pack(side="left")
        ttk.Button(navigation, text="ヘルプ", command=self.show_help).pack(side="right")
        self.status = tk.StringVar()
        self.status_label = ttk.Label(frame, textvariable=self.status, wraplength=950)
        self.status_label.pack(side="bottom", fill="x", pady=(8, 0))
        self.pages = {key: ttk.Frame(frame) for key in self.step_buttons}

        # Data: choose a source, inspect the preview, then continue.
        page = self.pages["data"]
        toolbar = ttk.Frame(page)
        toolbar.pack(fill="x", pady=(0, 12))
        for label, command in (("CSVを開く", self.open_csv), ("サンプル", self.use_demo), ("プロジェクトを開く", self.open_project)):
            button = ttk.Button(toolbar, text=label, command=command)
            button.pack(side="left", padx=(0, 8))
            self.widgets.append(button)
        ttk.Button(toolbar, text="CSV形式", command=self.show_csv_guide).pack(side="right")
        self.input_label = ttk.Label(page, text="")
        self.input_label.pack(anchor="w", pady=(0, 8))
        bottom = ttk.Frame(page)
        bottom.pack(side="bottom", fill="x", pady=(12, 0))
        self.next_button = ttk.Button(bottom, text="次へ", command=lambda: self.show_step("settings"), state="disabled")
        self.next_button.pack(side="right")
        self.preview = ttk.Frame(page)
        self.preview.pack(fill="both", expand=True)

        # Settings: the scoring rubric owns the range; observed scores are context.
        page = self.pages["settings"]
        ttk.Button(page, text="設定の説明", command=self.show_settings_help).pack(anchor="e", pady=(0, 8))
        columns = ttk.LabelFrame(page, text="列", padding=16)
        columns.pack(fill="x")
        self.person, self.score = tk.StringVar(), tk.StringVar()
        for index, (label, variable) in enumerate((("対象者", self.person), ("得点", self.score))):
            ttk.Label(columns, text=label).grid(row=0, column=index, sticky="w", padx=(0, 18))
            box = ttk.Combobox(columns, textvariable=variable, state="readonly", width=21)
            box.grid(row=1, column=index, sticky="w", padx=(0, 18), pady=(4, 16))
            self.widgets.append(box)
            if index == 0:
                self.person_box = box
            else:
                self.score_box = box
        ttk.Label(columns, text="評価者・課題など").grid(row=2, column=0, sticky="w")
        facet_frame = ttk.Frame(columns)
        facet_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(4, 0))
        self.facets = tk.Listbox(facet_frame, selectmode="multiple", exportselection=False, height=5)
        self.facets.pack(side="left", fill="x", expand=True)
        scroll = ttk.Scrollbar(facet_frame, orient="vertical", command=self.facets.yview)
        scroll.pack(side="right", fill="y")
        self.facets.configure(yscrollcommand=scroll.set)
        self.widgets.append(self.facets)
        rubric = ttk.LabelFrame(page, text="採点基準", padding=16)
        rubric.pack(fill="x", pady=16)
        self.minimum, self.maximum = tk.StringVar(), tk.StringVar()
        self.range_entries = []
        for index, (label, variable) in enumerate((("最低点", self.minimum), ("最高点", self.maximum))):
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
        self.input_error = tk.StringVar()
        ttk.Label(page, textvariable=self.input_error, foreground="#a12b2b", wraplength=900).pack(anchor="w")
        actions = ttk.Frame(page)
        actions.pack(side="bottom", fill="x", pady=(12, 0))
        self.back_button = ttk.Button(actions, text="戻る", command=lambda: self.show_step("data"))
        self.back_button.pack(side="left")
        self.run_button = ttk.Button(actions, text="分析", command=self.run, state="disabled")
        self.run_button.pack(side="right")
        self.cancel_button = ttk.Button(actions, text="中止", command=self.cancel, state="disabled")

        # Results: only result actions and outputs remain on this page.
        page = self.pages["results"]
        actions = ttk.Frame(page)
        actions.pack(fill="x", pady=(0, 12))
        self.result_label = ttk.Label(actions, text="")
        self.result_label.pack(side="left")
        self.project_button = ttk.Button(actions, text="プロジェクトを保存", command=self.save, state="disabled")
        self.project_button.pack(side="right", padx=(8, 0))
        self.report_button = ttk.Button(actions, text="レポートを保存", command=self.export_report, state="disabled")
        self.report_button.pack(side="right", padx=(8, 0))
        ttk.Button(actions, text="分析情報", command=self.show_analysis_info).pack(side="right", padx=(8, 0))
        self.qualification = tk.StringVar()
        self.qualification_label = ttk.Label(page, textvariable=self.qualification, wraplength=950, foreground="#8f3f16")
        self.draft_notice = tk.StringVar()
        self.draft_label = ttk.Label(page, textvariable=self.draft_notice, wraplength=950, foreground="#8f3f16")
        self.notebook = ttk.Notebook(page)
        self.notebook.pack(fill="both", expand=True)
        self.show_step("data")
        root.after(100, self.poll)

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
            settings = self.analysis["fit_settings"]
            same = (self.person.get() == settings["person_col"] and self.score.get() == settings["score_col"]
                    and [self.facets.get(i) for i in self.facets.curselection()] == list(settings["facet_cols"])
                    and self.minimum.get() == str(settings["rating_min"]) and self.maximum.get() == str(settings["rating_max"]))
            self.draft_notice.set("設定は変更されています。表示中の結果は前回の分析です。" if not same else "")
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
        self.observed_scores.set("データ内の得点：" + (sample + (" …" if len(values) > 12 else "") if values else "数値なし"))

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

    def show_settings_help(self):
        from mfrm_app.desktop_guidance import SETTINGS_HELP
        return self.guide_window("設定", SETTINGS_HELP)

    def show_analysis_info(self):
        from mfrm_app.basic_analysis import qualification_snapshot
        a = self.analysis
        text = "RSM / JMLE\n\n" + "\n".join(f"{key}: {value}" for key, value in a["fit_settings"].items())
        text += "\n\n分析ID: " + a["analysis_id"]
        text += "\n\nプロジェクトファイルには入力データと識別子を含みます。"
        text += "\n\n" + str(qualification_snapshot(a["result"]))
        return self.guide_window("分析情報", text)

    def discard(self):
        from tkinter import messagebox
        return not self.unsaved or messagebox.askyesno("未保存の結果", "現在の結果を保存せずに進みますか？", parent=self.root)

    def error(self, exc):
        from tkinter import messagebox
        messagebox.showerror("処理できませんでした", str(exc), parent=self.root)

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
        self.ttk.Button(buttons, text="閉じる", command=window.destroy).pack(side="right")
        return window

    def show_csv_guide(self):
        from mfrm_app.desktop_guidance import CSV_GUIDE
        return self.guide_window("CSVの形式と見本", CSV_GUIDE, "見本CSVを保存", self.save_example)

    def save_example(self):
        from tkinter import filedialog
        from mfrm_app.desktop_guidance import EXAMPLE_PATH
        from mfrm_app.project import write_bytes_atomic
        name = filedialog.asksaveasfilename(parent=self.root, initialfile="作文_見本.csv", defaultextension=".csv", filetypes=[("CSV UTF-8", "*.csv")])
        if name:
            try:
                if Path(name).suffix.lower() != ".csv":
                    raise ValueError("見本は .csv 拡張子で保存してください。")
                write_bytes_atomic(name, EXAMPLE_PATH.read_bytes())
                self.status.set(f"作文の人工データを保存しました：{name}")
            except Exception as exc:
                self.error(exc)

    def show_help(self):
        from mfrm_app.desktop_guidance import QUICK_HELP, WEB_APP_URL
        return self.guide_window("使い方・理論", QUICK_HELP + "\n" + WEB_APP_URL, "Web版の理論・ヘルプ", self.open_web_help)

    def open_web_help(self):
        import webbrowser
        from mfrm_app.desktop_guidance import WEB_APP_URL
        try:
            if not webbrowser.open_new_tab(WEB_APP_URL):
                raise RuntimeError("ブラウザーを開けませんでした。ヘルプ欄のURLをコピーして開いてください。")
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
        self.update_facet_choices(reset=True)
        self.input_label.configure(text=f"{label} · {len(data):,}行 × {len(columns)}列")
        self.run_button.configure(state="normal")
        self.next_button.configure(state="normal")
        self.input_error.set("")
        self.clear_result()
        for widget in self.preview.winfo_children():
            widget.destroy()
        self.table_view(self.preview, data.head(10))
        self.ttk.Label(self.preview, text="先頭10行").pack(anchor="w")
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
                self.set_data(writing_example(), "作文サンプル（人工データ）")
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

    def run(self):
        if self.busy:
            return
        try:
            from tkinter import messagebox
            from mfrm_app.desktop_guidance import input_review
            self.input_error.set("")
            for value, widget, name in ((self.person.get(), self.person_box, "対象者"), (self.score.get(), self.score_box, "得点"),
                                         (self.facets.curselection(), self.facets, "評価者・課題など")):
                if not value:
                    widget.focus_set()
                    raise ValueError(f"「{name}」の列を選択してください。")
            bounds = []
            for label, variable, entry in zip(("最低点", "最高点"), (self.minimum, self.maximum), self.range_entries):
                try:
                    bounds.append(int(variable.get()))
                except ValueError:
                    entry.focus_set()
                    raise ValueError(f"採点基準の{label}を整数で入力してください。") from None
            minimum, maximum = bounds
            settings = dict(person_col=self.person.get(), score_col=self.score.get(),
                            facet_cols=[self.facets.get(i) for i in self.facets.curselection()],
                            rating_min=minimum, rating_max=maximum)
            review = input_review(self.data, settings)
            if review["messages"] and not messagebox.askyesno("入力データの確認", "\n\n".join(review["messages"]) +
                    f"\n\n計算対象は{review['included']:,}行です。この内容で分析しますか？", parent=self.root):
                return
            if not self.discard():
                return
            self.worker.start(self.data, settings)
            self.clear_result()
            self.set_busy(True)
            self.status.set(f"計算中 · {review['included']:,}行")
        except Exception as exc:
            self.input_error.set(str(exc))

    def cancel(self):
        self.worker.cancel()
        self.set_busy(False)
        self.status.set("中止しました。")

    def poll(self):
        for status, detail in self.worker.poll():
            if status == "progress":
                self.status.set({"estimating": "推定しています…", "diagnostics": "診断と利用条件を確認しています…",
                                 "tables": "表を準備しています…"}.get(detail, detail))
            elif status == "complete":
                self.set_busy(False)
                try:
                    self.show_analysis(detail)
                    self.unsaved = True
                except Exception as exc:
                    self.error(exc)
            elif status == "failed":
                self.set_busy(False)
                self.status.set("計算に失敗しました。入力と設定を確認してやり直してください。")
                self.error(detail)
        self.root.after(100, self.poll)

    def show_analysis(self, analysis):
        from mfrm_app.basic_analysis import current_tables, qualification_snapshot
        from mfrm_app.visual_maps import MAP_GUIDES, pathway_map_values
        tables = current_tables(analysis)
        self.clear_result()
        self.analysis = analysis
        graph = self.ttk.Frame(self.notebook, padding=6)
        self.notebook.add(graph, text="図")
        controls = self.ttk.Frame(graph)
        controls.pack(fill="x")
        self.map_choice = self.tk.StringVar(value=MAP_GUIDES["wright"][0])
        selector = self.ttk.Combobox(controls, textvariable=self.map_choice, values=[v[0] for v in MAP_GUIDES.values()], state="readonly", width=20)
        selector.pack(side="left")
        selector.bind("<<ComboboxSelected>>", self.show_map)
        self.map_subset = self.tk.StringVar(value="すべて")
        facets = list(analysis["diagnostics"]["measures"]["Facet"].drop_duplicates())
        self.subset_box = self.ttk.Combobox(controls, textvariable=self.map_subset, values=["すべて", *facets], state="disabled", width=14)
        self.statistic = self.tk.StringVar(value="Infit MnSq")
        self.statistic_box = self.ttk.Combobox(controls, textvariable=self.statistic, values=["Infit MnSq", "近似Infit ZSTD"], state="readonly", width=18)
        self.statistic_box.bind("<<ComboboxSelected>>", self.show_map)
        if facets:
            self.map_subset.set(next((f for f in facets if f != "Person"), facets[0]))
        self.subset_box.bind("<<ComboboxSelected>>", self.show_map)
        self.ttk.Button(controls, text="図を保存…", command=self.save_figure).pack(side="right")
        self.ttk.Button(controls, text="読み方", command=self.show_map_help).pack(side="right", padx=6)
        self.hover_text = self.tk.StringVar()
        self.ttk.Label(graph, textvariable=self.hover_text, wraplength=1000).pack(side="bottom", fill="x", pady=(4, 0))
        self.plot_frame = self.ttk.Frame(graph)
        self.plot_frame.pack(fill="both", expand=True)
        titles = {"Analysis Design and Retained Ratings": "データ概要", "Measures and Fit": "推定値", "Observed Category Use": "得点"}
        for title, table, note in tables:
            frame = self.ttk.Frame(self.notebook, padding=6)
            self.notebook.add(frame, text=titles.get(title, title))
            self.ttk.Button(frame, text="注記", command=lambda t=title, n=note: self.guide_window(t, n)).pack(anchor="e", pady=(0, 6))
            self.table_view(frame, table)
        self.show_map()
        current = qualification_snapshot(analysis["result"])
        message = ("推論は保留です：" + current["calibration_hold_message"] if current["calibration_hold_reason"] else "")
        if current != analysis["qualification_at_save"]:
            message += " 保存時から利用条件が変わったため、現在の判定を適用しています。"
        self.qualification.set(message)
        if message:
            self.qualification_label.pack(fill="x", before=self.notebook, pady=(0, 8))
        self.result_label.configure(text=self.data_label)
        self.project_button.configure(state="normal")
        self.report_button.configure(state="normal")
        self.show_step("results")

    def show_map_help(self):
        from mfrm_app.visual_maps import MAP_GUIDES
        kind = next(key for key, value in MAP_GUIDES.items() if value[0] == self.map_choice.get())
        return self.guide_window(MAP_GUIDES[kind][0], MAP_GUIDES[kind][1])

    def show_map(self, event=None):
        from mfrm_app.basic_analysis import basic_figure
        from mfrm_app.visual_maps import MAP_GUIDES, attach_map_hover
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        kind = next(key for key, value in MAP_GUIDES.items() if value[0] == self.map_choice.get())
        self.subset_box.configure(state="readonly" if kind == "pathway" else "disabled")
        if kind == "pathway":
            self.subset_box.pack(side="left", padx=6)
            self.statistic_box.pack(side="left")
        else:
            self.subset_box.pack_forget()
            self.statistic_box.pack_forget()
        facet = self.map_subset.get() if kind == "pathway" and self.map_subset.get() != "すべて" else None
        statistic = "Infit" if self.statistic.get() == "Infit MnSq" else "InfitZSTD"
        key = (kind, facet, statistic if kind == "pathway" else None)
        if key not in self.figure_cache:
            self.figure_cache[key] = basic_figure(self.analysis, kind, facet=facet, statistic=statistic)
            self.figure_cache[key].suptitle("")  # The selector already names the displayed figure.
        for widget in self.plot_frame.winfo_children():
            widget.destroy()
        for canvas in self.canvases:
            if hasattr(canvas, "_mfrm_hover_id"):
                canvas.mpl_disconnect(canvas._mfrm_hover_id)
        self.canvases.clear()
        self.current_figure = self.figure_cache[key]
        canvas = FigureCanvasTkAgg(self.current_figure, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas._mfrm_hover_id = attach_map_hover(canvas, self.current_figure, self.hover_text.set)
        self.canvases.append(canvas)
        self.hover_text.set("")

    def save_figure(self):
        from tkinter import filedialog
        from mfrm_app.basic_analysis import figure_bytes
        from mfrm_app.visual_maps import MAP_GUIDES
        from mfrm_app.project import write_bytes_atomic
        kind = next(key for key, value in MAP_GUIDES.items() if value[0] == self.map_choice.get())
        name = filedialog.asksaveasfilename(parent=self.root, initialfile=f"{kind}_map.png", defaultextension=".png", filetypes=[("PNG", "*.png"), ("SVG", "*.svg")])
        if name:
            try:
                format = Path(name).suffix.lower().lstrip(".")
                facet = self.map_subset.get() if kind == "pathway" and self.map_subset.get() != "すべて" else None
                statistic = "Infit" if self.statistic.get() == "Infit MnSq" else "InfitZSTD"
                write_bytes_atomic(name, figure_bytes(self.analysis, kind, format=format, facet=facet, statistic=statistic))
                self.status.set(f"表示中の図を保存しました：{name}")
            except Exception as exc:
                self.error(exc)

    def open_project(self):
        from tkinter import filedialog
        from mfrm_app.project import load_project
        if not self.discard():
            return
        name = filedialog.askopenfilename(parent=self.root, filetypes=[("MFRMプロジェクト", "*.mfrm")])
        if name:
            try:
                analysis = load_project(name)
                self.set_data(analysis["input_data"], Path(name).name)
                settings = analysis["fit_settings"]
                self.person.set(settings["person_col"])
                self.score.set(settings["score_col"])
                self.minimum.set(str(settings["rating_min"]))
                self.maximum.set(str(settings["rating_max"]))
                self.facets.selection_clear(0, "end")
                for index in range(self.facets.size()):
                    column = self.facets.get(index)
                    if column in settings["facet_cols"]:
                        self.facets.selection_set(index)
                self.show_analysis(analysis)
            except Exception as exc:
                self.error(exc)

    def save(self):
        from tkinter import filedialog
        from mfrm_app.project import save_project
        name = filedialog.asksaveasfilename(parent=self.root, defaultextension=".mfrm", filetypes=[("MFRMプロジェクト", "*.mfrm")])
        if name:
            try:
                if Path(name).suffix.lower() != ".mfrm":
                    raise ValueError("プロジェクトは .mfrm 拡張子で保存してください。")
                save_project(name, self.analysis)
                self.unsaved = False
                self.status.set(f"プロジェクトを保存しました：{name}（入力データを含む）")
            except Exception as exc:
                self.error(exc)

    def export_report(self):
        from tkinter import filedialog
        from mfrm_app.basic_analysis import basic_report_html
        from mfrm_app.project import write_bytes_atomic
        name = filedialog.asksaveasfilename(parent=self.root, defaultextension=".html", filetypes=[("HTML", "*.html")])
        if name:
            try:
                if Path(name).suffix.lower() != ".html":
                    raise ValueError("レポートは .html 拡張子で保存してください。")
                statistic = "Infit" if self.statistic.get() == "Infit MnSq" else "InfitZSTD"
                facet = self.map_subset.get() if self.map_subset.get() != "すべて" else None
                write_bytes_atomic(name, basic_report_html(self.analysis, pathway_statistic=statistic, pathway_facet=facet).encode("utf-8"))
                self.status.set(f"レポートを保存しました：{name}")
            except Exception as exc:
                self.error(exc)

    def close(self):
        if self.discard():
            self.worker.close()
            self.root.destroy()


def main():
    import argparse
    parser = argparse.ArgumentParser(description="MFRM Desktop prototype")
    parser.add_argument("--smoke-test", type=Path, help="Run the worker/save smoke check in this output directory")
    parser.add_argument("--check-project", type=Path, help="Restore without any estimation (requires --output-dir)")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.smoke_test:
        smoke_check(args.smoke_test)
        return
    if args.check_project:
        if not args.output_dir:
            parser.error("--check-project requires --output-dir")
        smoke_reopen(args.check_project, args.output_dir)
        return
    import tkinter as tk
    root = tk.Tk()
    DesktopApp(root)
    root.mainloop()


def smoke_check(output_dir):
    """Exercise the real spawned worker, cancellation, failure and restart."""
    import json
    import time
    from mfrm_app.project import save_project
    from mfrm_app.basic_analysis import basic_report_html, figure_bytes
    from mfrm_app.desktop_guidance import writing_example
    output_dir.mkdir(parents=True, exist_ok=True)
    worker = AnalysisWorker()
    example = writing_example()
    assert example.shape == (960, 5)
    data = demonstration_data()
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
        save_project(output_dir / "作文.mfrm", analysis)
        report = basic_report_html(analysis)
        assert report.count("<svg") == 3
        (output_dir / "report.html").write_text(report, encoding="utf-8")
        (output_dir / "wright.png").write_bytes(figure_bytes(analysis, "wright", format="png"))
        record = dict(analysis_id=analysis["analysis_id"], events=events, worker_result="complete",
                      optimizer_success=bool(analysis["result"]["opt"].success),
                      bundled_example_rows=len(example), report_svg_count=3)
        (output_dir / "smoke.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    finally:
        worker.close()


def smoke_reopen(path, output_dir):
    import json
    from mfrm_app import engine, project, basic_analysis
    def forbidden(*args, **kwargs):
        raise AssertionError("Reopening must not estimate or recompute diagnostics")
    engine.mfrm_estimate = engine.mfrm_diagnostics = forbidden
    analysis = project.load_project(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "report.html").write_text(basic_analysis.basic_report_html(analysis), encoding="utf-8")
    (output_dir / "reopen.json").write_text(json.dumps({"analysis_id": analysis["analysis_id"], "estimation_calls": 0}), encoding="utf-8")


if __name__ == "__main__":
    mp.freeze_support()
    main()
