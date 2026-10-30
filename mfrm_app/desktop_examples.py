"""Curated synthetic example datasets and presets for MFRM-Desktop.

Every dataset is strictly synthetic, deterministic, and accompanied by clear
design metadata (facet structure, missingness/sparsity, anchor presets, and
what diagnostic behavior to inspect).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mfrm_app import engine

EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"


def _sample_rsm_score(rng: np.random.Generator, eta: float, steps: np.ndarray, rating_min: int) -> int:
    cum = np.r_[0.0, np.cumsum(np.asarray(steps, dtype=float))]
    probs = engine.category_prob_rsm(np.array([float(eta)]), cum)[0]
    probs = np.clip(probs, 0.0, 1.0)
    probs = probs / probs.sum()
    return int(rating_min + rng.choice(len(probs), p=probs))


def generate_clinical_4facet_df() -> pd.DataFrame:
    """24 Examinees x 3 Raters x 2 Stations x 2 Competencies x 2 Sessions = 576 rows (1..4)."""
    rng = np.random.default_rng(2026103101)
    examinees = [f"M{i:02d}" for i in range(1, 25)]
    raters = [("Dr_A", -0.35), ("Dr_B", 0.05), ("Dr_C", 0.30)]
    stations = [("History", -0.25), ("Counseling", 0.25)]
    competencies = [("Communication", -0.20), ("Reasoning", 0.20)]
    sessions = [("AM", -0.10), ("PM", 0.10)]
    steps = np.array([-1.15, 0.0, 1.15])
    theta = np.linspace(-1.45, 1.45, len(examinees)) + rng.normal(0.0, 0.12, len(examinees))

    rows = []
    for pi, ex in enumerate(examinees):
        for r_name, r_sev in raters:
            for st_name, st_diff in stations:
                for comp_name, comp_diff in competencies:
                    for ses_name, ses_eff in sessions:
                        eta = theta[pi] - r_sev - st_diff - comp_diff - ses_eff
                        score = _sample_rsm_score(rng, eta, steps, rating_min=1)
                        rows.append((ex, r_name, st_name, comp_name, ses_name, str(score)))
    return pd.DataFrame(
        rows,
        columns=["Examinee", "Rater", "Station", "Competency", "Session", "Rating"],
    )


def generate_anchor_equating_df() -> pd.DataFrame:
    """Two cohorts (Spring S01..S20, Fall F01..F20) linked via bank rater/task anchors (600 rows, 0..4)."""
    rng = np.random.default_rng(2026103102)
    rater_true = {"R1": -0.50, "R2": 0.35, "R3": -0.15, "R4": 0.25}
    task_true = {"Task1": -0.25, "Task2": 0.25}
    crit_true = {"Content": -0.20, "Organization": 0.05, "Language": 0.15}
    steps = np.array([-1.35, -0.40, 0.40, 1.35])

    rows = []
    # Spring cohort: S01..S20 (mean 0.0), rated by R1 and R2
    theta_spring = np.linspace(-1.3, 1.3, 20)
    for pi, th in enumerate(theta_spring, start=1):
        person = f"S{pi:02d}"
        for r in ("R1", "R2"):
            for t in ("Task1", "Task2"):
                for c in ("Content", "Organization", "Language"):
                    eta = th - rater_true[r] - task_true[t] - crit_true[c]
                    rows.append((person, r, t, c, str(_sample_rsm_score(rng, eta, steps, 0))))

    # Fall cohort: F01..F20 (mean +0.55), rated by R2, R3, R4 (2-3 raters/person = 360 rows)
    theta_fall = np.linspace(-0.75, 1.85, 20)
    fall_rater_sets = [("R2", "R3", "R4")] * 10 + [("R2", "R3"), ("R2", "R4")] * 5
    for pi, (th, r_set) in enumerate(zip(theta_fall, fall_rater_sets), start=1):
        person = f"F{pi:02d}"
        for r in r_set:
            for t in ("Task1", "Task2"):
                for c in ("Content", "Organization", "Language"):
                    eta = th - rater_true[r] - task_true[t] - crit_true[c]
                    rows.append((person, r, t, c, str(_sample_rsm_score(rng, eta, steps, 0))))

    return pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Criterion", "Score"])


ANCHOR_EQUATING_PRESETS = {
    "anchors": [
        {"Facet": "Rater", "Level": "R1", "Anchor": -0.50},
        {"Facet": "Rater", "Level": "R2", "Anchor": 0.35},
        {"Facet": "Task", "Level": "Task1", "Anchor": -0.25},
        {"Facet": "Task", "Level": "Task2", "Anchor": 0.25},
    ],
    "group_anchors": [
        {"Facet": "Criterion", "Level": "Content", "Group": "Rubric", "GroupValue": 0.0},
        {"Facet": "Criterion", "Level": "Organization", "Group": "Rubric", "GroupValue": 0.0},
        {"Facet": "Criterion", "Level": "Language", "Group": "Rubric", "GroupValue": 0.0},
    ],
}


def generate_missing_ratings_df() -> pd.DataFrame:
    """45 Persons x 4 Raters x 2 Tasks x 3 Criteria with ~17% absent rows + 12 blank Score cells."""
    rng = np.random.default_rng(2026103103)
    persons = [f"P{i:02d}" for i in range(1, 46)]
    raters = [("R1", -0.45), ("R2", -0.10), ("R3", 0.15), ("R4", 0.40)]
    tasks = [("EssayA", -0.20), ("EssayB", 0.20)]
    criteria = [("Ideas", -0.25), ("Coherence", 0.05), ("Mechanics", 0.20)]
    steps = np.array([-1.40, -0.45, 0.45, 1.40])
    theta = np.linspace(-1.55, 1.55, len(persons))

    full_rows = []
    for pi, p in enumerate(persons):
        person_rows = []
        for r_name, r_sev in raters:
            for t_name, t_diff in tasks:
                for c_name, c_diff in criteria:
                    eta = theta[pi] - r_sev - t_diff - c_diff
                    sc = _sample_rsm_score(rng, eta, steps, 0)
                    person_rows.append([p, r_name, t_name, c_name, str(sc)])
        # Retain ~82% of cells per person while guaranteeing >=12 ratings/person
        keep_idx = sorted(rng.choice(len(person_rows), size=20, replace=False))
        for idx in keep_idx:
            full_rows.append(person_rows[idx])

    # Mark 12 deterministic rows with explicit blank Score ("") to illustrate CSV blank handling
    blank_positions = np.linspace(5, len(full_rows) - 10, 12, dtype=int)
    for pos in blank_positions:
        full_rows[int(pos)][4] = ""
    return pd.DataFrame(full_rows, columns=["Person", "Rater", "Task", "Criterion", "Score"])


def generate_sparse_spiral_df() -> pd.DataFrame:
    """Connected sparse spiral design: 64 Persons x 8 Raters (2 raters/person = 25% density) x 2 Tasks x 3 Criteria."""
    rng = np.random.default_rng(2026103104)
    persons = [f"E{i:03d}" for i in range(1, 65)]
    raters = [f"R{i}" for i in range(1, 9)]
    rater_sev = np.linspace(-0.65, 0.65, len(raters))
    tasks = [("Prompt1", -0.25), ("Prompt2", 0.25)]
    criteria = [("Fluency", -0.25), ("Accuracy", 0.10), ("Interaction", 0.15)]
    steps = np.array([-1.60, -0.80, 0.0, 0.80, 1.60])
    theta = np.linspace(-1.65, 1.65, len(persons)) + rng.normal(0.0, 0.10, len(persons))

    rows = []
    for pi, p in enumerate(persons):
        r_idx_a = pi % len(raters)
        r_idx_b = (pi + 1) % len(raters)
        for ri in (r_idx_a, r_idx_b):
            r_name = raters[ri]
            for t_name, t_diff in tasks:
                for c_name, c_diff in criteria:
                    eta = theta[pi] - rater_sev[ri] - t_diff - c_diff
                    sc = _sample_rsm_score(rng, eta, steps, 0)
                    rows.append((p, r_name, t_name, c_name, str(sc)))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Criterion", "Score"])


def generate_disconnected_shifts_df() -> pd.DataFrame:
    """Two disjoint rater shifts (M01..M16 x R1,R2 vs A01..A16 x R3,R4) sharing Task1,Task2 (128 rows)."""
    rng = np.random.default_rng(2026103105)
    steps = np.array([-1.30, -0.35, 0.35, 1.30])
    rater_sev = {"R1": -0.35, "R2": 0.25, "R3": -0.20, "R4": 0.40}
    task_diff = {"Task1": -0.25, "Task2": 0.25}

    rows = []
    for shift_prefix, rater_pair, mu in (("M", ("R1", "R2"), -0.15), ("A", ("R3", "R4"), 0.25)):
        thetas = np.linspace(-1.2, 1.2, 16) + mu
        for i, th in enumerate(thetas, start=1):
            p = f"{shift_prefix}{i:02d}"
            for r in rater_pair:
                for t in ("Task1", "Task2"):
                    eta = th - rater_sev[r] - task_diff[t]
                    sc = _sample_rsm_score(rng, eta, steps, 0)
                    rows.append((p, r, t, str(sc)))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Score"])


def generate_wide_spreadsheet_df() -> pd.DataFrame:
    """Wide spreadsheet layout: 20 Persons x 2 Tasks with R1..R4 as score columns (1 blank per row)."""
    rng = np.random.default_rng(2026103106)
    steps = np.array([-1.50, -0.70, 0.0, 0.70, 1.50])
    rater_sev = [-0.45, -0.10, 0.15, 0.40]
    task_diff = {"Presentation": -0.20, "Q_and_A": 0.20}
    thetas = np.linspace(-1.4, 1.4, 20)

    rows = []
    for pi, th in enumerate(thetas, start=1):
        person = f"S{pi:02d}"
        for ti, (t_name, t_d) in enumerate(task_diff.items()):
            omitted_rater = (pi + ti) % 4
            row_scores = []
            for ri in range(4):
                if ri == omitted_rater:
                    row_scores.append("")
                else:
                    eta = th - rater_sev[ri] - t_d
                    row_scores.append(str(_sample_rsm_score(rng, eta, steps, 0)))
            rows.append((person, t_name, *row_scores))
    return pd.DataFrame(rows, columns=["Person", "Task", "R1", "R2", "R3", "R4"])


SCENARIO_CATALOG = [
    {
        "key": "writing_3facet",
        "filename": "desktop_writing.csv",
        "badge_ja": "標準 3相",
        "badge_en": "Standard 3-Facet",
        "title_ja": "作文評定（3相・完全交差・960行）",
        "title_en": "Essay Writing (3 Facets, Crossed, 960 rows)",
        "summary_ja": "30人 × 評価者4 × 課題2 × 観点4（0〜5点）。MFRM（RSM/PCM・JMLE/MML）・バイアス・指定区間・3相G/D条件付けの基本確認用。",
        "summary_en": "30 Persons × 4 Raters × 2 Tasks × 4 Criteria (0–5). Baseline dataset for RSM/PCM, JMLE/MML, bias, intervals, and conditioned G/D.",
        "workflow": "MFRM",
        "person_col": "Person",
        "score_col": "Score",
        "facet_cols": ["Rater", "Task", "Criterion"],
        "rating_min": 0,
        "rating_max": 5,
    },
    {
        "key": "clinical_4facet",
        "filename": "clinical_osce_4facet.csv",
        "badge_ja": "複数相（4相）",
        "badge_en": "Multi-Facet (4)",
        "title_ja": "臨床OSCE評価（4相・1〜4点・576行）",
        "title_en": "Clinical OSCE (4 Facets, 1–4 Rubric, 576 rows)",
        "summary_ja": "受験者24 × 評価医3 × ステーション2 × コンピテンシー2 × 実施枠2（最低点1〜最高点4）。4つの構造相を同時推定する例。",
        "summary_en": "24 Examinees × 3 Raters × 2 Stations × 2 Competencies × 2 Sessions (1–4 scale). Illustrates simultaneous 4-facet calibration.",
        "workflow": "MFRM",
        "person_col": "Examinee",
        "score_col": "Rating",
        "facet_cols": ["Rater", "Station", "Competency", "Session"],
        "rating_min": 1,
        "rating_max": 4,
    },
    {
        "key": "anchor_equating",
        "filename": "anchor_equating_cohorts.csv",
        "badge_ja": "アンカー等化",
        "badge_en": "Anchor Equating",
        "title_ja": "2コホート尺度連結（アンカー自動設定・540行）",
        "title_en": "2-Cohort Scale Equating (Preset Anchors, 540 rows)",
        "summary_ja": "春期(S)・秋期(F)の計40人。既存バンクの評価者(R1,R2)・課題(Task1,Task2)の固定アンカーと観点のグループアンカー(平均0)を自動設定します。",
        "summary_en": "40 examinees across Spring/Fall cohorts. Pre-loads fixed rater/task bank anchors (R1, R2, Task1, Task2) and a Criterion group anchor (mean 0).",
        "workflow": "MFRM",
        "person_col": "Person",
        "score_col": "Score",
        "facet_cols": ["Rater", "Task", "Criterion"],
        "rating_min": 0,
        "rating_max": 4,
        "anchors": ANCHOR_EQUATING_PRESETS["anchors"],
        "group_anchors": ANCHOR_EQUATING_PRESETS["group_anchors"],
    },
    {
        "key": "missing_ratings",
        "filename": "ratings_with_missing.csv",
        "badge_ja": "評定欠損あり",
        "badge_en": "Missing Ratings",
        "title_ja": "欠損・空欄セル混在データ（900行・空欄12行）",
        "title_en": "Ratings with Missing Rows & Blank Cells (900 rows)",
        "summary_ja": "45人 × 評価者4 × 課題2 × 観点3から約17%の未採点行を除き、さらに空欄得点セル12行を含む実務想定データ。除外確認と不完全設計の推定を体験できます。",
        "summary_en": "45 Persons × 4 Raters × 2 Tasks × 3 Criteria with ~17% absent rows plus 12 blank Score cells. Demonstrates missing-row audit and incomplete-data fit.",
        "workflow": "MFRM",
        "person_col": "Person",
        "score_col": "Score",
        "facet_cols": ["Rater", "Task", "Criterion"],
        "rating_min": 0,
        "rating_max": 4,
    },
    {
        "key": "sparse_spiral",
        "filename": "sparse_rater_spiral.csv",
        "badge_ja": "スパース（連結）",
        "badge_en": "Sparse (Linked)",
        "title_ja": "評価者スパイラル割当（充填率25%・768行）",
        "title_en": "Spiral Rater Assignment (25% Density, 768 rows)",
        "summary_ja": "評価者8名のうち各受験者(64名)を2名ずつ循環ペア(R1-R2, R2-R3, …, R8-R1)で採点。75%が未観測でも単一尺度へ完全連結される実務設計例。",
        "summary_en": "64 examinees rated by overlapping pairs from an 8-rater pool (75% unobserved cells). Demonstrates a connected, identified sparse judging plan.",
        "workflow": "MFRM",
        "person_col": "Person",
        "score_col": "Score",
        "facet_cols": ["Rater", "Task", "Criterion"],
        "rating_min": 0,
        "rating_max": 5,
    },
    {
        "key": "sparse_disconnected",
        "filename": "disconnected_shifts.csv",
        "badge_ja": "未接続警告例",
        "badge_en": "Disconnected Check",
        "title_ja": "午前・午後シフト分断例（未接続サブセット診断・128行）",
        "title_en": "Unlinked AM/PM Shifts (Disconnected Subset Warning, 128 rows)",
        "summary_ja": "午前班(M×R1,R2)と午後班(A×R3,R4)が共通課題を受けたが共通評価者がいない例。設定タブで「未接続サブセット（階数欠損1）」警告とアンカーによる解消を確認できます。",
        "summary_en": "AM shift (R1,R2) and PM shift (R3,R4) share tasks but zero raters. Triggers the disconnected-subset warning (nullity = 1) and lets you test anchor linking.",
        "workflow": "MFRM",
        "person_col": "Person",
        "score_col": "Score",
        "facet_cols": ["Rater", "Task"],
        "rating_min": 0,
        "rating_max": 4,
    },
    {
        "key": "wide_spreadsheet",
        "filename": "wide_rubric_sheet.csv",
        "badge_ja": "横長CSV変換",
        "badge_en": "Wide → Long",
        "title_ja": "横長Excel形式シート（評価者列 R1〜R4 を縦長変換）",
        "title_en": "Wide Spreadsheet Layout (Convert R1–R4 Columns to Long)",
        "summary_ja": "Person, Task の横に R1〜R4 の得点列が並ぶ横長データ（未採点空欄あり）。読込時に「横長→縦長変換プレビュー」を開いてそのまま変換できます。",
        "summary_en": "Wide table with R1–R4 as separate columns and blank unrated cells. Opens the Wide-to-Long conversion preview so you can melt rater columns in one click.",
        "workflow": "MFRM",
        "is_wide": True,
        "wide_defaults": {
            "id_cols": ["Person", "Task"],
            "score_cols": ["R1", "R2", "R3", "R4"],
            "new_facet_name": "Rater",
            "score_col_name": "Score",
        },
        "person_col": "Person",
        "score_col": "Score",
        "facet_cols": ["Rater", "Task"],
        "rating_min": 0,
        "rating_max": 5,
    },
    {
        "key": "gd_2facet",
        "filename": "gd_two_facet.csv",
        "badge_ja": "G/D-Study 2相",
        "badge_en": "G/D-Study 2-Facet",
        "title_ja": "G/D-Study 2相完全均衡サンプル（240行）",
        "title_en": "G/D-Study Balanced 2-Facet Sample (240 rows)",
        "summary_ja": "30人 × 評価者4 × 課題2（0〜5点・欠測なし）。観測得点の分散成分分解とD-Study（評価者数×課題数の信頼性予測）へ自動で切り替わります。",
        "summary_en": "30 Persons × 4 Raters × 2 Tasks (0–5, fully crossed). Automatically switches to standalone G/D-Study for variance components and D-Study projections.",
        "workflow": "G/D-Study",
        "person_col": "Person",
        "score_col": "Score",
        "facet_cols": ["Rater", "Task"],
        "rating_min": 0,
        "rating_max": 5,
    },
]


def get_scenario(key: str) -> dict:
    for item in SCENARIO_CATALOG:
        if item["key"] == key:
            return dict(item)
    raise ValueError(f"Unknown example scenario: {key}")


def load_scenario_dataframe(key: str) -> pd.DataFrame:
    from mfrm_app.desktop_guidance import gd_two_facet_example, read_desktop_csv, writing_example

    item = get_scenario(key)
    csv_path = EXAMPLES_DIR / item["filename"]
    if csv_path.exists():
        return read_desktop_csv(csv_path)
    generators = {
        "writing_3facet": writing_example,
        "clinical_4facet": generate_clinical_4facet_df,
        "anchor_equating": generate_anchor_equating_df,
        "missing_ratings": generate_missing_ratings_df,
        "sparse_spiral": generate_sparse_spiral_df,
        "sparse_disconnected": generate_disconnected_shifts_df,
        "wide_spreadsheet": generate_wide_spreadsheet_df,
        "gd_2facet": gd_two_facet_example,
    }
    return generators[key]()
