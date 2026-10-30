# MFRM Desktop for macOS (Apple Silicon) — START HERE / はじめに

This desktop preview is for macOS on Apple Silicon (`arm64`: M1 / M2 / M3 / M4).
Intel Macs have not been verified. Python is bundled inside `MFRM-Desktop.app`;
installing Python, R, FACETS, ConQuest or Xcode tools is not required.

---

## First launch on macOS Sequoia 15+ (Gatekeeper notice)

MFRM Desktop is distributed free of charge as open-source academic software
under the MIT License and does **not** use paid Apple Developer Program
notarization. On first launch, macOS displays:

> **"MFRM-Desktop" Not Opened** — Apple could not verify "MFRM-Desktop" is free of malware...

1. Click **Done** (grey button). **Do NOT click the blue "Move to Trash" button.**
2. Open **System Settings → Privacy & Security**, scroll down to **Security**,
   and click **Open Anyway** next to `"MFRM-Desktop" was blocked to protect your Mac`.
3. Confirm **Open Anyway** and authenticate with Touch ID or your Mac password.
   Subsequent launches open normally.
4. Optional integrity check before opening (Terminal):
   `shasum -a 256 MFRM-Desktop-macOS-arm64.zip` — compare with `SHA256SUMS.txt`
   on the GitHub Release page. Analysis runs 100% offline; you can turn Wi-Fi
   off before opening confidential rating data.

### 【日本語】macOS Sequoia 15以降での初回起動手順（重要）

本アプリは無償公開の学術オープンソース（MITライセンス）であり、有償の Apple Developer Program 公証は付与していません。初回起動時に *「“MFRM-Desktop” は開いていません — マルウェアが含まれていないことを検証できませんでした」* と表示された場合は、以下の手順で起動してください（2回目以降は通常どおり開きます）。

1. **青い「ゴミ箱に入れる」は押さず、下の灰色の「完了」をクリック**します。
2. **「システム設定」→「プライバシーとセキュリティ」** を開き、画面下部の「セキュリティ」欄にある **「このまま開く」** をクリックして Touch ID／パスワードで承認します。
3. 必要に応じて起動前にターミナルで `shasum -a 256 MFRM-Desktop-macOS-arm64.zip` を実行し、GitHub Release 記載の SHA-256 と一致することを確認できます。計算・保存はすべて完全オフラインで動作するため、機密データを扱う際は Wi-Fi を切って実行することも可能です。

---

## Quick trial (3 clicks)

- **MFRM (3-facet sample):** Select **Data → Sample → Next → Analyze** to fit
  synthetic writing ratings (`0–5`). Inspect **Figure** and **Table**; use
  **Reading guide** and table notes to interpret outputs.
- **Univariate G/D-Study (2-facet sample):** Select **Data → G/D 2-facet sample → Next → Analyze**
  to run a balanced `Person × Rater × Task` generalizability and decision study
  immediately. For 3+ facet datasets, use **Fixed-facet handling** on the
  **G/D-Study plan** tab to condition on a single level or compute a verified
  balanced sum/mean composite.
- **Language switch:** Select **日本語** or **English** in the top-right selector
  at any time without refitting.

## Save projects, figures and reports

- **Save project (`.mfrm`):** Stores input data, identifiers, settings and
  fitted results in a single file that reopens via **Data → Open project** (or
  Finder open) without re-estimation.
- **Save report (`.html`):** Exports a self-contained offline HTML report with
  all tables and vector SVG figures.
- **Save figure (`.png` / `.svg`):** Switch to **Manuscript** mode to move
  titles and explanatory notes into an accompanying `.caption.txt` file.

## Scope and interpretation reminders

- RSM and PCM support JMLE and fixed-population-SD MML. Outputs include Wright
  and Pathway maps, category probabilities, conditional bias diagnostics,
  specified intervals and joint interaction analysis (fixed `N(0,1)` MML).
- Conditional bias diagnostics hold main-effect calibration fixed and do not
  establish significance or fairness. Dedicated intervals and joint interactions
  check target-specific numerical conditions; pointwise approximations do not
  guarantee finite-sample coverage.
- G/D-Study uses original observed scores under complete balanced crossing with
  1–2 random facets, unit weights and one observation per cell. Multivariate
  G-theory and G/Phi intervals are not implemented.

Online guide and source code:
https://ryuya-dot-com.github.io/MFRM_Python_Application/
https://github.com/Ryuya-dot-com/MFRM_Python_Application
