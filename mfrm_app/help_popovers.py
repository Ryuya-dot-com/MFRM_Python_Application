"""Localized just-in-time help popover text."""

from __future__ import annotations

from collections.abc import Iterable, Mapping


HELP_POPOVER_REQUIRED_FIELDS = ("title", "what", "how", "watch")


HELP_POPOVER_TRANSLATIONS: dict[str, dict[str, dict[str, str]]] = {
    "ja": {
        "classical_dif": {
            "title": "古典的DIFスクリーニング",
            "what": "全体成績で対応づけたあと、外部グループ間で項目得点が異なるかを確認します。MFRMのfacetバイアス表とは別の入口です。",
            "how": "- Mantelのカイ二乗とETS標準化平均差で、対応層内の焦点群と参照群を比べます。\n- 順序ロジスティックLR検定で一様DIFと非一様DIFを見ます。\n- Holm/BH補正とpurificationは、多数項目での過剰検出を抑えるための確認です。",
            "watch": "これはスクリーニングです。フラグは内容・翻訳・採点設計の確認を促すもので、自動的な項目削除や公平性判定ではありません。",
        },
        "mml_person_sd": {
            "title": "MMLのPerson母集団SD",
            "what": "MMLが積分するPerson能力分布の標準偏差です。Person測定値の尺度を決めます。",
            "how": "- 固定SDでは、指定したSDの尺度上で測定値を読みます。\n- Free SDでは、EMエンジンがデータからSDを推定します。\n- Free SDのときは、収束パネルの推定SD、profile SE、95% CIも確認します。",
            "watch": "profile SEは他パラメータ固定の近似なのでやや楽観的です。Person分布が歪んだり二峰性の場合、SDだけでは解釈が不十分です。",
        },
        "wright_map": {
            "title": "Wright Map",
            "what": "Personとfacet要素を同じlogit軸に置いた図です。上ほどPersonは高能力、項目・課題は高難度として読みます。",
            "how": "- 良いターゲティングでは、項目や課題がPerson範囲を広く覆います。\n- 項目がPerson平均より上に偏ると床効果、下に偏ると天井効果の疑いがあります。\n- 狭い塊は識別範囲の狭さを示します。",
            "watch": "ラベルの密集だけで結論にしないでください。測定値表、SE、適合度、カテゴリ情報と合わせて確認します。",
        },
        "category_probability": {
            "title": "カテゴリ確率曲線",
            "what": "traitの各位置で、どの評定カテゴリが選ばれやすいかを示します。",
            "how": "- 各カテゴリには、最も選ばれやすいtrait範囲があるのが理想です。\n- 曲線のピークはカテゴリ番号順に並ぶべきです。\n- 支配領域のないカテゴリは、統合候補として確認します。",
            "watch": "カテゴリを統合する前に、実測カウント、閾値、業務上の意味を必ず確認します。",
        },
        "pathway_map": {
            "title": "Pathway Map",
            "what": "facet要素の測定値と近似Infit ZSTDを同時に見る図です。",
            "how": "- |ZSTD| < 2なら概ね許容範囲です。\n- 正のZSTDは予測より変動が大きい可能性を示します。\n- 負のZSTDは過度に予測可能、または依存を示すことがあります。\n- 同じ測定値付近で外れが集まる場合は構造的な問題を疑います。",
            "watch": "ZSTDはアプリ側の近似です。FACETS固有の自由度処理を完全再現するものではありません。",
        },
        "scree": {
            "title": "残差PCAのScree plot",
            "what": "主要なRasch次元を取り除いた後の残差相関行列の固有値を示します。",
            "how": "- 第1固有値が2未満なら強い二次元信号は弱いです。\n- 2から3程度は軽い残差次元の可能性があります。\n- 3以上は内容クラスターや設計要因を重点確認します。",
            "watch": "残差PCAは探索的です。負荷パターン、安定性、内容的意味を確認するまで次元性の結論にしません。",
        },
        "fit_scatter": {
            "title": "Infit / Outfit散布図",
            "what": "各要素のInfit MnSqとOutfit MnSqを同時に見ます。Infitは中心的なずれ、Outfitは外れ反応に敏感です。",
            "how": "- 両方が0.5から1.5なら概ね許容範囲です。\n- Outfitだけ高い場合は一部の外れ反応を確認します。\n- 両方が高い場合は項目・課題・評定者訓練を見直します。",
            "watch": "MnSqは1.0が期待値です。0を中心に読む指標ではありません。",
        },
        "bias_heatmap": {
            "title": "バイアス / 交互作用ヒートマップ",
            "what": "2つのfacet水準の組み合わせについて、主効果を取り除いた後の局所的なずれを示します。",
            "how": "- 0に近いセルは加法的期待に近いです。\n- |t|が2以上のセルは重点確認します。\n- セル数が少ないとtが不安定になります。\n- 多数セルでは多重補正後の結果も見ます。",
            "watch": "1つのfacetペアだけで「バイアスなし」とは言えません。全ペアの画面と内容上の根拠を合わせます。",
        },
        "forest_measures": {
            "title": "測定値のForest plot",
            "what": "各要素の推定値、50% CI、95% CIをlogit尺度で表示します。",
            "how": "- 95% CIが明確に離れると差の根拠が強くなります。\n- 50% CIが重なる場合は過度に順位づけしません。\n- CIが広い要素は観測数や設計を確認します。",
            "watch": "区間は近似SEに基づきます。高リスク判断では信頼性、適合、設計根拠を併用してください。",
        },
        "qq_residuals": {
            "title": "標準化残差のQ-Q plot",
            "what": "標準化残差を理論正規分布の分位点と比較します。",
            "how": "- 点が直線付近なら残差分布は概ね想定に近いです。\n- 端で大きく外れる場合は外れ反応や局所的不適合を確認します。\n- S字形は残差構造の可能性があります。",
            "watch": "単独では診断を確定しません。PCA、Pathway Map、fit表と合わせて読みます。",
        },
        "ecdf_measures": {
            "title": "測定値ECDF",
            "what": "Personやfacet要素の測定値分布を累積分布で示します。",
            "how": "- 平坦な区間は尺度上の空白を示します。\n- 急な立ち上がりは測定値の集中を示します。\n- Personと項目・課題の曲線を比べてターゲティングを見ます。",
            "watch": "Wright Mapの補助です。ビン幅に左右されず、狭いカバレッジ不足を見つけやすくします。",
        },
        "posterior_trace": {
            "title": "Posterior trace plot",
            "what": "選択したパラメータについて、各チェーンのサンプリング軌跡を表示します。",
            "how": "- チェーンが重なり、横ばいなら混合は良好です。\n- チェーンが離れている場合はRhatを確認します。\n- スパイクや途切れはdivergenceやモデル設定の再確認につなげます。",
            "watch": "目視はRhat/ESSの補助です。収束診断の代替ではありません。",
        },
        "posterior_rhat_ess": {
            "title": "Rhat / ESSバー",
            "what": "選択したパラメータごとのRhatと有効サンプルサイズを基準線付きで示します。",
            "how": "- Rhat < 1.01ならチェーン混合は概ね良好です。\n- ESSが小さいパラメータは推定精度に注意します。\n- 赤いバーや基準超過を最初に確認します。",
            "watch": "多数のパラメータで問題が出る場合は、反復数、adapt_delta、モデル指定を見直してから引用します。",
        },
        "coverage_heatmap": {
            "title": "観測カバレッジヒートマップ",
            "what": "Personとfacet水準の組み合わせに観測があるかを示します。",
            "how": "- 観測済みセルが多いほど設計は安定しやすいです。\n- 大きな空白ブロックはスパース性や非連結の疑いです。\n- 特定水準に空白が集中する場合はサンプリングやanchorを確認します。",
            "watch": "表示範囲で空セルが多い場合は、測定値を読む前に連結性レビューを優先します。",
        },
        "category_usage": {
            "title": "カテゴリ使用頻度",
            "what": "各評定カテゴリが実際にどれだけ使われたかを表示します。",
            "how": "- すべてのカテゴリに十分な観測があるかを見ます。\n- 0または極端に少ないカテゴリは統合候補です。\n- facet水準ごとの偏りは評定者バイアスや内容効果の入口です。",
            "watch": "低頻度カテゴリは下流の閾値 disorder につながりやすいです。",
        },
        "misfit_ranking": {
            "title": "Misfitランキング",
            "what": "近似|ZSTD|が大きいfacet要素を優先順に並べます。",
            "how": "- 上位5から10件を最初の確認対象にします。\n- 正のZSTDは反応が予測よりばらつく可能性を示します。\n- 負のZSTDは過度に予測可能な反応を示します。",
            "watch": "観測数が少ない行ではZSTDが不安定です。Pathway Mapと併用して測定値との関係を見ます。",
        },
        "facet_distribution": {
            "title": "Facet要素分布",
            "what": "facet内の要素測定値の分布を箱ひげと点で表示します。",
            "how": "- 広がりは、そのfacetがPerson順位にどれだけ影響するかの手がかりです。\n- 狭い分布は識別力の低さを示すことがあります。\n- 外れ点は分割、anchor、内容確認の候補です。",
            "watch": "信頼性表と合わせて見ます。広がりが大きく信頼性も高いfacetは情報量が高い傾向があります。",
        },
        "rater_agreement": {
            "title": "評定者一致",
            "what": "同じPersonと項目・課題に対する評定者間の一致率、相関、Krippendorff alphaを確認します。",
            "how": "- 完全一致率と隣接一致率をまず見ます。\n- alphaは偶然一致を考慮した一貫性の目安です。\n- 低い一致でもMFRM適合が良い場合、severity差をモデルが吸収している可能性があります。",
            "watch": "MFRMは高い生一致を必須としません。低一致と不適合が同時に出る場合が重要な警告です。",
        },
        "threshold_map": {
            "title": "閾値マップ",
            "what": "評定カテゴリ間のstep閾値が順序どおり並ぶかを確認します。",
            "how": "- 閾値はカテゴリ順に単調増加するのが理想です。\n- 後のstepが低い場合はカテゴリ使用の混乱を疑います。\n- 隣接距離が極端に狭い、または広い場合はレビューします。",
            "watch": "disorderだけで即統合とはしません。カテゴリ確率曲線と実測カウントを確認してから判断します。",
        },
        "zstd_distribution": {
            "title": "ZSTD分布",
            "what": "全facet要素の近似Infit/Outfit ZSTDの分布を示します。",
            "how": "- 大部分が±2内なら概ね安定しています。\n- 右裾が重い場合はノイズの多い要素が多い可能性があります。\n- 左裾が重い場合は過度な予測可能性や範囲制限を疑います。",
            "watch": "要素数が多いと一部の|ZSTD|>2は偶然でも起こります。個数だけでなく裾の大きさを見ます。",
        },
    },
}


def localized_help_popover_topic(
    topic_key: str,
    *,
    lang: str,
    fallback: Mapping[str, str],
) -> dict[str, str]:
    """Overlay localized text onto a fallback popover topic."""
    topic = {field: str(fallback.get(field, "")) for field in HELP_POPOVER_REQUIRED_FIELDS}
    localized = HELP_POPOVER_TRANSLATIONS.get(str(lang), {}).get(str(topic_key), {})
    for field in HELP_POPOVER_REQUIRED_FIELDS:
        value = str(localized.get(field, "")).strip()
        if value:
            topic[field] = value
    return topic


def missing_translation_fields(topic_keys: Iterable[str], *, lang: str = "ja") -> list[tuple[str, str]]:
    """Return missing localized popover fields for test and release checks."""
    localized_topics = HELP_POPOVER_TRANSLATIONS.get(str(lang), {})
    missing: list[tuple[str, str]] = []
    for topic_key in topic_keys:
        topic = localized_topics.get(str(topic_key), {})
        for field in HELP_POPOVER_REQUIRED_FIELDS:
            if not str(topic.get(field, "")).strip():
                missing.append((str(topic_key), field))
    return missing
