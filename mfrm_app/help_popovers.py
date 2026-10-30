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
            "how": "- 固定SDでは、Person分布のSDを指定します。\n- 自由SDでは、EMエンジンがデータからSDを推定します。\n- 推定SDは技術的な点推定値として確認します。SE・信頼区間は検証が完了するまで保留します。\n- 他のエンジンとの比較には、モデル・同定制約・尺度をそろえる必要があります。",
            "watch": "他パラメータを固定した曲率は、プロファイル情報ではありません。全パラメータの停留性と、他の推定パラメータを考慮した不確実性評価が未検証です。母集団の共変量を指定した場合、σは正規分布を仮定した条件付き残差SDです。",
        },
        "wright_map": {
            "title": "Wright Map",
            "what": "Personと各相の要素を同じlogit軸に置いた図です。上ほどPersonは高能力、項目・課題は高難度として読みます。",
            "how": "- 各相内の順序と間隔を読みます。\n- 得点分布、カテゴリー境界、利用可能な測定精度を合わせて確認します。",
            "watch": "相の原点は同定制約に依存します。位置の重なりだけで天井・床効果や妥当性を結論にせず、保留されたSEは利用しません。",
        },
        "category_probability": {
            "title": "カテゴリ確率曲線",
            "what": "traitの各位置で、どの評定カテゴリが選ばれやすいかを示します。",
            "how": "- 中間カテゴリーが最も選ばれやすい範囲を見ます。\n- 実際の度数、閾値の順序、採点基準の意味と合わせて確認します。",
            "watch": "ピークの欠如や閾値の逆転だけではカテゴリーを統合しません。読み方ガイドにLinacre（2002）の出典があります。",
        },
        "pathway_map": {
            "title": "測定値と適合度の図",
            "what": "縦軸に測定値、横軸に選択した適合度を示します。初期表示はInfit MnSqで、1がモデルの基準です。近似ZSTDにも切り替えられます。",
            "how": "同じ測定値付近でずれが集まるかを確認し、該当する評定・観測数・内容に戻って点検します。",
            "watch": "基準線は確認の目安です。範囲内であることは妥当性の証明ではなく、範囲外であることも除外の根拠にはなりません。近似したZSTDは、他ソフトと同じ値になるとは限りません。",
        },
        "scree": {
            "title": "残差の主成分分析：固有値",
            "what": "残差の相関に、どのような構造が残っているかを探す図です。",
            "how": "固有値だけでなく、負荷量、関係する要素、観測数、結果の安定性と内容上の意味を確認します。",
            "watch": "固有値の目安だけで次元数を決めたり、一次元性を証明したりはできません。",
        },
        "fit_scatter": {
            "title": "Infit・Outfitの散布図",
            "what": "各要素について、モデルの予測と評定のずれを2種類の指標で示します。",
            "how": "点を選び、元の評定、観測数、カテゴリの使用、課題内容を確認します。計算できない値を「適合している」と読まないでください。",
            "watch": "MNSQは1付近を基準に読みます。表示範囲だけで評価者の能力を判断したり、対象者・項目を除外したりしません。",
        },
        "bias_heatmap": {
            "title": "バイアス・交互作用のヒートマップ",
            "what": "要因の水準の組合せごとに、加法モデルの予測からのずれを示します。",
            "how": "セルの値、観測数、担当の割付と元の評定を確認します。観測が少ないセルでは値が不安定になることがあります。",
            "watch": "現在のセル計算は記述用です。検定・区間の解釈は保留中であり、多重比較の補正を行うだけでは推論に使えるようになりません。学校別DFFには、比較したい差を分離できる割付かどうかの確認が別に必要です。Helpで交互作用のトピックを選び、詳細を開くと、比較を計画する際の説明を読めます。",
        },
        "forest_measures": {
            "title": "要素ごとの推定値と区間",
            "what": "各要素の推定値と、対応する不確実性の方法で計算した区間を示します。",
            "how": "区間の計算法、固定した条件、保留状態を確認します。幅が広い要素は観測数や評定設計も見直してください。",
            "watch": "区間の重なりだけで差の検定にはなりません。近似区間や事後SDを、検証済みの信頼区間・標準誤差と読み替えないでください。",
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
