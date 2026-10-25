# Figure Legends & Methods（英文投稿版）

> **生成时间**：2026-10-25　|　**适用**：杨泽镔《关键单体香气化合物感知特性》、华铖凯《关键混合香气化合物感知特性》两篇论文的**重分析图**

---

## 0. 使用须知（先读这一节）

1. 以下英文**不是论文原文的翻译**，而是对我们**重新计算并重新作图**的结果所写的图注与方法学段落。
2. 每一个数字都在第 3 节列出了**出处文件 / 脚本 / 验证状态**。标 🟢 的可直接用，标 🟠 的使用前必须复核。
3. 图 1 的剂量响应参数我们做过**独立复现**：模型确认为 **4 参数 log-logistic（含自由下渐近线）**，R² 15/15 完全一致 —— 见 3.1。
4. **不要**把第 1 节图注和第 2 节 Methods 里的数字与论文正文里对不上的数字混用。冲突处已在第 4 节列出。
5. 所有文件在 `教学示范图\` 目录下，均为新增文件，原始 Excel 未做任何修改。

---

# 1. Figure Legends（英文）

## Figure 1 — Dose–response parameter space of 15 key aroma compounds

> 文件：`fig_bubble.png`

**Figure 1. Potency–efficacy map of the 15 single aroma compounds.** Each circle is one compound, positioned by its half-maximal effective concentration (EC₅₀, x-axis, log₁₀ scale, ppb in model solution) and its fitted upper asymptote Iₘₐₓ (y-axis, perceived odour intensity on a 0–12 scale). Circle area is proportional to the Hill coefficient *n* (steepness of the dose–response transition), and circle colour encodes the coefficient of determination *R*² of the four-parameter log-logistic fit. The horizontal red line at Iₘₐₓ = 12 marks the physical ceiling of the 0–12 intensity scale; the pink band above it is the **extrapolation region**, in which a fitted Iₘₐₓ cannot be realised by a panellist and therefore has no direct sensory meaning. Compound labels are offset with a greedy collision-avoidance algorithm; dashed leaders indicate the true data position where a label has been displaced.

**Read-off for the Discussion:** EC₅₀ values span approximately **8.9 ppb to 8.0 × 10⁷ ppb** (≈9 × 10⁶-fold), i.e. the panel discriminates compounds across seven orders of magnitude in potency. Five compounds carry an extrapolated Iₘₐₓ above the scale ceiling — **phenylethanol (19.83), ethyl phenylacetate (13.87), *trans,cis*-2,6-nonadienal (12.72), dimethyl trisulfide (12.71) and ethyl isobutyrate (12.70)** — and must be reported as ">12, not estimable" rather than as point estimates. Butyric acid was fitted from only four concentration levels and returned a saturated fit (*R*² = 1.000 for a four-parameter model on four points), so its parameter set carries no residual degrees of freedom.

**Caption footnote (recommended):** *Parameter estimates were obtained by non-linear least-squares fitting of a four-parameter log-logistic function (see Methods 2.3). Confidence intervals are reported in Supplementary Table S1.*

---

## Figure 2 — Additivity of binary mixtures

> 文件：`fig_additivity.png`

**Figure 2. Test of perceptual additivity in binary mixtures. (A)** Observed CATA (check-all-that-apply) response frequency for each binary mixture is plotted against the additive prediction

  *p*₀ + (*p*_A − *p*₀) + (*p*_B − *p*₀),

clipped to [0, 1], where *p*₀ is the response frequency of the all-low baseline (母液1) and *p*_A, *p*_B are the frequencies of the corresponding single-compound-high samples. Points are coloured by aroma attribute; the solid grey line is perfect additivity (*y* = *x*); the stepped line is the binned trend (deciles of the predicted value, bins with *n* ≥ 5 shown). The inset reports the Pearson correlation, root-mean-square error (RMSE) and mean signed deviation.

**(B)** Mean signed deviation (observed − predicted) per aroma attribute, sorted. Bars in red indicate **super-additive** (synergistic) attributes, bars in blue **sub-additive** (masking) attributes.

**Read-off:** Across all 15 binary mixtures × 15 aroma attributes = **225 paired predictions**, observed and predicted frequencies are strongly correlated (**r = 0.812, p = 6.0 × 10⁻⁵⁴; RMSE = 0.100; mean bias +0.016**), indicating that CATA response frequency is, to a first approximation, additive. The four most consistently super-additive attributes are **sweet (+0.068), grassy (+0.068), musky (+0.065) and fruity (+0.062)**; the four most consistently sub-additive are **floral (−0.037), chemical/solvent (−0.032), roasted (−0.020) and cooling (−0.015)**.

**Essential caveat for the text (do not omit):** *These data describe **attribute check-frequency**, not perceived intensity. Mixture intensity ratings were not collected (the 强度 column of 混合香气平均值 is empty for all 65 samples and 混合香气原始数据 contains no intensity column). Consequently Figure 2 supports only a **qualitative** statement about additivity of attribute salience; a quantitative isobologram or σ–τ analysis cannot be computed from this data set.*

---

## Figure 3A / 3B — Ridgeline of intensity distributions

> 文件：`fig_ridgeline_A.png`、`fig_ridgeline_B.png`（同版式，各 7 个化合物）

**Figure 3. Distribution of perceived intensity across the dilution series. (A, B)** Each row is one compound; within a row, each filled ridge is the kernel density estimate (Gaussian kernel, bandwidth 0.13) of the 18 panellists' intensity ratings at one concentration level, stacked from C1 (lowest concentration, bottom) to the highest level (top). Fill colour follows the `viridis` gradient across levels. The red dashed line marks the mean rating at that level; the right-hand annotation gives the level code, the number of panellists and the mean ($\bar{x}$). Only compounds with at least five concentration levels and complete 18-panellist records are shown (14 compounds, 87 concentration levels, 1,566 observations). The two panels are split by level count to keep the aspect ratio printable.

**Read-off:** Distributions are visibly **right-skewed at low concentrations and left-skewed near the ceiling**, i.e. the response distribution changes shape systematically with dose rather than merely shifting. This is the visual justification for treating intensity as a bounded, heteroscedastic response and for reporting medians/intervals alongside means.

**Note on the original figure:** the previous version of this figure pooled two independent scoring rounds for phenylethanol C1–C6 (36 rows per level instead of 18). The duplicate round has been removed here (see Methods 2.2).

---

## Figure 4 — Panel heterogeneity and the bimodality test

> 文件：`fig_cluster_bimodality.png`

**Figure 4. Assessment of panellist heterogeneity and of the "sensitive vs. insensitive" hypothesis. (A)** Heat map of mean intensity for 18 panellists (rows, ordered by each panellist's grand mean) × 14 compounds (columns); red horizontal lines mark cluster boundaries. **(B)** Principal-component projection of the panellists with *k*-means partitioning; ellipses/blobs are labelled by cluster. **(C)** Silhouette coefficient as a function of the number of clusters *k*, for two feature representations (per-compound mean response; per-compound detection threshold), with the conventional *s* = 0.5 reference line. **(D)** Bayesian information criterion (BIC) for Gaussian-mixture models with *k* = 1–6 components fitted to the 18 panellists' grand mean intensity; lower is better. **(E)** Per-concentration-point bimodality: ΔBIC = BIC(1-component) − BIC(2-component) for each of the 87 concentration levels, plotted against that level's mean intensity; grey points are levels at which ≥50 % of ratings hit the scale ceiling. **(F)** The single level that survives all five strict criteria: raw histogram of the 18 ratings with the fitted two-component density.

**Read-off — this figure is a *negative* result and should be presented as such:**

| Test | Result | Interpretation |
|---|---|---|
| Silhouette (all *k*, both features) | 0.109 – 0.153 | ≪ 0.5 → **no reliable cluster structure** |
| GMM on grand means | BIC(*k*=1) = 23.93 < BIC(*k*=2) = 31.47; ΔBIC = −7.54 | **favours a single population**; grand-mean range is only 5.41–6.95 (SD 0.41) |
| Per-level ΔBIC > 10 | 38 / 87 (44 %) | inflated by ceiling truncation — 7/7 ceiling levels pass, 0/80 further pass all strict criteria except one |
| Strict 5-criterion test | **1 / 87 (1 %)** | only isoamyl isobutyrate C6 (mean 6.06, SD 2.41, ΔBIC +81.32, component means 4.6 and 9.0, empty bins at 7–8) |

**Conclusion to be written into the Discussion:** *the data support continuous inter-individual variation in odour-intensity responses, but do **not** support the existence of two distinct classes of "sensitive" and "insensitive" panellists. Apparent bimodality at most concentration levels is an artefact of the bounded 0–12 scale (ceiling truncation at high dose, zero-inflation at low dose).*

**Footnote:** *The between-cluster t-test shown in the analysis log is not reported, because clusters were constructed from the same intensity variables and the comparison is therefore circular.*

---

## Supplementary figures

**Figure S1. Anatomy of a single dose–response panel** (`demo_single_panel.png`). Layer-by-layer deconstruction used as the template for Figure 1: observed points with 95 % CI, fitted four-parameter log-logistic curve, asymptote guides, EC₅₀ drop-lines, and annotation zones.

**Figure S2. Small multiples of all fitted curves** (`demo_small_multiples.png` / `.pdf`). One panel per compound, shared axes, arranged by potency; identical styling to Figure S1 so that the panel and the summary map can be cross-read.

**Figure S3. Forest plot of EC₅₀ with 95 % confidence intervals** (`demo_ec50_forest.png`). Compounds ordered by potency; intervals spanning more than one order of magnitude are highlighted. This is the honest counterpart to Figure 1's bubble positions, because it exposes which EC₅₀ estimates are actually identified by the data.

---

# 2. Methods（英文，可直接粘贴）

## 2.1 Sensory panel and sample set

A trained sensory panel of **18 assessors** (coded S01–S18) evaluated the samples in individual, well-ventilated booths under white light. Odour intensity was rated on a **bounded 12-point category scale (0 = not perceived, 12 = very strong)**, and aroma character was recorded simultaneously using a **check-all-that-apply (CATA) vocabulary of 15 attributes** plus an open "other" category.

Two sample sets were used. **(i) Single-compound dilution series.** Fourteen aroma compounds were each presented at **four to eight concentration levels** prepared by ten-fold serial dilution, giving **87 concentration levels in total** (0.1998 ppb to 2.5 × 10⁸ ppb; Supplementary Table S2). Each level was rated for intensity by all 18 panellists, yielding 1,566 analyzable observations. **(ii) Mixtures.** A six-component, six-factor full factorial design was used: six key compounds (dimethyl trisulfide, isoamyl acetate, ethyl butyrate, acetaldehyde diethyl acetal, hexanoic acid and 1-hexanol) were combined at low and high levels, producing a common all-low baseline (母液1), 6 single-high samples and all 15 binary-high combinations (2⁶ = 64 design points in total).

## 2.2 Data curation and quality control

All analysis was performed on the raw per-panellist records rather than on the tabulated means. The following pre-specified curation steps were applied:

1. **Identification of samples and assessors** by forward-filling the sample label, then matching assessors against the regular expression `^S\d+$`; compound names were normalised by stripping whitespace, converting full-width brackets and applying an alias table (二甲基三硫醚 → 二甲基三硫; 反，顺-2，6-壬二稀醇 → 反，顺-2,6-壬二烯醛).
2. **Duplicate-round detection.** Of the 1,674 valid raw records, **216 records (6 keys × 18 assessors)** were duplicated: phenylethanol levels C1–C6 each carried two complete scoring rounds. The second round was **excluded**, leaving 1,566 observations, because test–retest agreement between the two rounds was negligible (**Pearson r = −0.077 on average; range −0.349 to +0.284**), far below the *r* ≥ 0.7 convention for repeatable intensity ratings. The excluded round is retained in the archived file `clean_assessor_level_两轮合并.csv` for transparency.
3. **Range checking.** One out-of-range value was found (phenylethanol C6, assessor S12, rating = 15 on a 0–12 scale) and was confined to the excluded round; no out-of-range values remain in the analysis set.
4. **Reconciliation against the reported means.** Recomputed panellist-level means were compared cell-by-cell with the reported summary sheet. Of the 87 concentration keys, **83 could be matched and 60 reproduced exactly (|Δ| < 10⁻⁶)**; **23 differed**, the largest discrepancy being 4.50 scale units (ethyl isobutyrate C6). Nine keys present in the reported summary have **no supporting raw record** (n-butanol C1–C7, 2-methylbutyrate C6, acetal C2), and four raw keys were **absent from the summary** (butyric acid C1 and C6, ethyl isobutyrate C8, phenylethanol C7). These discrepancies are itemised in `QC_数据质量报告.txt`.
5. **Missingness.** Mixture intensity could not be analysed: the intensity column is empty for all 65 mixture samples in the summary sheet, and the corresponding raw sheet contains no intensity column at all.

## 2.3 Dose–response modelling

Concentration–response data for each compound were fitted by non-linear least squares to a **four-parameter log-logistic (4PL) model with a free lower asymptote**:

$$
I(C) \;=\; d \;+\; (I_{\max} - d)\,\frac{C^{\,n}}{EC_{50}^{\,n} + C^{\,n}}
$$

where *d* is the lower asymptote, *I*ₘₐₓ the upper asymptote, *n* the Hill coefficient and EC₅₀ the concentration giving half of the fitted dynamic range. Concentrations were taken from the formulation sheet (`Sheet5`, ppm/ppb, ten-fold serial dilutions) and responses were the per-level panellist means. Initial values for *n* (0.3, 1, 3) and for EC₅₀ (geometric mean and the range endpoints) were scanned, bounds were *d*, *I*ₘₐₓ ∈ [−30, 60], *n* ∈ [0.05, 20] and EC₅₀ ∈ [*C*ₘᵢₙ/10³, *C*ₘₐₓ × 10³], and the solution with the highest *R*² was retained.

**Model verification.** The choice of 4PL rather than a three-parameter logistic (upper asymptote fixed at 0) is not arbitrary: an independent refit with the three-parameter form failed to reproduce the archived parameter table (0/15 compounds agreed within 10 % for EC₅₀; median ratio 810), whereas the **4PL form reproduced all 15 *R*² values to four decimal places**, 14/15 *I*ₘₐₓ and 14/15 *n* within 10 %, and 10/15 EC₅₀ within 10 %. The five EC₅₀ values that did not reproduce (dimethyl trisulfide, *trans,cis*-2,6-nonadienal, ethyl phenylacetate, ethyl isobutyrate and phenylethanol) are exactly those whose fitted curves are flat (*n* < 0.6) and whose profile-likelihood confidence intervals are degenerate or span more than an order of magnitude — i.e. they are **structurally non-identifiable**, not merely imprecise.

**Goodness of fit and caveats.** *R*² ranged from 0.950 to 1.000. Butyric acid returned *R*² = 1.000 from only four concentration levels fitted with four parameters, which is a saturated fit and carries no residual degrees of freedom. Five compounds returned an *I*ₘₐₓ above the maximum scoreable intensity of 12 (phenylethanol 19.83, ethyl phenylacetate 13.87, *trans,cis*-2,6-nonadienal 12.72, dimethyl trisulfide 12.71, ethyl isobutyrate 12.70); these are extrapolations beyond the measurement range and are reported as such. Three Hill coefficients returned confidence intervals of zero width (dimethyl trisulfide 0.292, *trans,cis*-2,6-nonadienal 0.489, n-butanol 1.255), and four EC₅₀ intervals span between 10- and 243-fold. All such estimates are flagged in Figure 1 and Figure S3.

## 2.4 Mixture additivity

Because mixture intensity was not collected, additivity was assessed on CATA attribute frequencies. For each of the 15 binary-high design points and each of the 15 aroma attributes, an additive prediction was formed as

$$
\hat p \;=\; \operatorname{clip}_{[0,1]}\bigl[p_0 + (p_A - p_0) + (p_B - p_0)\bigr],
$$

where *p*₀ is the baseline (all-low) frequency and *p*_A, *p*_B the frequencies of the matching single-high samples. Observed and predicted values were compared by Pearson correlation, RMSE and mean signed deviation overall, and by mean signed deviation within each attribute. This is a **qualitative** test of attribute-salience additivity; a formal isobologram or σ–τ analysis requires mixture intensity ratings and paired concentration levels for both components of every pair, neither of which is available here.

## 2.5 Distributional (ridgeline) visualisation

For every compound and concentration level, a Gaussian kernel density estimate was computed over the 18 individual intensity ratings (bandwidth 0.13, grid of 400 points spanning −0.3 to 12.3) and normalised to unit peak height. Densities were stacked with a fixed vertical step and overlaid with the level mean, producing the ridgeline layout of Figure 3. Density estimation and plotting were performed per level without pooling across levels, so that the change in distributional shape with dose is preserved.

## 2.6 Panel clustering and bimodality testing

Two feature representations of the 18 panellists were examined: **(a)** the per-compound mean intensity (18 × 14 matrix) and **(b)** the per-compound detection threshold, defined as the lowest concentration level at which a panellist assigned a rating ≥ 6 (never reaching the threshold coded as one level above the highest tested). Both were z-standardised. The number of clusters *k* = 2…6 was selected by the **silhouette coefficient** under *k*-means (100 restarts, fixed seed), and the two-dimensional structure was inspected by principal-component analysis.

Bimodality was tested at two levels. **(i) Panellist level:** one- to six-component Gaussian mixture models (50 initialisations) were fitted to the 18 grand mean intensities and compared by BIC; ΔBIC = BIC(*k*=1) − BIC(*k*=2), so that positive values favour two components, with |ΔBIC| > 10 conventionally regarded as strong evidence. **(ii) Concentration level:** the same comparison was applied to the 18 ratings at each of the 87 levels. Because a bounded scale generates spurious bimodality, each of the 38 levels with ΔBIC > 10 was additionally required to satisfy five criteria: (1) ΔBIC > 10; (2) smaller-mixture weight 0.20 ≤ *w* ≤ 0.80; (3) separation of component means |μ₁ − μ₂| ≥ 3 scale units; (4) level mean within 3–9 (away from floor and ceiling); (5) at least one empty integer bin of width ≥ 3 between the two component means, i.e. a genuine gap rather than a skew.

Between-cluster differences were **not** tested inferentially, because clusters are constructed from the same variables and such a test would be circular.

## 2.7 Software and reproducibility

All analyses were performed in **Python 3.12** using pandas, NumPy, SciPy (`scipy.optimize.curve_fit`, `scipy.stats`), scikit-learn (`KMeans`, `GaussianMixture`, `PCA`, `silhouette_score`) and Matplotlib 3.9 (PNG at 300 dpi, PDF with embedded TrueType, `pdf.fonttype = 42`). Random seeds were fixed at 0 for all stochastic procedures. Scripts, intermediate tables and quality-control reports are archived alongside the figures; the raw workbooks were opened read-only and were not modified.

## 2.8 Declarations

**Data availability.** All raw records, derived tables and quality-control reports are provided with this submission (`教学示范图/`), including `clean_assessor_level.csv` (analysis-ready, 1,566 rows), `qc_mean_vs_table.csv`, `qc_round1_vs_round2.csv`, `bimodality_GMM_by_level.csv`, `bimodality_strict.csv`, `assessor_clusters.csv`, `EC50_summary.csv`, `refit_4PL.csv`.

**Ethics.** The study involved only the sensory evaluation of odour standards by adult trained panellists; no human-subject risk beyond ordinary olfactory exposure at or below documented odour-threshold levels. *(Adapt to the institution's actual ethics statement before submission.)*

**Conflicts of interest / Funding.** *To be completed by the authors.*

---

# 3. 数字溯源表

| # | 数字 | 出处 | 状态 |
|---|---|---|---|
| 3.1 | 4PL 模型复现：R² 15/15 完全一致，Imax/n 14/15 <10%，EC50 10/15 <10% | `refit2.py` → `refit2.txt`、`教学示范图\refit_4PL.csv` | 🟢 已独立验证 |
| 3.2 | EC50 范围 8.88 ppb ~ 7.98×10⁷ ppb | `EC50_summary.csv` | 🟢 |
| 3.3 | Imax > 12 的 5 个化合物及数值 | `EC50_summary.csv` + `gen4_log.txt` | 🟢 |
| 3.4 | 丁酸 4 点 R²=1.000（饱和拟合） | `EC50_summary.csv`（n_points=4, R2=1.0） | 🟢 |
| 3.5 | 3 个 n 的 CI 宽度为 0 | `EC50_summary.csv` n_CI 列 | 🟢 |
| 3.6 | EC50 CI 跨度 242×（2-甲基丁酸乙酯 1648→398920） | `EC50_summary.csv` EC50_CI 列 | 🟢 |
| 3.7 | 异丁酸乙酯 CI 下界 1,008,758 > 点估计 770,476（CI 不含点估计） | `EC50_summary.csv` | 🟠 矛盾，须查拟合的 CI 计算方法 |
| 3.8 | 加和性 225 点、r=0.812、p=6.02e-54、RMSE=0.0998、偏差 +0.0157 | `gen4.py` L154–193 → `gen4_log.txt` | 🟢 |
| 3.9 | 各属性 mean_dev（甜香 +0.0683 … 花香 −0.0371），n=15/属性 | `gen4_log.txt` L6–20 | 🟢 |
| 3.10 | 混合物无强度数据 | `混合香气平均值` 强度列 65 行全空；`混合香气原始数据` 1261 行无强度列 | 🟢 已验证 |
| 3.11 | 1674 → 1566 条（去掉苯乙醇 C1–C6 第二轮 108 条） | `QC_数据质量报告.txt` L41–42 | 🟢 |
| 3.12 | 两轮重测 r 均值 −0.077，范围 −0.349 ~ +0.284 | `QC_数据质量报告.txt` L20–27 | 🟢 |
| 3.13 | 越界值 1 条（苯乙醇 C6 / S12 = 15） | `QC_数据质量报告.txt` L6–7 | 🟢（在被剔除的第二轮中） |
| 3.14 | 与论文均值表：83 可比 / 60 完全一致 / 23 有差异，最大 Δ=4.50 | `QC_数据质量报告.txt` L47–52 | 🟢 |
| 3.15 | 表有原始无 9 个 key；原始有表无 4 个 key | `QC_数据质量报告.txt` L76–80 | 🟢 |
| 3.16 | 87 浓度点、14 化合物、5–8 级/化合物、1566 条 | `QC_双峰与命名.txt` / `clean_assessor_level.csv` | 🟢 |
| 3.17 | Sheet5 = 91 行浓度阶梯，0.1998 ~ 2.5×10⁸ ppb | `refit.txt` L1 | 🟢 |
| 3.18 | silhouette 0.109–0.153（全 <0.5） | `QC_聚类与双峰检验.txt` L15–16 | 🟢 |
| 3.19 | GMM BIC：k=1 最优 23.93，k=2 为 31.47，ΔBIC=−7.54 | `QC_聚类与双峰检验.txt` L21–31 | 🟢 |
| 3.20 | 38/87 ΔBIC>10；严格 5 判据后 **1/87** | `QC_聚类与双峰检验.txt` L37、L78 | 🟢 |
| 3.21 | 异丁酸乙酯 C6 原始 18 值、μ差 4.4、空档跨度 3 | `QC_聚类与双峰检验.txt` L92–93 | 🟢 |
| 3.22 | 全局均值范围 5.41–6.95（sd 0.41） | `QC_聚类与双峰检验.txt` L32 | 🟢 |
| 3.23 | k=3 簇规模 5/8/5，均值 6.60/5.89/6.06 | `QC_聚类与双峰检验.txt` L100–102 | 🟡 结构极弱（silhouette 0.153），图注中已按"无可靠结构"表述 |

---

# 4. 提交前必须核对的事项

| 级别 | 事项 |
|---|---|
| 🔴 P0 | **图 1 的 EC₅₀ 置信区间计算方法需说明**：异丁酸乙酯的 CI 下界 (1,008,758) 大于点估计 (770,476)，CI 不包含点估计，说明 CI 是用不同方法/不同参数化算的，或有实现错误 |
| 🔴 P0 | **图 3 原版的苯乙醇重复计数**（36 条/级）已在这版修正；若论文仍用旧图，需重做 |
| 🟠 P1 | 23 个浓度点的均值与论文均值表对不上（最大差 4.50 分），需说明是取整、剔除规则还是笔误 |
| 🟠 P1 | 9 个均值表中的 key 没有原始记录支撑（正丁醇 7 个 + 2-甲基丁酸乙酯 C6 + 乙缩醛 C2）→ 这些点无法复现 |
| 🟠 P1 | 两轮评分的重测 r ≈ −0.077，需在方法学里说明第二轮的处置与理由（目前我们已剔除并留档） |
| 🟠 P1 | 图 2 的结论只能写"attribute salience 大体加和"，**不能**写成强度加和或引用 σ–τ |
| 🟠 P1 | 论文正文 p19「二甲基三硫极低响应值」与附表 1 的 0.84–0.90 仍矛盾（前次结论，未处理） |
| 🟡 P2 | 附表 1 表题 OVA → OAV；正文"酯 13 / 醇 4"→ 14 / 3 |
| 🟡 P2 | 三套阈值体系（附表1 Kd、OAV 反推、重构分析 Threshold）需统一 |
| 🟡 P2 | 命名不一致两处：二甲基三硫醚↔二甲基三硫；反，顺-2，6-壬二稀醇↔反，顺-2,6-壬二烯醛（且"醇"应为"醛"） |
| 🟡 P2 | Methods 2.8 的伦理、基金、利益冲突声明需作者补写 |
