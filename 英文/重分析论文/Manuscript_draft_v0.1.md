# Preliminary Manuscript Draft v0.1

> **Draft status.** This is an English-language working draft prepared from the files currently in PAPER. Square-bracketed items are author checks, not established facts. The dose–response parameters and mixture additivity statistics below are provisional until the raw data are reconciled and the final analysis is rerun. Do not submit this version.

## Title

**Concentration-dependent sensory responses to aroma compounds and CATA-based profiling of mixtures: a data-auditable reanalysis**

## Abstract

**Background:** Sensory responses to aroma compounds depend on concentration, whereas mixture descriptions may be recorded using different response formats. These outcomes should be analysed and interpreted separately.

**Objective:** This study re-examined archived sensory records to characterize concentration–intensity relationships for single aroma compounds and to assess whether check-all-that-apply (CATA) descriptor frequencies in mixtures were consistent with a simple additive expectation.

**Materials and methods:** Records from 18 assessors were reviewed. The cleaned single-compound data set contains 1,566 intensity observations from 14 compounds across 87 concentration levels. A separate mixture data set contains CATA responses for mixtures prepared from six aroma compounds. Single-compound responses were evaluated using a four-parameter log-logistic model as a candidate dose–response model. For 15 binary-mixture conditions and 15 descriptors, observed CATA frequencies were compared with predictions based on the baseline and corresponding single-high samples.

**Results:** An independent refit of the archived summary table reproduced the archived coefficients of determination under a four-parameter model, but several concentration-response parameters were weakly identified; five archived upper-asymptote estimates exceeded the 0–12 rating scale. The exploratory comparison of 225 mixture–descriptor pairs yielded a Pearson correlation of 0.812 and an RMSE of 0.100 between observed and predicted CATA frequencies. These are descriptive summaries; assessor-level uncertainty remains to be recalculated. Clustering diagnostics did not support a robust division of assessors into “sensitive” and “insensitive” groups (silhouette coefficients, 0.109–0.153; one-component GMM favored over two components for the assessor-level grand means).

**Conclusion:** The single-compound records are suitable for a transparent reanalysis of dose-dependent intensity, provided that weakly identified parameters and scale-boundary effects are reported. The mixture records support conclusions about descriptor selection frequency only; they do not establish intensity additivity, sensory synergy, or masking. The integrated interpretation remains provisional pending data reconciliation and confirmation of the original sensory protocol.

**Keywords:** aroma compounds; sensory evaluation; dose–response; CATA; mixture; assessor variability

## 1. Materials and methods

### 1.1. Data sources and study records

This work is a reanalysis of archived records associated with two undergraduate research projects on single and mixed aroma compounds. The source workbooks and analysis outputs are archived in the PAPER repository. Original workbooks were treated as read-only.

The single-compound intensity records include 14 compounds, 87 compound–concentration levels, and 18 assessor codes (S01–S18). The recorded response is an odour-intensity score on a bounded 0–12 scale. **[AUTHOR TO VERIFY: scale anchors, assessor training and recruitment, age range, and whether the same assessors evaluated all samples.]**

The mixture records concern six compounds: dimethyl trisulfide, isoamyl acetate, ethyl butyrate, acetaldehyde diethyl acetal, hexanoic acid, and 1-hexanol. The archived design sheet describes two concentration levels for each factor. The mixture summary contains 65 sample records, while the design notes also refer to a 2⁶ design. **[AUTHOR TO VERIFY: the exact number of unique samples, baseline/control rows, any repeated samples, the meaning of the two levels, and the mapping between design codes and raw CATA records.]**

The archived analysis compares 15 binary mixtures across 15 descriptors. **[AUTHOR TO VERIFY: the final descriptor list, descriptor translations, and whether the descriptor set is identical across the single and mixture tasks.]**

### 1.2. Sensory procedure and sample preparation

The archived intensity data use a 0–12 response scale, and the mixture data contain CATA selections. Details required for reproducibility were not fully established from the analysis files and must be recovered from the original thesis, ballot, and laboratory records:

- solvent or model matrix and ethanol concentration;
- compound purity, preparation, final concentration units, and dilution procedure;
- sample volume, vessel type, headspace/equilibration time, and serving temperature;
- sample randomization, blinding, presentation order, rest intervals, and palate/olfactory recovery procedure;
- assessor instructions, training, and intensity-scale anchors;
- ethics approval or exemption and participant consent.

No ethics, safety, or panel-training statement should be inserted until verified against the institutional record. **[AUTHOR TO COMPLETE.]**

### 1.3. Data curation and quality control

Raw observations were linked to sample and assessor identifiers, compound names were harmonized using a documented alias table, and concentration labels were checked against the formulation sheet. The analysis-ready single-compound file contains 1,566 assessor-level observations, corresponding to 18 assessors at each of 87 concentration levels.

The source workbook contains 1,674 valid rows. An additional 108 rows occur in phenylethanol C1–C6 and have provisionally been interpreted as a second scoring round. The average correlation between the two putative rounds was −0.077 (range, −0.349 to 0.284). The current analysis retains one record per assessor and level and excludes the additional round. **[AUTHOR TO VERIFY the row order, round identifiers, and exclusion rationale from original records before treating this as a prespecified data-cleaning decision.]**

One score of 15 was identified on a nominal 0–12 scale in the excluded records. No out-of-range score remains in the current analysis-ready set. Recomputed means were also compared with the archived summary table: 83 keys could be matched, 60 matched exactly, and 23 differed; the largest absolute difference was 4.50 scale units. Nine summary keys lacked corresponding raw records, and four raw keys were absent from the summary table. These discrepancies must be resolved or explicitly accounted for before the manuscript data set is frozen.

### 1.4. Single-compound dose–response analysis

For each compound, the candidate four-parameter log-logistic model was:

I(C) = d + (u − d) × Cⁿ / (EC50ⁿ + Cⁿ),

where C is concentration, d is the lower asymptote, u is the upper asymptote, n is the slope parameter, and EC50 is the concentration halfway between the fitted asymptotes.

The existing archived parameter table was independently checked against three- and four-parameter formulations. This check reproduced the archived R² values under the four-parameter formulation, but it did not constitute a final model-selection analysis using the reconciled assessor-level data. **[FINAL ANALYSIS REQUIRED: refit from the frozen raw-data set; examine residuals and parameter identifiability; compare candidate models using an appropriate criterion; estimate uncertainty by resampling assessors or a model that accounts for repeated observations within assessors.]**

The final dose–response figure should show the observed assessor-level data or concentration-level summaries, the fitted curve, and uncertainty bands. EC50 estimates should be reported only where the data identify them adequately. An upper asymptote beyond 12 is an extrapolation beyond the response scale, not an attainable sensory score.

### 1.5. Mixture CATA frequencies and additive predictions

For each mixture and descriptor, CATA frequency was calculated as the proportion of valid assessor responses selecting that descriptor. For a binary mixture, the archived additive prediction was:

p̂ = clip[0,1]{p₀ + (pA − p₀) + (pB − p₀)},

where p₀ is the all-low baseline frequency and pA and pB are the corresponding single-high frequencies.

Observed and predicted frequencies were summarized using Pearson correlation, root-mean-square error (RMSE), and mean signed deviation. Because the 225 mixture–descriptor pairs are not 225 independent assessors, inferential uncertainty should be estimated using the assessor as the resampling unit or an appropriate repeated-measures model. **[FINAL ANALYSIS REQUIRED: compute assessor-level intervals and confirm the correct mapping of baseline and single-high samples.]**

A correspondence-analysis map and a sample-by-descriptor frequency heat map are planned for the final manuscript. These analyses have not yet been run on the reconciled mixture table and should not be described as completed results.

### 1.6. Assessor heterogeneity

The archived exploratory analysis examined assessor-level mean responses and concentration-level response distributions. K-means silhouette coefficients and Gaussian mixture models (GMMs) were used to assess whether the observed data supported discrete assessor groups. Concentration-level bimodality was interpreted cautiously because the intensity scale is bounded and can produce floor or ceiling concentrations of ratings. These exploratory diagnostics were not treated as confirmatory evidence for latent assessor classes.

### 1.7. Software and reproducibility

The archived scripts use Python 3.12 with pandas, NumPy, SciPy, scikit-learn, and Matplotlib. The final manuscript should report the versions actually used for the final analysis and should archive the exact scripts, input-data checksums, and derived tables. **[AUTHOR TO VERIFY execution environment and package versions.]**

## 2. Results

### 2.1. Data quality and analysis set

The current cleaned file contains 1,566 single-compound intensity observations from 18 assessors, 14 compounds, and 87 concentration levels. The source records include a second set of phenylethanol ratings at six concentration levels; these ratings were excluded provisionally because the two rounds showed poor agreement. Comparison with the archived mean table identified 23 mismatches among 83 comparable keys. The largest discrepancy was for ethyl isobutyrate C6 (recalculated mean 6.056 versus archived mean 10.556). These differences prevent treating every archived summary value as a verified observation.

### 2.2. Single-compound concentration–intensity responses

The archived dose–response parameter table reports EC50 values spanning approximately 8.9 to 8.0 × 10⁷ ppb and R² values from 0.950 to 1.000. An independent refit of that summary table reproduced all archived R² values using a four-parameter log-logistic formulation, whereas a three-parameter formulation did not reproduce the archived EC50 values. This is a consistency check against the archived table, not independent confirmation that the archived parameter estimates are all reliable.

Five archived upper-asymptote estimates exceed the scale maximum of 12: phenylethanol, ethyl phenylacetate, trans,cis-2,6-nonadienal, dimethyl trisulfide, and ethyl isobutyrate. The archived results also identify several shallow curves for which EC50 is weakly identified. In addition, the stored 95% confidence interval for ethyl isobutyrate has a lower bound above its point estimate. **These parameter estimates and intervals must be recalculated before the final manuscript reports them.**

### 2.3. CATA additivity in binary mixtures

Across the archived 15 binary mixtures and 15 descriptors, the 225 observed–predicted pairs had a descriptive correlation of 0.812 and an RMSE of 0.0998; the mean signed deviation was +0.0157 and the mean absolute deviation was 0.0669. The archived attribute summaries showed positive mean deviations for sweet, grassy, musky, and fruity descriptors and negative deviations for floral, chemical/solvent, roasted, and cooling descriptors.

These values summarize CATA selection frequencies. Assessor-level uncertainty has not yet been recomputed, so the attribute-level directions should be treated as exploratory rather than as statistically established synergy or masking effects. Mixture intensity ratings are absent from the archived mixture sheets.

### 2.4. Assessor clustering and bimodality

The archived clustering analysis yielded silhouette coefficients from 0.109 to 0.153 across the tested feature sets and values of k. For the distribution of assessor grand means, the one-component GMM had a lower BIC than the two-component model (23.93 versus 31.47; BIC₁ − BIC₂ = −7.54). Of 87 concentration-level distributions, one met all five archived strict criteria for bimodality. These findings do not support a stable division of the panel into “sensitive” and “insensitive” groups.

## 3. Discussion

### 3.1. Interpret dose–response parameters in light of the observed data

A fitted dose–response curve can summarize how perceived intensity changes with concentration, but a high R² alone does not establish that EC50, slope, or the upper asymptote is well determined. This is particularly important for shallow curves and for fits with few concentration levels. Sensory dose–response work emphasizes the need to consider the response format, concentration range, and experimental design when interpreting the fitted function [1].

The upper asymptote is constrained by the 0–12 rating scale. Estimates above 12 therefore describe model extrapolation and should not be presented as observed sensory intensity. The final figure should show raw observations and confidence bands so that readers can see where the model is supported by data. The anomalous ethyl-isobutyrate EC50 interval and mismatches with the archived means require resolution before compound potency rankings are interpreted.

### 3.2. CATA frequencies describe descriptor selection, not intensity additivity

The mixture data reveal how often assessors selected each descriptor, but they do not contain mixture intensity scores. CATA frequency can be informative for sensory profiling, yet its relationship to perceived intensity is not a direct one-to-one conversion and may vary across assessors and attributes [2]. Thus, the observed–predicted comparison can be discussed as an exploratory test of additivity in descriptor-selection frequency. It cannot establish additive mixture intensity, a classical isobologram, or sensory synergy.

A correspondence-analysis map and frequency heat map would make the mixture profiles easier to interpret. For CATA data, correspondence analysis has been evaluated as a method for constructing sensory spaces [3]. Any inferential statements should be based on assessor-level uncertainty rather than treating each mixture–descriptor pair as independent.

### 3.3. The current data do not establish discrete sensitivity groups

The low silhouette values, preference for one GMM component in assessor grand means, and the single strict bimodal concentration point do not support a robust sensitive/insensitive split. The scale boundaries are also relevant: floor and ceiling responses can distort apparent distribution shapes. The current evidence is more consistent with continuous assessor variation, although stronger conclusions would require a design with additional mid-range concentrations and a larger panel.

### 3.4. Strengths and limitations

A strength of this reanalysis is the use of raw assessor-level records and explicit comparison with archived summaries. The work also separates intensity responses for single compounds from CATA selections for mixtures.

The main limitations are material. First, the raw records and archived means disagree for multiple keys, and some summary values lack raw support. Second, the second phenylethanol scoring round and the exact mixture-design row count require confirmation from source records. Third, the archived mixture data contain no intensity outcome, and the mixture compounds have limited overlap with the single-compound concentration series. These constraints weaken a unified mechanistic interpretation of single-compound potency and mixture interactions. Finally, panel characteristics, presentation protocol, ethics documentation, and the final sample-preparation details remain to be completed.

**Working conclusion:** The cleanest defensible story is a dose–response reanalysis of single-compound ratings accompanied by a separate, clearly bounded CATA analysis of mixtures. Whether these two datasets should remain in one article or become separate manuscripts should be decided after the data reconciliation and novelty review.

## References

1. Kershaw JC, Running CA. Dose–response functions and methodological insights for sensory tests with astringent stimuli. *Journal of Sensory Studies*. 2019;34(1):e12480. https://doi.org/10.1111/joss.12480
2. Jaeger SR, Chheang SL, Jin D, Ryan GS, Ares G. How do CATA questions work? Relationship between likelihood of selecting a term and perceived attribute intensity. *Journal of Sensory Studies*. 2023;38(4):e12833. https://doi.org/10.1111/joss.12833
3. Vidal L, Tárrega A, Antúnez L, Jaeger SR, Ares G. Comparison of correspondence analysis based on Hellinger and chi-square distances to obtain sensory spaces from check-all-that-apply (CATA) questions. *Food Quality and Preference*. 2015;43:106–112. https://doi.org/10.1016/j.foodqual.2015.03.003
4. Llobell F, Bonnet L, Giacalone D. Assessment of panel performance in CATA and RATA experiment. *Journal of Sensory Studies*. 2024;39(4):e12941. https://doi.org/10.1111/joss.12941

## Author checks required before v0.2

- Confirm the study protocol, panel training, scale anchors, sample matrix, concentration units, serving conditions, ethics statement, and consent.
- Reconcile the 64-design-point description with the 65 mixture records and map every sample code to its formulation.
- Resolve the 23 mismatches between raw-data means and archived means; document the nine unsupported and four unlisted keys.
- Confirm that the 108 additional phenylethanol records are a distinct second round and that excluding them is scientifically justified.
- Refit dose–response models from the frozen assessor-level data; correct the ethyl-isobutyrate interval and report uncertainty only for identifiable parameters.
- Recompute mixture additivity uncertainty with assessors as the resampling units; remove unsupported p-values and significance language.
- Run the planned CATA correspondence analysis and frequency heat map, then decide whether the two datasets form one coherent manuscript.
- Complete the literature review and choose a target journal before adapting formatting and claims.
