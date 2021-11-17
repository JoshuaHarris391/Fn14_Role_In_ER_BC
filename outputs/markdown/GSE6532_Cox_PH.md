# GSE6532 Cox Proportional Hazards Analysis
- Goal is to evaluate whether Fn14 is significantly associated with recurrence in patients treated with tamoxifen monotherapy.
- It would also be interesting to see Fn14's association with recurrence in patients without tamoxifen treatment.
- Tertiles were coded as semi-continuous variables

# Patients with Tamoxifen monotherapy

## Univariate Cox Models

- [ER IHC Status - RFS](../GSE6532/Cox/Tamoxifen/Uni_Cox_ER_IHC.csv)
- [ESR1 mRNA tertile expression - RFS](../GSE6532/Cox/Tamoxifen/Uni_Cox_ESR1_TERT.csv)
- [TNFRSF12A mRNA tertile expression - RFS](../GSE6532/Cox/Tamoxifen/Uni_Cox_TNFRSF12A_TERT.csv)

Note: ER and ESR1 are unlikely to discern patients with recurrence as almost all patients that reccieved tamoxifen monotherapy are ER IHC positive.

## Multivariate Cox Models
Covariates consist of age (continuous), node status (binary, 0,1), grade (semi-continuous, 1,2,3), and tumour size (continuous). PAM50 subtype is coded as a categorical variable with Normal-like as the baseline.

- [ER IHC Status - RFS](../GSE6532/Cox/Tamoxifen/Multi_Cox_ER_IHC.csv)
- [ESR1 mRNA tertile expression - RFS](../GSE6532/Cox/Tamoxifen/Multi_Cox_ESR1_TERT.csv)
- [ESR1 mRNA tertile expression + PAM50 - RFS](../GSE6532/Cox/Tamoxifen/Multi_Cox_ESR1_TERT_PAM50.csv)
- [TNFRSF12A mRNA tertile expression - RFS](../GSE6532/Cox/Tamoxifen/Multi_Cox_TNFRSF12A_TERT.csv)
- [TNFRSF12A mRNA tertile expression + PAM50 - RFS](../GSE6532/Cox/Tamoxifen/Multi_Cox_TNFRSF12A_TERT_PAM50.csv)

-------------------------

# Patients without Tamoxifen monotherapy

## Univariate Cox Models

- [ER IHC Status - RFS](../GSE6532/Cox/Untreated/Uni_Cox_ER_IHC.csv)
- [ESR1 mRNA tertile expression - RFS](../GSE6532/Cox/Untreated/Uni_Cox_ESR1_TERT.csv)
- [TNFRSF12A mRNA tertile expression - RFS](../GSE6532/Cox/Untreated/Uni_Cox_TNFRSF12A_TERT.csv)

Note: The ESR1 high expression tertile is not representative of ER IHC positivity as ESR1 tertiles were calculated within patients that did not receive Tamoxifen monotherapy.

## Multivariate Cox Models
Covariates consist of age (continuous), node status (binary, 0,1), grade (semi-continuous, 1,2,3), and tumour size (continuous). PAM50 subtype is coded as a categorical variable with Normal-like as the baseline.

- [ER IHC Status - RFS](../GSE6532/Cox/Untreated/Multi_Cox_ER_IHC.csv)
- [ESR1 mRNA tertile expression - RFS](../GSE6532/Cox/Untreated/Multi_Cox_ESR1_TERT.csv)
- [ESR1 mRNA tertile expression + PAM50 - RFS](../GSE6532/Cox/Untreated/Multi_Cox_ESR1_TERT_PAM50.csv)
- [TNFRSF12A mRNA tertile expression - RFS](../GSE6532/Cox/Untreated/Multi_Cox_TNFRSF12A_TERT.csv)
- [TNFRSF12A mRNA tertile expression + PAM50 - RFS](../GSE6532/Cox/Untreated/Multi_Cox_TNFRSF12A_TERT_PAM50.csv)