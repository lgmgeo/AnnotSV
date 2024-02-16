# SV ranking/classification

In order to assist the clinical interpretation of SV, AnnotSV provides on top of the annotations a ranking score to assess SV pathogenicity. This score is an adaptation of the work provided by the joint consensus recommendation of ACMG and ClinGen ([Riggs et al., 2020](https://www.nature.com/articles/s41436-019-0686-8)). We especially took attention in scoring as much as possible recessive SV observed in various dataset (NGS, array based...)

<br />

| Score | Class type | Class|
| :---: | :--- | :--- |
|   0.99 or more points | Pathogenic | Class 5 |
| 0.90 to 0.98 points   | Likely pathogenic	| Class 4 |
| 0.89 to −0.89 points | Variant of uncertain significance | Class 3 |
| −0.90 to −0.98 points | Likely benign	| Class 2 |
| −0.99 or fewer points | Benign | Class 1 |

<br />

# Method

The comprehensive and detailed scoring guidelines are available in the [Scoring_Criteria_AnnotSV](../../../Scoring_Criteria_AnnotSV_v3.4.xlsx) file (see Table1 for loss SV and Table2 for gain SV).
