# A-new-method-for-detecting-the-skewness-of-data-from-the-five-number-summary
Source codes for "A new method for detecting the skewness of data from the five-number summary"

simulated_meta.R provides the simulated meta-analysis to rigorously evaluate the performance and practical utility of the proposed skewness test in meta-analytic workflows.

realdata.R analyzes the Patient Health Questionnaire-9 dataset, explicitly demonstrating how the skewness test informs subgroup analysis and estimator selection.

density.R provides density plots for six skewed distributions.

type1error.R provides the type I error results of the newly proposed test, compared with the competing methods.

power_largelyskewed.R provides the statistical power of the newly proposed test under six largely skewed distributions, compared with the competing methods.

power_mildlyskewed.R provides the statistical power of the newly proposed test under six mildly skewed distributions, compared with the competing methods.

SkewnessTestS3.R provides the skewness test function that can provide the p-value, independent of the choice of significant level.
