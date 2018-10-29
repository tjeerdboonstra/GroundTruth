# GroundTruth
Code for analysing correlations between linguistic features and mental health scores of the GroundTruth study

Overview: To investigate the associations between linguistic features and symptoms of depression, generalised anxiety, and suicidal ideation, we extracted linguistic features from individuals’ blog content and correlated it with validated mental health data in a longitudinal study (n=38). Depressive symptoms were assessed using the self-report Patient Health Questionnaire (PHQ-9), anxiety symptoms using the self-report Generalised Anxiety Disorder Scale (GAD-7), and social media data was analysed using the Linguistic Inquiry and Word Count (LIWC) tool for linguistic features.  Bivariate and multivariate analyses were performed to investigate the correlations between the linguistic features and mental health scores between subjects. We then used the multivariate regression model to predict longitudinal changes in mood within subjects. 

Data availability: The linguistic features and mental health scores used in this study are available at Zenodo: ... 

Code: This repository contains the Matlab script to analyse the data and generate the results reported in: ... The following scripts are provided:

performBivariateGroupAnalysis.m: Computes the rank-order correlations between each of the 68 linguistic features and the 3 mental health scores. Correlations are computed between subjects using data averaged across repeated measures. Permutation testing is used to control the family-wise error rate.

performMultivariateGroupAnalysis.m: Performs partial least squares (PLS) regression between multiple linguistic features and the mental health scores at group level. 5-fold cross-validation is used to determine the number of components of the model and prediction accuracy is assessed using Mean Square Error (MSE). PLS regression is performed both on the full and a restricted feature set using the most robust features as determined by bootstrapping.

performLongitudinalAnalysis.m: Tests the PLS regression model constructed on group-level data on repeated measures of single participants. The group-level model is used to predict the mental scores at each time point at which linguistic features were extracted and mental health were assessed. The predicted and observed mental health scores are correlated across time points for each participant and correlation coefficients estimated for each participant are compared at the group level.
