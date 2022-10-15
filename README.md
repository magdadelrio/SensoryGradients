# Sensory Gradients


The main analysis script is 'allDistanceMetrics.m', which computes various integration/segregation metrics on the basis of the gradient decomposition results. These are found in 'allgradients_n370.mat', which contains a 3D matrix with the first 10 gradients per subject using the Schaefer atlas with 400 parcels (dimensions are number of subjects x number of parcels x number of gradients). The main script also contains code to plot the correlations with the AQ and GSQ. The main questionnaire scores (AQ total and subscale scores, GSQ total scores), age and gender are in a table in 'alldata_n370.mat'. The table in 'alldata_subscales_n370.mat' additionally contains all subscales for the GSQ, divided by sensory modality as well as hypo-/hypersensitivity. The files 'roimask.mat' and 'roinames.mat' are needed for the parcellation of the gradient scores.

The scripts 'rgs_0301_individualFC.m', 'rgs_0302_groupFC.m' and 'rgs_031_gradients.m' show the pipeline used to derive the gradient decomposition results from the timeseries data (unfortunately not allowed to share).

The script 'rgs_032_gradients_similarity.m' calculates the similarity between each group-level gradient and each individual-level gradient for gradients 1:3, and saves the results in 'mvrdata.mat'. Requires 'avg_corrmat_schaefer400_7n.mat', which contains the group-average correlation matrix (dimensions 400 x 400).

The markdown file 'gradientAnalysis.Rmd' contains code to:
1. plot the distribution of the median gradient scores segregated by the 7 networks of the Schaefer atlas for participants with high and low AQ/GSQ scores as per a median split
2. run the multivariate regression model based on gradient similarity
3. run a multiple regression model predicting the distance between the median gradient scores for default and visual networks from age, gender, AQ and GSQ
4. plot schematic illustrations of several integration/segregation metrics.

Several additional files are required to run the script without modifications:
1. 'mvrdata.mat' contains the gradient similarity scores for gradients 1:3
2. 'mediandata.mat' contains the questionnaire scores and distances between medians for each network combination
3. 'roinetworks_lhrh.mat' contains the parcellation scheme.
