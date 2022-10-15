%% Multivariate similarity analysis
% spearman rank correlation between each group-level gradient and each individual-level gradient for each participant,
% for gradients 1:3, resulting in 3 similarity scores per participant

clear all;close all;clc
addpath(genpath('.\BrainSpace-0.1.1\'))

load('allgradients_n370.mat')
load('alldata_n370.mat')
load('avg_corrmat_schaefer400_7n.mat');

Grgs10 = GradientMaps('kernel', 'cosine similarity', 'approach', 'pca', 'n_components', 10); 
Grgs10 = Grgs10.fit(avgcorrmat_r, 'sparsity', 95);

g_group = Grgs10.gradients{1};
g_similarity_all=nan(size(g_aligned_all,1),3); g_similarity_all_z=nan(size(g_aligned_all,1),3);

for pp=1:size(g_aligned_all,1)
    for gg=1:3
        g_similarity_all(pp,gg) = corr(g_group(:,gg), g_aligned_all(pp,:,gg)', 'Type', 'Spearman');
        g_similarity_all_z(pp,gg) = atanh(corr(g_group(:,gg), g_aligned_all(pp,:,gg)', 'Type', 'Spearman'));

    end
end

simTab=splitvars(table(g_similarity_all_z));
simTab.Properties.VariableNames={'simGradient1', 'simGradient2', 'simGradient3'};

mvrDataTab=[data_final, simTab];
mvrData=table2array(mvrDataTab);
% save('mvrData.mat','mvrData');

