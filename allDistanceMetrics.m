%% Distance metrics
% Magdalena del Rio, University of Sussex
% January 2022

clear all;close all;clc
cd('C:\Users\Magda\Documents\Matlab\rgs\n370');
addpath('C:\Users\Magda\Documents\Matlab\rgs\datafiles');
load('allgradients_n370.mat') % g_aligned_all: matrix with first 10 gradients per subject using Schaefer 400 (size nSub x nParcel x nGradients)
load('alldata_n370.mat') % data_final: table with subject ID, gender, age and questionnaire scores

%% Network median distance

% get data for all subjects per network
% load in labels
hemis={'lh','rh'};
atlasstr='Schaefer2018_400Parcels_7Networks_order'
surfsuffix='orig';
load('roimask.mat'); load('roinames.mat');
roimask_lh=roimask{1}; roimask_rh=roimask{2};

roinames=regexprep(roinames{1},['@',atlasstr],''); %strip superfluous information from roinames (labels) (chop off '@atlasname')

bgIdx = find(contains(roinames,'Background'));
roinames{bgIdx}=[]; roinames=roinames(~cellfun('isempty',roinames));

roinetworks=extractBetween(roinames, "LH_", "_");
roinetworks_lhrh=[roinetworks,roinetworks];
roinetworks_list=unique(extractBetween(roinames, "LH_", "_"));

for ii=1:length(roinetworks_list)
    roinetworks_lhrh_idx=find(strcmp(roinetworks_lhrh,roinetworks_list{ii}));
    roinetworks_sizes(ii)=numel(roinetworks_lhrh_idx);
    roinetworks_fillers{ii}=zeros(size(g_aligned_all,1),roinetworks_sizes(ii))
end

args=[roinetworks_list;roinetworks_fillers];
g_all_nw = struct(args{:});

for ii=1:length(roinetworks_list)
    roinetworks_lhrh_idx=find(strcmp(roinetworks_lhrh,roinetworks_list{ii}));
    thisNw=roinetworks_list{ii};
    for p=1:size(g_aligned_all,1)
        g_all_nw.(thisNw)(p,:)=g_aligned_all(p,roinetworks_lhrh_idx,1);
    end
end


% calculate median for each network per subject
for ii=1:length(roinetworks_list)
    thisNw=roinetworks_list{ii};
    for p=1:size(g_all_nw.(thisNw),1)
        median_all_nw(p,ii) = median(g_all_nw.(thisNw)(p,:));
    end
end

median_data_aqgsq = [data_final, splitvars(table(median_all_nw))];
median_data_aqgsq.Properties.VariableNames(end-6:end) = lower(cellfun(@(roinetworks_list)[roinetworks_list '_median'],roinetworks_list,'uni',false));

% planned comparison: correlation of AQ and GSQ with distance default-vis
median_data_aqgsq.defaultvis_median_diff=abs(median_data_aqgsq.default_median-median_data_aqgsq.vis_median);

% correlate distance with AQ and GSQ
[rho,pval]=corr(median_data_aqgsq.defaultvis_median_diff, median_data_aqgsq.AQtotal, 'Type', 'Spearman') 
[rho,pval]=corr(median_data_aqgsq.defaultvis_median_diff, median_data_aqgsq.GSQtotal, 'Type', 'Spearman')

% plot (Figure 4)
% 4 A (AQ)
figure
scatter(median_data_aqgsq.defaultvis_median_diff, median_data_aqgsq.AQtotal, 'filled',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlabel({'Difference in median gradient scores';'for visual and default networks'}); ylabel('AQ');
xlim([0 5]); ylim([0 50]);
axis square
hold on
p = polyfit(median_data_aqgsq.defaultvis_median_diff, median_data_aqgsq.AQtotal, 1);
px = [min(median_data_aqgsq.defaultvis_median_diff) max(median_data_aqgsq.defaultvis_median_diff)];
py = polyval(p, px);
plot(px, py, 'b', 'LineWidth', 2);
set(gcf,'color','w'); set(gca,'FontSize',12, 'FontName', 'Arial')

% 4 A (GSQ)
figure
scatter(median_data_aqgsq.defaultvis_median_diff, median_data_aqgsq.GSQtotal, 'filled',...
        'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlabel({'Difference in median gradient scores';'for visual and default networks'}); ylabel('GSQ');
xlim([0 5]); ylim([0 120]);
axis square
hold on
p = polyfit(median_data_aqgsq.defaultvis_median_diff, median_data_aqgsq.GSQtotal, 1);
px = [min(median_data_aqgsq.defaultvis_median_diff) max(median_data_aqgsq.defaultvis_median_diff)];
py = polyval(p, px);
plot(px, py, 'b', 'LineWidth', 2);
set(gcf,'color','w'); set(gca,'FontSize',12, 'FontName', 'Arial')


% follow-up correlations: correlate medians with GSQ
[rho,pval]=corr(median_data_aqgsq.default_median, median_data_aqgsq.GSQtotal, 'Type', 'Spearman') 
[rho,pval]=corr(median_data_aqgsq.vis_median, median_data_aqgsq.GSQtotal, 'Type', 'Spearman') 


% exploratory: distance between all other network combinations and AQ, GSQ
for n1 = 1:length(roinetworks_list)
    thisNw1=roinetworks_list{n1};
    for n2 = 1:length(roinetworks_list)
        thisNw2=roinetworks_list{n2};
        combilabels{n1,n2}=strcat(thisNw1,thisNw2);
        combis{n1,n2}={thisNw1,thisNw2};
    end
end

tmp = magic(size(combis));
tmp = tmp - diag(diag(tmp));

allcombilabels = combilabels(logical(tril(tmp)));
allcombis = combis(logical(tril(tmp)));
% median_data_aqgsq.Properties.VariableNames(end-6:end)=roinetworks_list;

for nn = 1:length(allcombis)   
    all_d_aqgsq(:,nn)=abs(eval(['median_data_aqgsq.',lower(allcombis{nn}{1}),'_median'])-eval(['median_data_aqgsq.',lower(allcombis{nn}{2}),'_median']));
    [rho_daq(nn),pval_daq(nn)]=corr(all_d_aqgsq(:,nn), median_data_aqgsq.AQtotal, 'Type', 'Spearman');
    [rho_dgsq(nn),pval_dgsq(nn)]=corr(all_d_aqgsq(:,nn), median_data_aqgsq.GSQtotal, 'Type', 'Spearman');
end

dCorrTab = table(allcombilabels, rho_daq', pval_daq', rho_dgsq', pval_dgsq');
dCorrTab.Properties.VariableNames={'Network combination', 'rho_AQ', 'p_AQ', 'rho_GSQ', 'p_GSQ'}

% exploratory: all other network medians with AQ, GSQ
for n = 1:length(roinetworks_list)   
    all_m_aqgsq(:,n)=eval(['median_data_aqgsq.',lower(roinetworks_list{n}),'_median']);
    [rho_maq(n),pval_maq(n)]=corr(all_m_aqgsq(:,n), median_data_aqgsq.AQtotal, 'Type', 'Spearman');
    [rho_mgsq(n),pval_mgsq(n)]=corr(all_m_aqgsq(:,n), median_data_aqgsq.GSQtotal, 'Type', 'Spearman');
end

mCorrTab = table(roinetworks_list', rho_maq', pval_maq', rho_mgsq', pval_mgsq');
mCorrTab.Properties.VariableNames={'Network', 'rho_AQ', 'p_AQ', 'rho_GSQ', 'p_GSQ'}


%% Network peak distance

% calculate peak for each network per subject
for ii=1:length(roinetworks_list)
    thisNw=roinetworks_list{ii};
    for p=1:size(g_all_nw.(thisNw),1)
        [f, xi] = ksdensity(g_all_nw.(thisNw)(p,:));% create histogram
        [pks, locs] = findpeaks(f); % find peaks
        peak_all_nw(p,ii)= double(pks(1)); % find height of first and last peak
    end
end


% calculate distance between peaks of interest (visual and default)
peakDistanceVisDefault = abs( peak_all_nw(:,2)- peak_all_nw(:,7))

% correlate with distance metric 
[rho,pval]=corr(abs(median_data_aqgsq.vis_median), peak_all_nw(:,7), 'Type', 'Spearman') 
[rho,pval]=corr(abs(median_data_aqgsq.default_median), peak_all_nw(:,2), 'Type', 'Spearman') 

[rho,pval]=corr(median_data_aqgsq.defaultvis_median_diff, peakDistanceVisDefault, 'Type', 'Spearman') 

% correlate with AQ, GSQ
[rho,pval]=corr(peakDistanceVisDefault, data_final.AQtotal, 'Type', 'Spearman') 
[rho,pval]=corr(peakDistanceVisDefault, data_final.GSQtotal, 'Type', 'Spearman') 
    


%% Gradient range and variation
% (Concept from Xia et al, 2020)

for g = 1:3
    for p = 1:size(g_aligned_all, 1)
        g_range(p,g) = max(g_aligned_all(p,:,g)) - min(g_aligned_all(p,:,g));
        g_variance(p,g) = std(g_aligned_all(p,:,g));
    end
end

% correspondence across metrics: distance vs range and variance
[rho,pval]=corr(median_data_aqgsq.defaultvis_median_diff, g_range(:,1), 'Type', 'Spearman') 
[rho,pval]=corr(median_data_aqgsq.defaultvis_median_diff, g_variance(:,1), 'Type', 'Spearman')

% AQ/GSQ correlations
[rho,pval]=corr(median_data_aqgsq.GSQtotal, g_range(:,1), 'Type', 'Spearman') 
[rho,pval]=corr(median_data_aqgsq.AQtotal, g_range(:,1), 'Type', 'Spearman')

[rho,pval]=corr(median_data_aqgsq.GSQtotal, g_variance(:,1), 'Type', 'Spearman') 
[rho,pval]=corr(median_data_aqgsq.AQtotal, g_variance(:,1), 'Type', 'Spearman') 

%% Eccentricity
% adapted from https://github.com/CNG-LAB/cngopen/tree/main/social_gradients (Valk et al, 2021)

% for each vertex, calculated its distance from the center of the gradient coordinate system formed by G1, G2, and G3
% for each individual. This eccentricity captured vertex-wise intrinsic functional integration (low eccentricity)
% and segregation (high eccentricity) in a single scalar value

for p=1:size(g_aligned_all,1)
    GGG(p,:) = sqrt((g_aligned_all(p,:,1).^2)+(g_aligned_all(p,:,2).^2)+(g_aligned_all(p,:,3).^2));
end

% average global eccentricity 
GGG_mean = mean(GGG(:,:),2);

% eccentricity by network
for ii=1:length(roinetworks_list)
    roinetworks_lhrh_idx=find(strcmp(roinetworks_lhrh,roinetworks_list{ii}));
    roinetworks_sizes(ii)=numel(roinetworks_lhrh_idx);
    roinetworks_fillers{ii}=zeros(size(g_aligned_all,1),roinetworks_sizes(ii))
end

args=[roinetworks_list;roinetworks_fillers];
g1_all_nw = struct(args{:});
g2_all_nw = struct(args{:});
g3_all_nw = struct(args{:});

for ii=1:length(roinetworks_list)
    roinetworks_lhrh_idx=find(strcmp(roinetworks_lhrh,roinetworks_list{ii}));
    thisNw=roinetworks_list{ii};
    for p=1:size(g_aligned_all,1)
        g1_all_nw.(thisNw)(p,:)=g_aligned_all(p,roinetworks_lhrh_idx,1);
        g2_all_nw.(thisNw)(p,:)=g_aligned_all(p,roinetworks_lhrh_idx,2);
        g3_all_nw.(thisNw)(p,:)=g_aligned_all(p,roinetworks_lhrh_idx,3);
    end
end


for ii=1:length(roinetworks_list)   
    thisNw=roinetworks_list{ii};
    for p=1:size(g_aligned_all,1)
        GGG_nw{ii}(p,:) = sqrt((g1_all_nw.(thisNw)(p,:).^2)+(g2_all_nw.(thisNw)(p,:).^2)+(g3_all_nw.(thisNw)(p,:).^2));
        GGG_nw_mean(p,ii) = mean(GGG_nw{ii}(p,:)); % average eccentricity by network
    end    
end

GGG_aqgsq = [data_final, splitvars(table(GGG_nw_mean))];
GGG_aqgsq.Properties.VariableNames(end-6:end) = cellfun(@(roinetworks_list)[roinetworks_list '_GGG'],roinetworks_list,'uni',false);

[rho_meanGSQ,pval_meanGSQ]=corr(GGG_mean, GGG_aqgsq.GSQtotal, 'Type', 'Spearman');
[rho_meanAQ,pval_meanAQ]=corr(GGG_mean, GGG_aqgsq.AQtotal, 'Type', 'Spearman');

corrMeanTab = table(rho_meanAQ', pval_meanAQ', rho_meanGSQ', pval_meanGSQ');
corrMeanTab.Properties.VariableNames = {'rho_AQ', 'pval_AQ', 'rho_GSQ', 'pval_GSQ'}
corrMeanTab.Properties.RowNames = {'Global'};

for ii=1:length(roinetworks_list)
    [rho_GSQ(ii),pval_GSQ(ii)]=corr(table2array(GGG_aqgsq(:,end-length(roinetworks_list)+ii)), GGG_aqgsq.GSQtotal, 'Type', 'Spearman');
    [rho_AQ(ii),pval_AQ(ii)]=corr(table2array(GGG_aqgsq(:,end-length(roinetworks_list)+ii)), GGG_aqgsq.AQtotal, 'Type', 'Spearman');

end

corrTab = table(rho_AQ', pval_AQ', rho_GSQ', pval_GSQ');
corrTab.Properties.VariableNames = {'rho_AQ', 'pval_AQ', 'rho_GSQ', 'pval_GSQ'}
corrTab.Properties.RowNames = roinetworks_list
corrTab=[corrMeanTab;corrTab]

%% Dispersion
% adapted from https://github.com/rb643/GradientDispersion (Bethlehem et al, 2020)

% each axis of this 3D space is defined by the values along the first three gradients. Within network dispersion is 
% quantified as sum squared Euclidean distance of network nodes to the network centroid at individual level. Between
% network dispersion is calculated as the Euclidean distance between network centroids.

% within-network dispersion (pdist2 between each network node and and the network median)
dist_to_center=0;

for n=1:length(roinetworks_list)   
    thisNw=roinetworks_list{n};
    for p=1:size(g_aligned_all,1)
        
            median_nw{n}(p,:) = median([g1_all_nw.(thisNw)(p,:);g2_all_nw.(thisNw)(p,:);g3_all_nw.(thisNw)(p,:)],2);
            clear dist_to_center
            for ii = 1:length(g1_all_nw.(thisNw)(p,:))
                dist_to_center(ii) = pdist2([g1_all_nw.(thisNw)(p,ii);g2_all_nw.(thisNw)(p,ii);g3_all_nw.(thisNw)(p,ii)]',...
                    median_nw{n}(p,:));
            end
            yeo_netSSD(p,n) = sumsqr(dist_to_center); % global metric
    end    
end

dispersion_aqgsq = [data_final, splitvars(table(yeo_netSSD))];
dispersion_aqgsq.Properties.VariableNames(end-6:end) = strcat('disp',roinetworks_list)


for ii=1:length(roinetworks_list)
    [rho_GSQ(ii),pval_GSQ(ii)]=corr(table2array(dispersion_aqgsq(:,end-length(roinetworks_list)+ii)), dispersion_aqgsq.GSQtotal, 'Type', 'Spearman');
    [rho_AQ(ii),pval_AQ(ii)]=corr(table2array(dispersion_aqgsq(:,end-length(roinetworks_list)+ii)), dispersion_aqgsq.AQtotal, 'Type', 'Spearman');

end

corrTab = table(rho_AQ', pval_AQ', rho_GSQ', pval_GSQ');
corrTab.Properties.VariableNames = {'rho_AQ', 'pval_AQ', 'rho_GSQ', 'pval_GSQ'};

corrTab.Properties.RowNames = roinetworks_list;

% between-network dispersion - ignoring overlap calculations
net_dist = zeros(7,7,size(g_aligned_all,1));
for n1 = 1:length(roinetworks_list)
    thisNw1=roinetworks_list{n1};
    for n2 = 1:length(roinetworks_list)
        thisNw2=roinetworks_list{n2};
        net_dist(n1,n2,:) = sqrt(sum((median_nw{n1}(:,:)' - median_nw{n2}(:,:)').^2));
    end
end

for s = 1:size(g_aligned_all,1)
    net_dist_long(s,:) = squareform(net_dist(:,:,s));
end

for n1=1:length(roinetworks_list)
    for n2=1:length(roinetworks_list)
        combi(n1,n2)=strcat(roinetworks_list(n1),roinetworks_list(n2))
    end
end


tmp = magic(size(combi));
tmp = tmp - diag(diag(tmp));

allcombis = combi(logical(tril(tmp)));

dispersion_between_aqgsq = [data_final, splitvars(table(net_dist_long))];
dispersion_between_aqgsq.Properties.VariableNames(end-20:end) = strcat('disp_',allcombis)

for ii=1:length(allcombis)
    [rho_dispbetween_GSQ(ii),pval_dispbetween_GSQ(ii)]=corr(table2array(dispersion_between_aqgsq(:,end-length(allcombis)+ii)), dispersion_between_aqgsq.GSQtotal, 'Type', 'Spearman');
    [rho_dispbetween_AQ(ii),pval_dispbetween_AQ(ii)]=corr(table2array(dispersion_between_aqgsq(:,end-length(allcombis)+ii)), dispersion_between_aqgsq.AQtotal, 'Type', 'Spearman');

end

corrTab = table(rho_dispbetween_AQ', pval_dispbetween_AQ', rho_dispbetween_GSQ', pval_dispbetween_GSQ');
corrTab.Properties.VariableNames = {'rho_AQ', 'pval_AQ', 'rho_GSQ', 'pval_GSQ'};

corrTab.Properties.RowNames = allcombis;

% correspondence across metrics: distance vs between network dispersion
[rho,pval]=corr(median_data_aqgsq.defaultvis_median_diff, dispersion_between_aqgsq.disp_VisDefault, 'Type', 'Spearman') 


%% Explained variance - control test
load('alllambdas_n370.mat') % g_lambda_all contains the explained variance for each gradient and subject (size nSub x nGradients)

T = [data_final,splitvars(table(g_lambda_all(:,1:3)))];
T.Properties.VariableNames(end-2:end)={'Lambda1', 'Lambda2', 'Lambda3'};


for nn = 1:3   
    [rho_lambdaaq(nn),pval_lambdaaq(nn)]=corr(eval(['T.Lambda',num2str(nn)]), T.AQtotal, 'Type', 'Spearman');
    [rho_lambdagsq(nn),pval_lambdagsq(nn)]=corr(eval(['T.Lambda',num2str(nn)]), T.GSQtotal, 'Type', 'Spearman');
end


