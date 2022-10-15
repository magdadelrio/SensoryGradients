%% Gradient derivation and alignment
% written by Magdalena del Rio
% 2021

clear all;close all;clc

% define pathways
gradientdatadir = '.\gradients\'; % FC matrices and gradients
addpath(genpath('.\BrainSpace-0.1.1\'))
datafilesdir = '.\datafiles'; %intermediate data files saved here, to bypass custom-made functions
cd(gradientdatadir)
fl = dir([gradientdatadir, 's*']);
fl = {fl.name}';

%define batch
subjectid_batch=cellfun(@(fl) fl(2:end),fl(1:end,1), 'UniformOutput', false);

% exclusions:
load([gradientdatadir,filesep,'\avg\allexclusions.mat']);
subjectid=str2double(subjectid_batch);
toexclude_index=find(ismember(subjectid,allexclusions));
toexclude_index=sort(toexclude_index,'desc');

for ii=1:length(toexclude_index)
    fl(toexclude_index(ii))=[];
end


% Gradient derivation from HCP data
conn_matrix_hcp = load_group_fc('schaefer',400);
conn_matrix_hcp = conn_matrix_hcp.schaefer_400;

Ghcp = GradientMaps('kernel', 'cosine similarity', 'approach', 'pca'); % kernel default normalized angle, approach default diffusion embedding,
% alignment default none, random state default none, n_componenents default 10 as arguments
Ghcp = Ghcp.fit(conn_matrix_hcp, 'sparsity', 95);

% Gradient derivation from our mean data
conn_matrix_rgs = load([gradientdatadir,filesep,'avg',filesep,'lhrh_corrmat_schaefer400_7n_clean.mat']);
conn_matrix_rgs=cell2mat(struct2cell(conn_matrix_rgs));

% Construct the reference gradient
Gref = GradientMaps('kernel', 'cosine similarity', 'approach', 'pca'); % kernel default normalized angle, approach default diffusion embedding,
% alignment default none, random state default none, n_componenents default 10 as arguments
Gref = Gref.fit(conn_matrix_rgs, 'sparsity', 95);

% Align each subject to the group average
for ii=1:numel(fl)
    subjectid=fl{ii};
    conn_matrix_sub = load([gradientdatadir,filesep,subjectid,filesep,'lhrh_corrmat_schaefer400_7n.mat']);
    conn_matrix_sub=cell2mat(struct2cell(conn_matrix_sub));
    

    Gp = GradientMaps('kernel', 'cosine similarity', 'approach', 'pca', 'alignment', 'pa');
    Gp = Gp.fit(conn_matrix_sub,'reference',Gref.aligned{1}, 'sparsity', 95);
    
    %save un/aligned gradients
    gradients_unaligned = Gp.gradients{1};
    gradients_aligned = Gp.aligned{1};
    
    save(fullfile([gradientdatadir,filesep,subjectid],'gradients_unaligned_schaefer400_7n.mat'), 'gradients_unaligned');
    save(fullfile([gradientdatadir,filesep,subjectid],'gradients_aligned_schaefer400_7n.mat'), 'gradients_aligned');

end


why


%% Plot 
% Figure 1 and Supplementary Figure 1
hemis={'lh','rh'};
atlasstr='Schaefer2018_400Parcels_7Networks_order'
surfsuffix='orig';
load(fullfile(datafilesdir, 'roimask.mat')); load(fullfile(datafilesdir, 'roinames.mat'));
roimask_lh=roimask{1}; roimask_rh=roimask{2};

roinames=regexprep(roinames{1},['@',atlasstr],''); %strip superfluous information from roinames (labels) (chop off '@atlasname')

bgIdx = find(contains(roinames,'Background'));
roinames{bgIdx}=[]; roinames=roinames(~cellfun('isempty',roinames));
roimask_lh=roimask_lh-1;
roimask_rh=roimask_rh-1;


roimask_rh=roimask_rh+200;
roimask_rh(roimask_rh==200)=0;

roimask=[roimask_lh;roimask_rh];

labeling = roimask;

load(fullfile(datafilesdir, 'surf_lh.mat')); load(fullfile(datafilesdir, 'surf_rh.mat'))

surf_lh = convert_surface(surf_lh); surf_rh = convert_surface(surf_rh);

% Plot gradients separately (note that the colormap is flipped within plot_hemispheres by using flipud(parula(256)))
fig = plot_hemispheres(Gref.gradients{1,1}(:,1:3),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'Gradient 1','Gradient 2', 'Gradient 3'});
                           
h=gca
print('-painters','-dsvg','fig1a')
exportgraphics(h, 'fig1a.svg','ContentType','vector')

% Scree plot: how many gradients do we want?         
lambdas=Gref.lambda{1};

h.figure = figure('Color','White');
h.axes = axes(); 
h.plot = plot(lambdas(1:30) ./ sum(lambdas),'o--','Color','k','MarkerFaceColor', 'k');
xlabel('Component number'); ylabel('Scaled eigenvalues');
ylim([0 0.18]);
set(h.axes,'box','off','FontName','DroidSans','FontSize',12)


% Plot more than one gradient on the brain
gradient_in_euclidean(Gref.gradients{1,1}(:,1:2),{surf_lh,surf_rh},labeling);
print('-painters','-dsvg','myVectorFile') 

%% Correlation with HCP data
% Mean similarity was calculated by averaging all participant’s R-to-Z transformed Spearman Rank correlation
% coefficients for each respective gradient.

%load individual aligned gradients
cd(gradientdatadir);

% load individual aligned gradients into n x 400 x 10 matrix
for ii=1:numel(fl)
    subjectid=fl{ii};
    load([gradientdatadir,filesep,subjectid,filesep,'gradients_aligned_schaefer400_7n.mat']);
    g_aligned_all(ii,:,:)=gradients_aligned(:,:);

end

% calculate (z-transformed) correlation between group average gradient 1:3 from he HCP data
% and each of the individual aligned gradients 1:3
for ii=1:size(g_aligned_all,1)
    for g=1:3
        rho(ii,g)=(corr(Ghcp.gradients{1}(:,g),g_aligned_all(ii,:,g)', 'Type','Spearman'));
        rho_z(ii,g)=atanh(corr(Ghcp.gradients{1}(:,g),g_aligned_all(ii,:,g)', 'Type','Spearman'));
    end
end

% get minimum, maximum, mean and SD values for the correlations
for g=1:3
    minrho(g)=min(abs(rho_z(:,g)));
    maxrho(g)=max(abs(rho_z(:,g)));
    meanrho(g)=mean(abs(rho_z(:,g)));
    sdrho(g)=std(rho_z(:,g));
end

% put into a table
rho_individual_z_tab=splitvars(table([minrho; maxrho; meanrho; sdrho]))
rho_individual_z_tab.Properties.VariableNames={'Gradient 1', 'Gradient 2', 'Gradient 3'}
rho_individual_z_tab.Properties.RowNames={'Minimum', 'Maximum', 'Mean', 'SD'}


% calculate (z-transformed) correlation between group average gradient 1:3 from he HCP data
% and the group average aligned gradients 1:3
for g=1:3
    rho_group(g)=(corr(Ghcp.gradients{1}(:,g),Gref.gradients{1}(:,g), 'Type','Spearman'));
    rho_group_z(g)=atanh(corr(Ghcp.gradients{1}(:,g),Gref.gradients{1}(:,g), 'Type','Spearman'));
end

%  put into a table
rho_group_z_tab=splitvars(table(rho_group_z));
rho_group_z_tab.Properties.VariableNames={'Gradient 1', 'Gradient 2', 'Gradient 3'}

% calculate (z-transformed) correlation between group average gradient 1:3 from he HCP data
% and the group average aligned gradients 1:3
for g=1:3
    rho_group(g)=(corr(Ghcp.gradients{1}(:,g),Gref.aligned{1}(:,g), 'Type','Spearman'));
    rho_group_z(g)=atanh(corr(Ghcp.gradients{1}(:,g),Gref.aligned{1}(:,g), 'Type','Spearman'));
end

%  put into a table
rho_group_z_tab=splitvars(table(rho_group_z));
rho_group_z_tab.Properties.VariableNames={'Gradient 1', 'Gradient 2', 'Gradient 3'}


