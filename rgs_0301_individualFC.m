%% Gradient analysis - calculate individual FC matrices
% written by Magdalena del Rio and Chris Racey
% 2021

clear all;close all;clc

% define pathways
gradientdatadir = '.\gradients\'; % where to save FC matrix
rawdatadir = '.\dicom\'; % where is raw data saved
fmridir = '.\fmri\'; % where is the timeseries data saved
datafilesdir = '.\datafiles'; %intermediate data files saved here, to bypass custom-made functions

cd(rawdatadir)
fl = dir([rawdatadir, 'GBB_s*']);
fl = {fl.name}';

%define batch
subjectid_batch=cellfun(@(fl) fl(5:end),fl(1:end,1), 'UniformOutput', false);

for pp=1:length(subjectid_batch)
    
    subjectid=subjectid_batch{pp}

    % Import data as timeseries   
    cd('.\scripts_batches\')

    if isfile([fmridir,filesep,subjectid,'\preprocessSURF\fsa_run01.mat'])
        timeseries_fsaverage=load([fmridir,filesep,subjectid,'\preprocessSURF\fsa_run01.mat']);

        % Confound regression (load in motion traces from PP, then regress out)
        load([fmridir,filesep,subjectid,'\preprocessSURFfigures\record.mat'])

        timeseries_fsaverage.data=DetrendDataByDesign(double(timeseries_fsaverage.data),mparams{1}(2:end,:));

        % Deal with atlas parcelation
        atlasstr='Schaefer2018_400Parcels_7Networks_order';
        hemis={'lh','rh'};
        surfsuffix='orig';
%         [roimask,roinames]=cvnroimask('fsaverage',hemis,atlasstr,[],surfsuffix,'collapsevals'); %extract atlas indices and names for a particular atlas and particular subject (in this case fsaverage)
        load(fullfile(datafilesdir, 'roimask.mat')); load(fullfile(datafilesdir, 'roinames.mat'));        
        roimask=catcell(1,roimask); %concatenate hemisphere roi mask indices into whole brain vector
        roinames=regexprep(roinames{1},['@',atlasstr],''); %strip superfluous information from roinames (labels) (chop off '@atlasname')
        
        %pre-assign empty variables of various types
        roimask_lhrh={};
        roimask_lh={};
        roimask_rh={};
        roinames_lh={};
        roinames_rh={};

        %Make a unique binary mask for each atlas roi (separate and joined hemis)
        for r=1:numel(roinames)
            roimask_lhrh{r}=roimask==r; % find indices belonging to each atlas roi
            roimask_lh{r}=roimask_lhrh{r}(1:timeseries_fsaverage.numlh); % split into hemispheres
            roimask_rh{r}=roimask_lhrh{r}(timeseries_fsaverage.numlh+1:end);
            roinames_lh{r}=sprintf('%s.%s',hemis{1},roinames{r}); % corresponding atlas label
            roinames_rh{r}=sprintf('%s.%s',hemis{2},roinames{r});
        end

        % get rid of the background+freesurfer defined medial wall ROI
        bgIdx = find(contains(roinames,'Background'));
        
        roimask_lh{bgIdx}=[]; roimask_lh=roimask_lh(~cellfun('isempty',roimask_lh));
        roimask_rh{bgIdx}=[]; roimask_rh=roimask_rh(~cellfun('isempty',roimask_rh));
        roinames_lh{bgIdx}=[]; roinames_lh=roinames_lh(~cellfun('isempty',roinames_lh));
        roinames_rh{bgIdx}=[]; roinames_rh=roinames_rh(~cellfun('isempty',roinames_rh));
        roinames{bgIdx}=[]; roinames=roinames(~cellfun('isempty',roinames));

        roimask_wholebrain=[roimask_lh,roimask_rh]; % and don't split into hemispheres
        roinames_wholebrain=[roinames_lh,roinames_rh];

        lhdata=timeseries_fsaverage.data(:,1:timeseries_fsaverage.numlh); % split timeseries by hemisphere
        rhdata=timeseries_fsaverage.data(:,timeseries_fsaverage.numlh+1:end);

        % Downsample to roi regions: mean timeseries for each label
        lhsubsample=zeros(size(lhdata,1),numel(roinames));
        rhsubsample=zeros(size(rhdata,1),numel(roinames));

        for r=1:numel(roinames)%loop over all rois
            lhsubsample(:,r)=mean(lhdata(:,roimask_lh{r}),2);%index the data using the atlas mask and take the average of all vertices in the region
            rhsubsample(:,r)=mean(rhdata(:,roimask_rh{r}),2);%index the data using the atlas mask and take the average of all vertices in the region
        end
        why

        wholebrainsubsample = [lhsubsample, rhsubsample];


        % Correlation matrix
        close all
        tic
        lh_corrmat=corr(lhsubsample);
        rh_corrmat=corr(rhsubsample);
        wholebrain_corrmat=corr(wholebrainsubsample);

        % Save correlation matrices
        cd(gradientdatadir);
        mkdirquiet(subjectid);
        cd([gradientdatadir,filesep,subjectid]);

        save(fullfile([gradientdatadir,filesep,subjectid],'lh_corrmat_schaefer400_7n.mat'), 'lh_corrmat');
        save(fullfile([gradientdatadir,filesep,subjectid],'rh_corrmat_schaefer400_7n.mat'), 'rh_corrmat');
        save(fullfile([gradientdatadir,filesep,subjectid],'lhrh_corrmat_schaefer400_7n.mat'), 'wholebrain_corrmat');

    end

end

why
