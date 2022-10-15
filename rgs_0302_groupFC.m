%% Gradient analysis - calculate mean FC matrix
% written by Magdalena del Rio
% 2021

clear all;close all;clc

% define pathways
gradientdatadir = '.\gradients\'; % where to save FC matrix
rawdatadir = '.\dicom\'; % where is raw data saved
fmridir = '.\fmri\'; % where is the timeseries data saved
ppdir = '.\ppresults\'; % where is the timeseries data saved

cd(rawdatadir)
fl = dir([rawdatadir, 'GBB_s*']);
fl = {fl.name}';

%define batch
subjectid_batch=cellfun(@(fl) fl(5:end),fl(1:end,1), 'UniformOutput', false);

lhrh_corrmat_all=zeros(400,400,length(subjectid_batch)); lhrh_corrmat_all_z=zeros(400,400,length(subjectid_batch));

for pp=1:length(subjectid_batch)

    subjectid=subjectid_batch{pp};
    corrmat=load([gradientdatadir,filesep,subjectid,filesep,'lhrh_corrmat_schaefer400_7n.mat']);
    corrmat=cell2mat(struct2cell(corrmat));
    corrmat_z=atanh(corrmat); %z transform

    %put into a nnodes x nnodes x n matrix
    corrmat_all_z(:,:,pp)=corrmat_z;
end

% exclusions:
load([ppdir,filesep,'\allsubjects\toexclude.mat']);
corrmat_all_z_clean = corrmat_all_z;
subjectid=cellfun(@(subjectid_batch) subjectid_batch(2:end),subjectid_batch(1:end,1), 'UniformOutput', false);
subjectid=str2double(subjectid);
toexclude_index=find(ismember(subjectid,toexclude));
toexclude_index=sort(toexclude_index,'desc');

for ii=1:length(toexclude_index)
    corrmat_all_z_clean(:,:,toexclude_index(ii))=[];
end

[r,c,v] = ind2sub(size(corrmat_all_z_clean),find(isnan(corrmat_all_z_clean)));
toexclude2=sort(unique(v),'desc');

for ii=1:length(toexclude2)
    corrmat_all_z_clean(:,:,toexclude2(ii))=[];
end

avgcorrmat_z=mean(corrmat_all_z_clean,3);

corrmat=tanh(avgcorrmat_z); %z to r-transform   

save(fullfile([gradientdatadir,filesep,'avg'],'lhrh_corrmat_schaefer400_7n_clean.mat'), 'corrmat');
 
