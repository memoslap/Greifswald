%% create struct for files 

spm('defaults','fmri');
spm_jobman('initcfg');

%% instead of deleting files use regressor
%https://www.nitrc.org/forum/message.php?msg_id=20752

% get als subjects in list
%path = "/media/Data03/Studies/MeMoSLAP/SUBJECTS";

%% find subjectfolder
first_level=pwd;
first_level=fileparts(first_level);
first_level=fileparts(first_level);
first_level=fileparts(first_level);
path=fullfile(first_level,'/SPM-HGW-P7SP-clean');
%path = "/media/MeMoSLAP_Subjects/derivatives/preprocessing/SPM-standard";
[SubFolderDir, SubFolderNames] = get_subdir(path);

for subfolder_subject=1:length(SubFolderNames)
      if length(SubFolderNames{subfolder_subject}) == 7 % subfolder must habe this length otherwise wrong folder

          SubFolderNames{subfolder_subject}=SubFolderNames{subfolder_subject};
          SubFolderDir{subfolder_subject}=SubFolderDir{subfolder_subject};
      else
          SubFolderNames{subfolder_subject}={};
          SubFolderDir{subfolder_subject}={};
      end
end

SubFolderNames=SubFolderNames(~cellfun('isempty',SubFolderNames));
SubFolderDir=SubFolderDir(~cellfun('isempty',SubFolderDir));

%filelist.SubFolderNames = SubFolderNames
% get all sessions in list
for subfolder_subject=1:length(SubFolderNames)
    %coregistration_path="/media/Data03/Studies/MeMoSLAP/derivatives/coregistration_coordinates";
    %path_preproc= fullfile(path,"/derivatives/preprocessing/SPM-standard");
        fprintf('Sub folder #%d = %s\n', subfolder_subject, SubFolderNames{subfolder_subject});
        path_preproc = fullfile(SubFolderDir{subfolder_subject}, SubFolderNames{subfolder_subject});
        [SesFolderDir, SesFolderNames] = get_subdir(path_preproc);
        sub_name = erase(SubFolderNames{subfolder_subject},"-");
        filelist.(sub_name).Session_list = SesFolderNames;
    
        %%% loop throught Session for each Subject
        % get all files in every session
        for subfolder_session = 1:length(SesFolderNames)
    
            ses_name = erase(SesFolderNames{subfolder_session},"-");
            path_sub=fullfile(SesFolderDir{subfolder_session}, SesFolderNames{subfolder_session});
            
            filelist.(sub_name).Session.(ses_name).T1 = dir(fullfile(path_sub, '**/sub*mprage_T1w.nii'));
           
            filelist.(sub_name).Session.(ses_name).BOLD = dir(fullfile(path_sub, '**/sub*_task-learning_dir-AP_run-01_bold.nii'));
            filelist.(sub_name).Session.(ses_name).MAT = dir(fullfile(path_sub, '**/sub*_task-learning_dir-AP_run-01_bold.mat'));

            %filelist.(sub_name).Session.(ses_name).Correction='No';
        %end
        end       
end



%% loop through all subjects and create matlabbatch
All_Subjects=fieldnames(filelist);
for sub =1:length(All_Subjects)

%% loop through all sessions

        all_session = {'1','2','3','4'};
        %all_session = {'3','4'}

        for session_num = 1:length(all_session)
            if (isfield(filelist.(All_Subjects{sub}).Session,(['ses',all_session{session_num}]))) 
                if (~isempty(filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD)) && (~isempty(filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).T1))  
                %try
                source_image_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).T1.folder;
                source_image_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).T1.name;
                source_image_path = fullfile(source_image_folder, source_image_name);
                source_image_path_cor=fullfile(source_image_folder, ['r',source_image_name]);
                source_image_path_cor_norm=fullfile(source_image_folder, ['y_r',source_image_name]);
                ref_image_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD.folder;
                ref_image_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD.name;
                ref_image_path = fullfile(ref_image_folder, ref_image_name);



                %if ~exist(fullfile(source_image_path,[SubFolderNames{sub}]),'file')
                if isempty(filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).MAT)

                    fprintf('It is running for session %s for subject %s\n', all_session{session_num},All_Subjects{sub})
                    %try
                    path_outlier=fullfile(ref_image_folder,[All_Subjects{sub},'_ses-',all_session{session_num},'_outlier_list.csv']);

                    number_images=get_images(ref_image_path, ref_image_folder,ref_image_name,[]);

                    matlabbatch{1} = create_batch_realign(number_images);
                    %number_images=get_images(ref_image_path, ref_image_folder,ref_image_name,'r');
                    matlabbatch{2} = create_batch_temp(number_images);
                    
                    matlabbatch{3} = create_batch_cor(ref_image_path, source_image_path);   
                    matlabbatch{4} = create_batch_segm(source_image_path_cor);
                    
                    number_images=get_images(ref_image_path, ref_image_folder,ref_image_name,'a');
                    %path_y_T1_image=source_image_path_cor_norm;
                    path_y_T1_image=source_image_path_cor;
                    matlabbatch{5} = create_batch_normalise(path_y_T1_image,number_images);
                    
                    number_images=get_images(ref_image_path, ref_image_folder,ref_image_name,'wa');
                    matlabbatch{6} = create_batch_smooth(number_images);
                    
                    
                    if ~exist(fullfile(path, SubFolderNames{sub},['ses-',all_session{session_num}],'code'),'dir')
                        mkdir(fullfile(path, SubFolderNames{sub},['ses-',all_session{session_num}],'code'))
                    end


                    %preprocessing_batch = fullfile(path, SubFolderNames{sub},['ses-',all_session{session_num}],'code','preprocessing_batch.mat');
                    preprocessing_batch = fullfile(path,SubFolderNames{sub},['ses-',all_session{session_num}],'code','preprocessing_batch.mat');
                    save(preprocessing_batch,'matlabbatch');
                    spm_jobman('run',matlabbatch);

                    
                else
                fprintf('Corregistration already done for session %s for subject %s\n', all_session{session_num},All_Subjects{sub})    
                end
                end
            else 
                fprintf('No session %s for subject %s\n', all_session{session_num},All_Subjects{sub})
            end
       end
end

%% starting defining function

function [SubFolderDir_func, SubFolderNames_func] = get_subdir(act_directory)
files_sub = dir(act_directory);
files_sub(~[files_sub.isdir])=[];
tf = ismember({files_sub.name}, {'.', '..'});
files_sub(tf) = [];  %remove current and parent directory.
SubFolderNames_func = {files_sub.name};
SubFolderDir_func = {files_sub.folder};
end


function [matlabbatch] = create_batch_cor(ref_image_path, source_image_path)
matlabbatch.spm.spatial.coreg.estwrite.ref = {[ref_image_path,',1']};
matlabbatch.spm.spatial.coreg.estwrite.source = {[source_image_path,',1']};
matlabbatch.spm.spatial.coreg.estwrite.other = {''};
matlabbatch.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
%matlabbatch.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch.spm.spatial.coreg.estwrite.eoptions.sep = [2 1];
matlabbatch.spm.spatial.coreg.estwrite.eoptions.tol = [0.01 0.01 0.01 0.0005 0.0005 0.0005 0.005 0.005 0.005 0.0005 0.0005 0.0005];
%matlabbatch.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch.spm.spatial.coreg.estwrite.roptions.interp = 7;
%matlabbatch.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

end


function [matlabbatch] = create_batch_segm(ref_image_path)

matlabbatch.spm.spatial.preproc.channel.vols = {[ref_image_path,',1']};
%matlabbatch.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate & Reslice: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch.spm.spatial.preproc.channel.write = [0 0];
matlabbatch.spm.spatial.preproc.tissue(1).tpm = {'/home/niemannf/Documents/MATLAB/spm12/tpm/TPM.nii,1'};
matlabbatch.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch.spm.spatial.preproc.tissue(2).tpm = {'/home/niemannf/Documents/MATLAB/spm12/tpm/TPM.nii,2'};
matlabbatch.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch.spm.spatial.preproc.tissue(3).tpm = {'/home/niemannf/Documents/MATLAB/spm12/tpm/TPM.nii,3'};
matlabbatch.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch.spm.spatial.preproc.tissue(4).tpm = {'/home/niemannf/Documents/MATLAB/spm12/tpm/TPM.nii,4'};
matlabbatch.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch.spm.spatial.preproc.tissue(5).tpm = {'/home/niemannf/Documents/MATLAB/spm12/tpm/TPM.nii,5'};
matlabbatch.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch.spm.spatial.preproc.tissue(6).tpm = {'/home/niemannf/Documents/MATLAB/spm12/tpm/TPM.nii,6'};
matlabbatch.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch.spm.spatial.preproc.warp.mrf = 1;
matlabbatch.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch.spm.spatial.preproc.warp.samp = 3;
matlabbatch.spm.spatial.preproc.warp.write = [1 1];
matlabbatch.spm.spatial.preproc.warp.vox = NaN;
matlabbatch.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];

end

function [matlabbatch] = create_batch_realign_ets_write(number_images)
%% % matlabbatch{1}
a=number_images(:);
matlabbatch.spm.spatial.realign.estwrite.data = {a};
matlabbatch.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch.spm.spatial.realign.estwrite.roptions.which = [2 0];
matlabbatch.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch.spm.spatial.realign.estwrite.roptions.prefix = 'r';
end


function [matlabbatch] = create_batch_realign(number_images)
a=number_images(:);
matlabbatch.spm.spatial.realign.estimate.data = {a};
matlabbatch.spm.spatial.realign.estimate.eoptions.quality = 0.9;
matlabbatch.spm.spatial.realign.estimate.eoptions.sep = 4;
matlabbatch.spm.spatial.realign.estimate.eoptions.fwhm = 5;
matlabbatch.spm.spatial.realign.estimate.eoptions.rtm = 0;
matlabbatch.spm.spatial.realign.estimate.eoptions.interp = 2;
matlabbatch.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
matlabbatch.spm.spatial.realign.estimate.eoptions.weight = '';

end

function [matlabbatch] = create_batch_temp(number_images)

%% there is some issue here only 1 image is loaded
% error message 13-Dec-2023 11:40:47 - Done    'Realign: Estimate & Reslice'
% Item 'Image to Align', field 'val': Number of matching files larger than max allowed, keeping 1/1813 files.
% 13-Dec-2023 11:40:47 - Running 'Slice Timing'

%matlabbatch.spm.temporal.st.scans{1}(1) = cfg_dep('Realign: Estimate: Realigned Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
%matlabbatch.spm.temporal.st.scans{1}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
a=number_images(:);
matlabbatch.spm.temporal.st.scans = {a};

% matlabbatch{2}
matlabbatch.spm.temporal.st.nslices = 72;
matlabbatch.spm.temporal.st.tr = 1;
matlabbatch.spm.temporal.st.ta = 0.986111111111111;
matlabbatch.spm.temporal.st.so = [0 0.41 0.8175 0.245 0.655 0.0825 0.49 0.9 0.3275 0.735 0.165 0.5725 0 0.41 0.8175 0.245 0.655 0.0825 0.49 0.9 0.3275 0.735 0.165 0.5725 0 0.41 0.8175 0.245 0.655 0.0825 0.49 0.9 0.3275 0.735 0.165 0.5725 0 0.41 0.8175 0.245 0.655 0.0825 0.49 0.9 0.3275 0.735 0.165 0.5725 0 0.41 0.8175 0.245 0.655 0.0825 0.49 0.9 0.3275 0.735 0.165 0.5725 0 0.41 0.8175 0.245 0.655 0.0825 0.49 0.9 0.3275 0.735 0.165 0.5725];
matlabbatch.spm.temporal.st.refslice = 1;
matlabbatch.spm.temporal.st.prefix = 'a';
end


function [matlabbatch] = create_batch_normalise(path_y_T1_image,number_images)
% matlabbatch{5}

%matlabbatch.spm.spatial.normalise.estwrite.subj.vol(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
%matlabbatch.spm.spatial.normalise.estwrite.subj.resample(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));

a=number_images(:);
matlabbatch.spm.spatial.normalise.estwrite.subj.vol = {path_y_T1_image};
matlabbatch.spm.spatial.normalise.estwrite.subj.resample = a;

matlabbatch.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch.spm.spatial.normalise.estwrite.eoptions.tpm = {'/home/niemannf/Documents/MATLAB/spm12/tpm/TPM.nii'};
matlabbatch.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
matlabbatch.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
end

function [matlabbatch] = create_batch_smooth(number_images)
% matlabbatch{6}
a=number_images(:); 
%matlabbatch.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch.spm.spatial.smooth.data=a;
matlabbatch.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch.spm.spatial.smooth.dtype = 0;
matlabbatch.spm.spatial.smooth.im = 0;
matlabbatch.spm.spatial.smooth.prefix = 's';
end

function [number_images]=get_images(ref_image_path, ref_image_folder,ref_image_name,pre)
V=spm_vol(ref_image_path);
for i=1:length(V)
number_images{i}=fullfile(ref_image_folder,[pre,ref_image_name,',',num2str(i)]);
end
end





   