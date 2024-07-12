%path = "/media/Data03/Studies/MeMoSLAP/SUBJECTS";

first_level=pwd;
first_level=fileparts(first_level);
first_level=fileparts(first_level);
first_level=fileparts(first_level);
first_level=fileparts(first_level);
path_gr=first_level;

path = fullfile(path_gr,'derivatives/preprocessing/motion_corrected');

path_anaylsis=fullfile(path_gr,'derivatives/analysis/first_level/task_pmod');
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
    coregistration_path=fullfile(path_gr,'/derivatives/preprocessing/motion_corrected');
        fprintf('Sub folder #%d = %s\n', subfolder_subject, SubFolderNames{subfolder_subject});
        path = fullfile(SubFolderDir{subfolder_subject}, SubFolderNames{subfolder_subject});
        [SesFolderDir, SesFolderNames] = get_subdir(path);
        sub_name = erase(SubFolderNames{subfolder_subject},"-");
        filelist.(sub_name).Session_list = SesFolderNames;
%% loop throught Session for each Subject
% get all files in every session

        for subfolder_session = 1:length(SesFolderNames)
    
            ses_name = erase(SesFolderNames{subfolder_session},"-");
            path=fullfile(SesFolderDir{subfolder_session}, SesFolderNames{subfolder_session});
            
            filelist.(sub_name).Session.(ses_name).REGRESSOR_realign = dir(fullfile(path, '**/rp*.txt'));
            %filelist.(sub_name).Session.(ses_name).REGRESSOR_outlier = dir(fullfile(path, '**/*Regressor_Outlier.csv'));
            
            %filelist.(sub_name).Session.(ses_name).BOLD = dir(fullfile(path, '**/swa*task*run-*_bold.nii'));
            filelist.(sub_name).Session.(ses_name).BOLD = dir(fullfile(path, '**/swa*task*_bold.nii'));

            filelist.(sub_name).Session.(ses_name).ONSETS = dir(fullfile(path_anaylsis,SubFolderNames{subfolder_subject},SesFolderNames{subfolder_session}, '**/Onsets_Durations_Names.mat'));
            filelist.(sub_name).Session.(ses_name).ONSETS_PMOD = dir(fullfile(path_anaylsis,SubFolderNames{subfolder_subject},SesFolderNames{subfolder_session}, '**/Onsets_Durations_Names_pmod.mat'));
            filelist.(sub_name).Session.(ses_name).LOGFILE = dir(fullfile(path, '**/*task-learning.log'));
        end       
end
%% loop through all subjects and create matlabbatch

All_Subjects=fieldnames(filelist);
for sub =1:length(All_Subjects)
%% loop through all sessions

all_session = {'3','4'};
    if isfield(filelist.(All_Subjects{sub}).Session,('ses3')) && ...
        (~isempty(filelist.(All_Subjects{sub}).Session.('ses3').BOLD)) && ...
    (~isempty(filelist.(All_Subjects{sub}).Session.('ses3').REGRESSOR_realign)) && ...
    (~isempty(filelist.(All_Subjects{sub}).Session.(('ses3')).ONSETS)) && ...
    (~isempty(filelist.(All_Subjects{sub}).Session.(('ses3')).ONSETS_PMOD)) && ...
    isfield(filelist.(All_Subjects{sub}).Session,('ses4')) && ...
    (~isempty(filelist.(All_Subjects{sub}).Session.('ses4').BOLD)) && ...
    (~isempty(filelist.(All_Subjects{sub}).Session.('ses4').REGRESSOR_realign)) && ...
    (~isempty(filelist.(All_Subjects{sub}).Session.(('ses4')).ONSETS)) && ...
    (~isempty(filelist.(All_Subjects{sub}).Session.(('ses3')).ONSETS_PMOD))
  
        for session_num = 1:length(all_session)
            if (isfield(filelist.(All_Subjects{sub}).Session,(['ses',all_session{session_num}]))) 
               %session_num=3;
               %sub=10;
%% load log file to get experiment name

                    Logfile_name= filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).LOGFILE.name;
                        if contains(Logfile_name, 'p1')
                            Experiment = 'P1';
                        elseif contains(Logfile_name,'p3')
                            Experiment = 'P3';
                        else
                            Experiment = 'P7';
                        end
%% concatenate session 3 and 4

                    if all_session{session_num}=='3'
                        [Regressor_real_path_ses3,BOLD_path_ses3,Onsets_path_ses3, Onsets_path_mod_ses3 ,Logifle_path_ses3,number_images_ses3,Regressor_real_folder_ses3]=create_path(filelist,All_Subjects, sub,all_session, session_num);
                        [error_var]=plot_motion(Regressor_real_path_ses3,Regressor_real_folder_ses3);
                        if exist('error_var') && error_var==1
                            fprintf('Motion error for Subject %s session 3',SubFolderNames{sub})
                            error_file= fopen(fullfile(Regressor_real_folder_ses3, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!There_is_a_motion_error!!!!!!!!!!!!!!!!!!!!!!!!!!!!.pdf'),'w');
                            fclose(error_file);
                            clear error_file
                            clear error_var
                        end
                    elseif all_session{session_num}=='4'
                        [Regressor_real_path_ses4,BOLD_path_ses4,Onsets_path_ses4, Onsets_path_mod_ses4,Logifle_path_ses4,number_images_ses4,Regressor_real_folder_ses4]=create_path(filelist,All_Subjects, sub,all_session, session_num);
                        [error_var]=plot_motion(Regressor_real_path_ses4,Regressor_real_folder_ses4);
                        if exist('error_var') && error_var==1
                            fprintf('Motion error for Subject %s session 4',SubFolderNames{sub})
                            error_file= fopen(fullfile(Regressor_real_folder_ses4,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!There_is_a_motion_error!!!!!!!!!!!!!!!!!!!!!!!!!!!!.pdf'),'w');
                            fclose(error_file);
                            clear error_file
                            clear error_var
                        end
                     end              
            end           
        end

        %if we want to run from the start? just run the second line and comment the first line 
        %if we want to skip the already done subjects, then run the first line and comment the second one
        if (~exist(fullfile(path_anaylsis,SubFolderNames{sub},'SPM.mat'),'file')) 
        %if exist(BOLD_path_ses3,'file') && exist(BOLD_path_ses4,'file') 
        Output_dir=fullfile(path_anaylsis, SubFolderNames{sub});    
           
            if ~exist(Output_dir,'dir')
                mkdir(Output_dir)
            end
            if ~exist(fullfile(Output_dir,'code'), 'dir')
                mkdir(fullfile(Output_dir,'code'))
            end
        fprintf('SPM.mat should run for subject %s',All_Subjects{sub})
        matlabbatch{1} = create_first_level_batch(Output_dir,number_images_ses3, Onsets_path_ses3, Regressor_real_path_ses3,number_images_ses4, Onsets_path_ses4, Regressor_real_path_ses4);
        matlabbatch{2} = create_estimate();
        matlabbatch{3} = create_contrasts(Experiment);
    
        firstlevel_batch = fullfile(path_anaylsis, SubFolderNames{sub},'code','firstlevel_batch.mat');
        save(firstlevel_batch,'matlabbatch');
        %%uncomment if you want the batches run automatically
        spm_jobman('run',matlabbatch);
        
        else 
        fprintf('SPM.mat already exists %s for subject ',All_Subjects{sub})
    
        end
    else
    fprintf('ses3 and ses4 not complete subject %s\n or not all needed files \n',All_Subjects{sub})    
   end 
end


%%
function [number_images]=get_images(BOLD_path, BOLD_folder,BOLD_name)
V=spm_vol(BOLD_path);
for i=1:length(V)
number_images{i}=fullfile(BOLD_folder,[BOLD_name,',',num2str(i)]);
end
end

function [SubFolderDir_func, SubFolderNames_func] = get_subdir(act_directory)
files_sub = dir(act_directory);
files_sub(~[files_sub.isdir])=[];
tf = ismember({files_sub.name}, {'.', '..'});
files_sub(tf) = [];  %remove current and parent directory.
SubFolderNames_func = {files_sub.name};
SubFolderDir_func = {files_sub.folder};
end


function [matlabbatch] = create_first_level_batch(Output_dir,number_images_ses3, Onsets_path_ses3, Regressor_path_ses3,number_images_ses4, Onsets_path_ses4, Regressor_path_ses4)
dummy_number_images_ses3=number_images_ses3(:);
dummy_number_images_ses4=number_images_ses4(:);
matlabbatch.spm.stats.fmri_spec.dir = {Output_dir};
matlabbatch.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch.spm.stats.fmri_spec.timing.RT = 1;
matlabbatch.spm.stats.fmri_spec.timing.fmri_t = 12;
matlabbatch.spm.stats.fmri_spec.timing.fmri_t0 = 7;
matlabbatch.spm.stats.fmri_spec.sess(1).scans = dummy_number_images_ses3;
matlabbatch.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch.spm.stats.fmri_spec.sess(1).multi = {Onsets_path_ses3};
matlabbatch.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch.spm.stats.fmri_spec.sess(1).multi_reg = {Regressor_path_ses3};
matlabbatch.spm.stats.fmri_spec.sess(1).hpf = 128;
matlabbatch.spm.stats.fmri_spec.sess(2).scans = dummy_number_images_ses4;
matlabbatch.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch.spm.stats.fmri_spec.sess(2).multi = {Onsets_path_ses4};
matlabbatch.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch.spm.stats.fmri_spec.sess(2).multi_reg = {Regressor_path_ses4};
matlabbatch.spm.stats.fmri_spec.sess(2).hpf = 128;
matlabbatch.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch.spm.stats.fmri_spec.volt = 1;
matlabbatch.spm.stats.fmri_spec.global = 'None';
matlabbatch.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch.spm.stats.fmri_spec.mask = {''};
matlabbatch.spm.stats.fmri_spec.cvi = 'AR(1)';
end

function [matlabbatch] = create_estimate()
matlabbatch.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch.spm.stats.fmri_est.write_residuals = 0;
matlabbatch.spm.stats.fmri_est.method.Classical = 1;
end

function [matlabbatch] = create_contrasts(Experiment)
    if (contains(Experiment,'P1')) || (contains(Experiment,'P3'))
%matlabbatch.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch.spm.stats.con.consess{1}.tcon.name = 'L > C';
matlabbatch.spm.stats.con.consess{1}.tcon.weights = [1 0 1 0 1 0 1 0 -1 0 -1 0 -1 0 -1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 -1 0 -1 0 -1 0 -1 0];
matlabbatch.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{2}.tcon.name = 'L < C';
matlabbatch.spm.stats.con.consess{2}.tcon.weights = [-1 0 -1 0 -1 0 -1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 -1 0 -1 0 -1 0 -1 0 1 0 1 0 1 0 1 0];
matlabbatch.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{3}.tcon.name = 'pmod L > C';
matlabbatch.spm.stats.con.consess{3}.tcon.weights = [1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1];
matlabbatch.spm.stats.con.consess{3}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{4}.tcon.name = 'pmod L < C';
matlabbatch.spm.stats.con.consess{4}.tcon.weights = [-1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1];
matlabbatch.spm.stats.con.consess{4}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{3}.tcon.name = 'L_st1/2 > C';
% matlabbatch.spm.stats.con.consess{3}.tcon.weights = [1 1 0 0 -1 -1 -1 -1 0 0 0 0 0 0 1 1 0 0 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{3}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{4}.tcon.name = 'L_st1/2 < C';
% matlabbatch.spm.stats.con.consess{4}.tcon.weights = [-1 -1 0 0 1 1 1 1 0 0 0 0 0 0  -1 -1 0 0 1 1 1 1];
% matlabbatch.spm.stats.con.consess{4}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{5}.tcon.name = 'L_st1/2 > L_st3/4';
% matlabbatch.spm.stats.con.consess{5}.tcon.weights = [1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1];
% matlabbatch.spm.stats.con.consess{5}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{6}.tcon.name = 'L_st1/2 < L_st3/4';
% matlabbatch.spm.stats.con.consess{6}.tcon.weights = [-1 -1 1 1 0 0 0 0 0 0 0 0 0 0 -1 -1 1 1];
% matlabbatch.spm.stats.con.consess{6}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{7}.tcon.name = 'L_st1/2';
% matlabbatch.spm.stats.con.consess{7}.tcon.weights = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1  ];
% matlabbatch.spm.stats.con.consess{7}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{8}.tcon.name = 'L_st3/4';
% matlabbatch.spm.stats.con.consess{8}.tcon.weights = [0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1];
% matlabbatch.spm.stats.con.consess{8}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{9}.tcon.name = 'L_st3/4 > C';
% matlabbatch.spm.stats.con.consess{9}.tcon.weights = [0 0 1 1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 1 1 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{9}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{10}.tcon.name = 'L_st3/4 < C';
% matlabbatch.spm.stats.con.consess{10}.tcon.weights = [0 0 -1 -1 1 1 1 1 0 0 0 0 0 0 0 0 -1 -1 1 1 1 1];
% matlabbatch.spm.stats.con.consess{10}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{11}.tcon.name = 'L';
% matlabbatch.spm.stats.con.consess{11}.tcon.weights = [1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1];
% matlabbatch.spm.stats.con.consess{11}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{12}.tcon.name = 'C';
% matlabbatch.spm.stats.con.consess{12}.tcon.weights = [0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1];
% matlabbatch.spm.stats.con.consess{12}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{13}.tcon.name = 'S3_L > S4_L';
% matlabbatch.spm.stats.con.consess{13}.tcon.weights = [1 1 1 1 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 0 0 0 0];
% matlabbatch.spm.stats.con.consess{13}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{14}.tcon.name = 'S3_L < S4_L';
% matlabbatch.spm.stats.con.consess{14}.tcon.weights = [-1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0];
% matlabbatch.spm.stats.con.consess{14}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{15}.tcon.name = 'S3_C > S4_C';
% matlabbatch.spm.stats.con.consess{15}.tcon.weights = [0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{15}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{16}.tcon.name = 'S3_C < S4_C';
% matlabbatch.spm.stats.con.consess{16}.tcon.weights = [0 0 0 0 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 1 1 1 1];
% matlabbatch.spm.stats.con.consess{16}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{17}.tcon.name = 'S3_L_st1/2 > S4_L_st1/2';
% matlabbatch.spm.stats.con.consess{17}.tcon.weights = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0];
% matlabbatch.spm.stats.con.consess{17}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{18}.tcon.name = 'S3_L_st1/2 < S4_L_st1/2';
% matlabbatch.spm.stats.con.consess{18}.tcon.weights = [-1 -1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0];
% matlabbatch.spm.stats.con.consess{18}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{19}.tcon.name = 'S3_L_st3/4 > S4_L_st3/4';
% matlabbatch.spm.stats.con.consess{19}.tcon.weights = [0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 0 0 0 0];
% matlabbatch.spm.stats.con.consess{19}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{20}.tcon.name = 'S3_L_st3/4 < S4_L_st3/4';
% matlabbatch.spm.stats.con.consess{20}.tcon.weights = [0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0];
% matlabbatch.spm.stats.con.consess{20}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{21}.tcon.name = 'S3_only_L > S4_only_L';
% matlabbatch.spm.stats.con.consess{21}.tcon.weights = [1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0];
% matlabbatch.spm.stats.con.consess{21}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{22}.tcon.name = 'S3_only_C < S4_only_C';
% matlabbatch.spm.stats.con.consess{22}.tcon.weights = [0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1];
% matlabbatch.spm.stats.con.consess{22}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{23}.tcon.name = 'S3_L > S3_C & S4_L > S4_C';
% matlabbatch.spm.stats.con.consess{23}.tcon.weights = [1 1 1 1 -1 -1 -1 -1 0 0 0 0 0 0 1 1 1 1 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{23}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{24}.tcon.name = 'S3_L < S3_C & S4_L < S4_C';
% matlabbatch.spm.stats.con.consess{24}.tcon.weights = [-1 -1 -1 -1 1 1 1 1 0 0 0 0 0 0 -1 -1 -1 -1 1 1 1 1];
% matlabbatch.spm.stats.con.consess{24}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{25}.tcon.name = 'S3_L_st1/2 > S3_C & S4_L_st1/2 > S4_C';
% matlabbatch.spm.stats.con.consess{25}.tcon.weights = [1 1 0 0 -1 -1 -1 -1 0 0 0 0 0 0 1 1 0 0 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{25}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{26}.tcon.name = 'S3_L_st1/2 < S3_C & S4_L_st1/2 < S4_C';
% matlabbatch.spm.stats.con.consess{26}.tcon.weights = [-1 -1 0 0 1 1 1 1 0 0 0 0 0 0 -1 -1 0 0 1 1 1 1];
% matlabbatch.spm.stats.con.consess{26}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{27}.tcon.name = 'S3_L_st3/4 > S3_C & S4_L_st3/4 > S4_C';
% matlabbatch.spm.stats.con.consess{27}.tcon.weights = [0 0 1 1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 1 1 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{27}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{28}.tcon.name = 'S3_L_st3/4 < S3_C & S4_L_st3/4 < S4_C';
% matlabbatch.spm.stats.con.consess{28}.tcon.weights = [0 0 -1 -1 1 1 1 1 0 0 0 0 0 0 0 0 -1 -1 1 1 1 1];
% matlabbatch.spm.stats.con.consess{28}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{29}.tcon.name = 'S3_L_st1/2 > S3_L_st3/4 & S4_L_st1/2 > S4_L_st3/4';
% matlabbatch.spm.stats.con.consess{29}.tcon.weights = [1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 0 0 0 0];
% matlabbatch.spm.stats.con.consess{29}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{30}.tcon.name = 'S3_L_st1/2 < S3_L_st3/4 & S4_L_st1/2 < S4_L_st3/4';
% matlabbatch.spm.stats.con.consess{30}.tcon.weights = [-1 -1 1 1 0 0 0 0 0 0 0 0 0 0 -1 -1 1 1 0 0 0 0];
% matlabbatch.spm.stats.con.consess{30}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{31}.tcon.name = 'S3_L > S3_C';
% matlabbatch.spm.stats.con.consess{31}.tcon.weights = [1 1 1 1 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{31}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{32}.tcon.name = 'S3_L < S3_C';
% matlabbatch.spm.stats.con.consess{32}.tcon.weights = [-1 -1 -1 -1 1 1 1 1];
% matlabbatch.spm.stats.con.consess{32}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{33}.tcon.name = 'S4_L > S4_C';
% matlabbatch.spm.stats.con.consess{33}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{33}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{34}.tcon.name = 'S4_L < S4_C';
% matlabbatch.spm.stats.con.consess{34}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 1 1 1 1 ];
% matlabbatch.spm.stats.con.consess{34}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{35}.tcon.name = 'L > C2,3,4';
% matlabbatch.spm.stats.con.consess{35}.tcon.weights = [1 1 1 1 0 -1 -1 -1 0 0 0 0 0 0 1 1 1 1 0 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{35}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{36}.tcon.name = 'L > C3,4';
% matlabbatch.spm.stats.con.consess{36}.tcon.weights = [1 1 1 1 0 0 -1 -1 0 0 0 0 0 0 1 1 1 1 0 0 -1 -1];
% matlabbatch.spm.stats.con.consess{36}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{37}.tcon.name = 'L3,4 > C3,4';
% matlabbatch.spm.stats.con.consess{37}.tcon.weights = [0 0 1 1 0 0 -1 -1 0 0 0 0 0 0 0 0 1 1 0 0 -1 -1];
% matlabbatch.spm.stats.con.consess{37}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{38}.fcon.name = 'ses3 comp ses4 contrast L > C';
% matlabbatch.spm.stats.con.consess{38}.fcon.weights = [1 1 1 1 -1 -1 -1 -1 zeros(1,20)
%                                                       zeros(1,14) 1 1 1 1 -1 -1 -1 -1 zeros(1,6)];                                                 
% matlabbatch.spm.stats.con.consess{38}.fcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{39}.tcon.name = 's3_L > C2,3,4';
% matlabbatch.spm.stats.con.consess{39}.tcon.weights = [1 1 1 1 0 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{39}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{40}.tcon.name = 's4_L > C2,3,4';
% matlabbatch.spm.stats.con.consess{40}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{40}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{31}.tcon.name = 'S3_L st1/2 > S3_C';
% matlabbatch.spm.stats.con.consess{31}.tcon.weights = [1 1 0 0 -1 -1 -1 -1];
% matlabbatch.spm.stats.con.consess{31}.tcon.sessrep = 'none';
% 
% matlabbatch.spm.stats.con.consess{31}.tcon.name = 'S3_L st1/2 < S3_C';
% matlabbatch.spm.stats.con.consess{31}.tcon.weights = [-1 -1 0 0 1 1 1 1];
% matlabbatch.spm.stats.con.consess{31}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.delete = 0;

    elseif contains(Experiment,'P7')
matlabbatch.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

%matlabbatch.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%matlabbatch.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

matlabbatch.spm.stats.con.consess{1}.tcon.name = 'ses3 II > CI';
matlabbatch.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{2}.tcon.name = 'ses4 II > CI';
matlabbatch.spm.stats.con.consess{2}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1 -1 0 0];
matlabbatch.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{3}.tcon.name = 'ses3 II,IC > CI,CC';
matlabbatch.spm.stats.con.consess{3}.tcon.weights = [1 -1 1 -1];
matlabbatch.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{4}.tcon.name = 'ses4 II,IC > CI,CC';
matlabbatch.spm.stats.con.consess{4}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1 -1 1 -1];
matlabbatch.spm.stats.con.consess{4}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{5}.tcon.name = 'both-ses II > CI';
matlabbatch.spm.stats.con.consess{5}.tcon.weights = [1 -1 0 0 0 0 0 0 0 0 1 -1 0 0];
matlabbatch.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{6}.tcon.name = 'both-ses II,IC > CI,CC';
matlabbatch.spm.stats.con.consess{6}.tcon.weights = [1 -1 1 -1 0 0 0 0 0 0 1 -1 1 -1];
matlabbatch.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{7}.tcon.name = 'both-ses II > all';
matlabbatch.spm.stats.con.consess{7}.tcon.weights = [1 -1 -1 -1 0 0 0 0 0 0 1 -1 -1 -1];
matlabbatch.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{8}.tcon.name = 'both-ses CI > all';
matlabbatch.spm.stats.con.consess{8}.tcon.weights = [-1 1 -1 -1 0 0 0 0 0 0 -1 1 -1 -1];
matlabbatch.spm.stats.con.consess{8}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{9}.tcon.name = 'both-ses IC > all';
matlabbatch.spm.stats.con.consess{9}.tcon.weights = [-1 -1 1 -1 0 0 0 0 0 0 -1 -1 1 -1];
matlabbatch.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{10}.tcon.name = 'both-ses CC > all';
matlabbatch.spm.stats.con.consess{10}.tcon.weights = [-1 -1 -1 1 0 0 0 0 0 0 -1 -1 -1 1];
matlabbatch.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{11}.tcon.name = 'ses3 II < CI';
matlabbatch.spm.stats.con.consess{11}.tcon.weights = [-1 1];
matlabbatch.spm.stats.con.consess{11}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{12}.tcon.name = 'ses4 II < CI';
matlabbatch.spm.stats.con.consess{12}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 -1 1 0 0];
matlabbatch.spm.stats.con.consess{12}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{13}.tcon.name = 'ses3 II,IC < CI,CC';
matlabbatch.spm.stats.con.consess{13}.tcon.weights = [-1 1 -1 1];
matlabbatch.spm.stats.con.consess{13}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{14}.tcon.name = 'ses4 II,IC < CI,CC';
matlabbatch.spm.stats.con.consess{14}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 -1 1 -1 1];
matlabbatch.spm.stats.con.consess{14}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{15}.tcon.name = 'both-ses II < CI';
matlabbatch.spm.stats.con.consess{15}.tcon.weights = [-1 1 0 0 0 0 0 0 0 0 -1 1 0 0];
matlabbatch.spm.stats.con.consess{15}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{16}.tcon.name = 'both-ses II,IC < CI,CC';
matlabbatch.spm.stats.con.consess{16}.tcon.weights = [-1 1 -1 1 0 0 0 0 0 0 -1 1 -1 1];
matlabbatch.spm.stats.con.consess{16}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{17}.tcon.name = 'ses3 II > ses4 II';
matlabbatch.spm.stats.con.consess{17}.tcon.weights = [1 0 0 0 0 0 0 0 0 0 -1 0 0 0];
matlabbatch.spm.stats.con.consess{17}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{18}.tcon.name = 'ses3 II < ses4 II';
matlabbatch.spm.stats.con.consess{18}.tcon.weights = [-1 0 0 0 0 0 0 0 0 0 1 0 0 0];
matlabbatch.spm.stats.con.consess{18}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{19}.tcon.name = 'ses3 II,IC > ses4 II,IC';
matlabbatch.spm.stats.con.consess{19}.tcon.weights = [1 0 1 0 0 0 0 0 0 0 -1 0 -1 0];
matlabbatch.spm.stats.con.consess{19}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{20}.tcon.name = 'ses3 II,IC < ses4 II,IC';
matlabbatch.spm.stats.con.consess{20}.tcon.weights = [-1 0 -1 0 0 0 0 0 0 0 1 0 1 0];
matlabbatch.spm.stats.con.consess{20}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{21}.tcon.name = 'ses3 CI > ses4 CI';
matlabbatch.spm.stats.con.consess{21}.tcon.weights = [0 1 0 0 0 0 0 0 0 0 0 -1 0 0];
matlabbatch.spm.stats.con.consess{21}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{22}.tcon.name = 'ses3 CI < ses4 CI';
matlabbatch.spm.stats.con.consess{22}.tcon.weights = [0 -1 0 0 0 0 0 0 0 0 0 1 0 0];
matlabbatch.spm.stats.con.consess{22}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{23}.tcon.name = 'ses3 CI,CC > ses4 CI,CC';
matlabbatch.spm.stats.con.consess{23}.tcon.weights = [0 1 0 1 0 0 0 0 0 0 0 -1 0 -1];
matlabbatch.spm.stats.con.consess{23}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{24}.tcon.name = 'ses3 CI,CC < ses4 CI,CC';
matlabbatch.spm.stats.con.consess{24}.tcon.weights = [0 -1 0 -1 0 0 0 0 0 0 0 1 0 1];
matlabbatch.spm.stats.con.consess{24}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{25}.tcon.name = 'ses3 II,CI > IC,CC (Stroop)';
matlabbatch.spm.stats.con.consess{25}.tcon.weights = [1 1 -1 -1];
matlabbatch.spm.stats.con.consess{25}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{26}.tcon.name = 'ses4 II,CI > IC,CC (Stroop)';
matlabbatch.spm.stats.con.consess{26}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1 1 -1 -1];
matlabbatch.spm.stats.con.consess{26}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{27}.tcon.name = 'both-ses II,CI > IC,CC (Stroop)';
matlabbatch.spm.stats.con.consess{27}.tcon.weights = [1 1 -1 -1 0 0 0 0 0 0 1 1 -1 -1];
matlabbatch.spm.stats.con.consess{27}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{28}.tcon.name = 'ses3 II,CI < IC,CC (Revstroop)';
matlabbatch.spm.stats.con.consess{28}.tcon.weights = [-1 -1 1 1];
matlabbatch.spm.stats.con.consess{28}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{29}.tcon.name = 'ses4 II,CI < IC,CC (Revstroop)';
matlabbatch.spm.stats.con.consess{29}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 -1 -1 1 1];
matlabbatch.spm.stats.con.consess{29}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{30}.tcon.name = 'both-ses II,CI < IC,CC (Revstroop)';
matlabbatch.spm.stats.con.consess{30}.tcon.weights = [-1 -1 1 1 0 0 0 0 0 0 -1 -1 1 1];
matlabbatch.spm.stats.con.consess{30}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{31}.fcon.name = 'ses3 comp ses4 contrast CI,CC';
matlabbatch.spm.stats.con.consess{31}.fcon.weights = [0 1 0 1 0 0 0 0 0 0 0 0 0 0
                                                      0 0 0 0 0 0 0 0 0 0 0 1 0 1];                                                 
matlabbatch.spm.stats.con.consess{31}.fcon.sessrep = 'none';

matlabbatch.spm.stats.con.delete = 0;  
    end
end

function [Regressor_real_path,BOLD_path,Onsets_path,Onsets_path_mod,Logifle_path,number_images,Regressor_real_folder] = create_path(filelist,All_Subjects, sub,all_session, session_num)
Regressor_real_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).REGRESSOR_realign.folder;
Regressor_real_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).REGRESSOR_realign.name;
Regressor_real_path = fullfile(Regressor_real_folder, Regressor_real_name);

BOLD_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD.folder;
BOLD_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD.name;
BOLD_path = fullfile(BOLD_folder, BOLD_name);

Onsets_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).ONSETS.folder;
Onsets_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).ONSETS.name;
Onsets_path = fullfile(Onsets_folder, Onsets_name);

Onsets_folder_mod = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).ONSETS_PMOD.folder;
Onsets_name_mod = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).ONSETS_PMOD.name;
Onsets_path_mod = fullfile(Onsets_folder_mod, Onsets_name_mod);

Logfile_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).LOGFILE.folder;
Logfile_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).LOGFILE.name;
Logifle_path = fullfile(Logfile_folder, Logfile_name);
%% get number of images

number_images=get_images(BOLD_path, BOLD_folder,BOLD_name);
%% names in onsets are cells, change to char

load(Onsets_path);
        for i=1:length(names)
            if ischar(names{i})
                continue
            else
                names{i}=names{i}{:};
            end
        end

load(Onsets_path_mod);
    %for i=1:length(pmod_wrong.name)
    %        if ischar(pmod_wrong.name{i})
    %            continue
    %        else
    %            pmod_wrong.name{i}=pmod_wrong.name{i}{:};
    %        end
    %end

    for i=1:8
        %i=int32(i);
        %pmod{i} = struct('name',{''},'param',{},'poly',{});
        name = pmod_wrong.name(i);
        poly=cell(1,1);
        poly{1} = pmod_wrong.poly(i);
        param = pmod_wrong.param(i);
    

        pmod(i) = struct('name',name, 'param',{param}, 'poly',{poly});
    end

Onsets_path=fullfile(Onsets_folder,'Onsets_final.mat');
save(Onsets_path, 'names', 'onsets', 'durations','pmod')
end

function [error_var]=plot_motion(Regressor_file, Regressor_folder)
    data=readmatrix(Regressor_file);
    gcf=figure;
    subplot(2,1,1)
    plot(data(:,1:3))
    
    lgd=legend('X','Y','Z');
    lgd.Location = 'eastoutside';
    subplot(2,1,2)
    plot(data(:,4:6))
    
    lgd2=legend('trans1','trans2','trans3'); 
    lgd2.Location = 'eastoutside';
    exportgraphics(gcf, fullfile(Regressor_folder,'Plot_motion.pdf'));
    close(gcf)

    if ~isempty(find(abs(diff(data(:,1:3))) > 1));
        error_var=1;
    else
        error_var=0;
    end    
  

end