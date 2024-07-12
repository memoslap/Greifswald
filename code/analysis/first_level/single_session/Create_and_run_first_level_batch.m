%path = "/media/Data03/Studies/MeMoSLAP/SUBJECTS";

first_level=pwd;
first_level=fileparts(first_level);
first_level=fileparts(first_level);
first_level=fileparts(first_level);
first_level=fileparts(first_level);
path_gr=first_level;

path = fullfile(path_gr,'derivatives/preprocessing/motion_corrected');

path_anaylsis=fullfile(path_gr,'derivatives/analysis/first_level/onset_corrected');
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
    
        %%% loop throught Session for each Subject
        % get all files in every session
        for subfolder_session = 1:length(SesFolderNames)
    
            ses_name = erase(SesFolderNames{subfolder_session},"-");
            path=fullfile(SesFolderDir{subfolder_session}, SesFolderNames{subfolder_session});
            
            filelist.(sub_name).Session.(ses_name).REGRESSOR_realign = dir(fullfile(path, '**/rp*.txt'));
            %filelist.(sub_name).Session.(ses_name).REGRESSOR_outlier = dir(fullfile(path, '**/*Regressor_Outlier.csv'));
            
            %filelist.(sub_name).Session.(ses_name).BOLD = dir(fullfile(path, '**/swa*task*run-*_bold.nii'));
            filelist.(sub_name).Session.(ses_name).BOLD = dir(fullfile(path, '**/swa*task*_bold.nii'));

            filelist.(sub_name).Session.(ses_name).ONSETS = dir(fullfile(path, '**/Onsets_Durations_Names.mat'));
            
            filelist.(sub_name).Session.(ses_name).LOGFILE = dir(fullfile(path, '**/*task-learning.log'));
        end       
end

%% loop through all subjects and create matlabbatch
All_Subjects=fieldnames(filelist);
for sub =1:length(All_Subjects)

%% loop through all sessions

    all_session = {'1','2','3','4'};
    %all_session = {'3','4'};
    for session_num = 1:length(all_session)
        if (isfield(filelist.(All_Subjects{sub}).Session,(['ses',all_session{session_num}]))) 
           %session_num=3;
           %sub=10;
            if (~isempty(filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD)) && ...
                (~isempty(filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).REGRESSOR_realign)) && ...
                (~isempty(filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).ONSETS))
                %try
                Regressor_real_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).REGRESSOR_realign.folder;
                Regressor_real_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).REGRESSOR_realign.name;
                Regressor_real_path = fullfile(Regressor_real_folder, Regressor_real_name);

                %Regressor_out_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).REGRESSOR_outlier.folder;
                %Regressor_out_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).REGRESSOR_outlier.name;
                %Regressor_out_path = fullfile(Regressor_out_folder, Regressor_out_name);

                BOLD_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD.folder;
                BOLD_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).BOLD.name;
                BOLD_path = fullfile(BOLD_folder, BOLD_name);

                Onsets_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).ONSETS.folder;
                Onsets_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).ONSETS.name;
                Onsets_path = fullfile(Onsets_folder, Onsets_name);
                
                Logfile_folder = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).LOGFILE.folder;
                Logfile_name = filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).LOGFILE.name;
                Logifle_path = fullfile(Logfile_folder, Logfile_name);
                
                %% names in onsets are cells, change to char
                load(Onsets_path);
                for i=1:length(names)
                    if ischar(names{i})
                        continue
                    else
                        names{i}=names{i}{:};
                    end
                end
                Onsets_path=fullfile(Onsets_folder,'Onsets_final.mat');
                save(Onsets_path, 'names', 'onsets', 'durations')

                %% load log file to get experiment name
                Logfile_name= filelist.(All_Subjects{sub}).Session.(['ses',all_session{session_num}]).LOGFILE.name;
                if contains(Logfile_name, 'p1')
                    Experiment = 'P1';
                elseif contains(Logfile_name,'p3')
                    Experiment = 'P3';
                else
                    Experiment = 'P7';
                end
                    %if we want to run from the start? just run the second line and comment the first line 
                    %if we want to skip the already done subjects, then run the first line and comment the second one
                    if (~exist(fullfile(path_anaylsis,SubFolderNames{sub},['ses-',all_session{session_num}],'SPM.mat'),'file')) && ...
                        (exist(BOLD_path,'file'))
                    %if exist(BOLD_path,'file')

                        
                        Output_dir=fullfile(path_anaylsis, SubFolderNames{sub},['ses-',all_session{session_num}]);    
                        
                        if ~exist(Output_dir,'dir')
                            mkdir(Output_dir)
                        end
                        if ~exist(fullfile(Output_dir,'code'), 'dir')
                            mkdir(fullfile(Output_dir,'code'))
                        end


                        %if ~exist(Regressor_out_path,'file')  % use only
                        %for motion correction
                        
                        number_images=get_images(BOLD_path, BOLD_folder,BOLD_name);
                        matlabbatch{1} = create_first_level_batch(Output_dir,number_images, Onsets_path, Regressor_real_path);
                        matlabbatch{2} = create_estimate();
                        matlabbatch{3} = create_contrasts(Experiment);
                        
                        %else
                        %% concatenate both regressor use only when outlier is need
                        real_reg=importdata(Regressor_real_path);

%                         Regressor_path=fullfile(Regressor_real_folder,'final_regressor.txt');
%                         writematrix(real_reg,Regressor_path);
%                         %save final_regressor.txt real_reg -ascii -double
%                         figure('Name',strcat(All_Subjects{sub},'_',['ses',all_session{session_num}]));
%                         plot(real_reg);
%                         f=gcf;
%                         path_image=fullfile(path_anaylsis,SubFolderNames{sub},['ses-',all_session{session_num}],'realignment_outlier.png');
%                         exportgraphics(f,path_image,'Resolution',300)
%                         
% 
%                         %% run firstlevel
%                         number_images=get_images(BOLD_path, BOLD_folder,BOLD_name);
%                         matlabbatch{1} = create_first_level_batch(Output_dir,number_images, Onsets_path, Regressor_path);
%                         matlabbatch{2} = create_estimate();
%                         matlabbatch{3} = create_contrasts(Experiment);
                        
                        %end
    
  
                        firstlevel_batch = fullfile(path_anaylsis, SubFolderNames{sub},['ses-',all_session{session_num}],'code','firstlevel_batch.mat');
                        save(firstlevel_batch,'matlabbatch');
                        spm_jobman('run',matlabbatch);
    
    
                    end

                else 
                    fprintf('first level anaylsis already done for session %s for subject %s\n or not all needed files', all_session{session_num},All_Subjects{sub})
            end
        
        end
    end
end
        
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


function [matlabbatch] = create_first_level_batch(Output_dir,number_images, Onsets_path, Regressor_path)
dummy_number_images=number_images(:);
matlabbatch.spm.stats.fmri_spec.dir = {Output_dir};
matlabbatch.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch.spm.stats.fmri_spec.timing.RT = 1;
matlabbatch.spm.stats.fmri_spec.timing.fmri_t = 12;
matlabbatch.spm.stats.fmri_spec.timing.fmri_t0 = 7;
matlabbatch.spm.stats.fmri_spec.sess.scans = dummy_number_images;
matlabbatch.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch.spm.stats.fmri_spec.sess.multi = {Onsets_path};
matlabbatch.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch.spm.stats.fmri_spec.sess.multi_reg = {Regressor_path};
matlabbatch.spm.stats.fmri_spec.sess.hpf = 128;
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
matlabbatch.spm.stats.con.consess{1}.tcon.name = 'Learning > Control';
matlabbatch.spm.stats.con.consess{1}.tcon.weights = [1 1 1 1 -1 -1 -1 -1];
matlabbatch.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{2}.tcon.name = 'Learning < Control';
matlabbatch.spm.stats.con.consess{2}.tcon.weights = [-1 -1 -1 -1 1 1 1 1 ];
matlabbatch.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{3}.tcon.name = 'Learning st1/2 > Control';
matlabbatch.spm.stats.con.consess{3}.tcon.weights = [1 1 0 0 1 1 1 1];
matlabbatch.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{4}.tcon.name = 'Learning st1/2 < Control';
matlabbatch.spm.stats.con.consess{4}.tcon.weights = [-1 -1 0 0 1 1 1 1];
matlabbatch.spm.stats.con.consess{4}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{5}.tcon.name = 'Learning st1/2 > Learning st3/4';
matlabbatch.spm.stats.con.consess{5}.tcon.weights = [1 1 -1 -1];
matlabbatch.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{6}.tcon.name = 'Learning st1/2 < Learning st3/4';
matlabbatch.spm.stats.con.consess{6}.tcon.weights = [-1 -1 1 1];
matlabbatch.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{7}.tcon.name = 'Learning st1/2';
matlabbatch.spm.stats.con.consess{7}.tcon.weights = [1 1 ];
matlabbatch.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{8}.tcon.name = 'Learning st3/4';
matlabbatch.spm.stats.con.consess{8}.tcon.weights = [0 0 1 1];
matlabbatch.spm.stats.con.consess{8}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{9}.tcon.name = 'Learning st1/2 > control';
matlabbatch.spm.stats.con.consess{9}.tcon.weights = [1 1 0 0 -1 -1 -1 -1];
matlabbatch.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{10}.tcon.name = 'Learning st1/2 < control';
matlabbatch.spm.stats.con.consess{10}.tcon.weights = [-1 -1 0 0 1 1 1 1];
matlabbatch.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{11}.tcon.name = 'Learning st3/4 > control';
matlabbatch.spm.stats.con.consess{11}.tcon.weights = [0 0 1 1 -1 -1 -1 -1];
matlabbatch.spm.stats.con.consess{11}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{12}.tcon.name = 'Learning st3/4 < control';
matlabbatch.spm.stats.con.consess{12}.tcon.weights = [0 0 -1 -1 1 1 1 1];
matlabbatch.spm.stats.con.consess{12}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{13}.tcon.name = 'Learning';
matlabbatch.spm.stats.con.consess{13}.tcon.weights = [1 1 1 1];
matlabbatch.spm.stats.con.consess{13}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{14}.tcon.name = 'contorl';
matlabbatch.spm.stats.con.consess{14}.tcon.weights = [0 0 0 0 1 1 1 1];
matlabbatch.spm.stats.con.consess{14}.tcon.sessrep = 'none';

%matlabbatch.spm.stats.con.consess{7}.tcon.name = 'motion';
%matlabbatch.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 0 0 1];
%matlabbatch.spm.stats.con.consess{7}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.delete = 0;

    elseif contains(Experiment,'P7')

matlabbatch.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
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
matlabbatch.spm.stats.con.consess{6}.tcon.name = 'both-sesII,IC > CI,CC';
matlabbatch.spm.stats.con.consess{6}.tcon.weights = [1 -1 1 -1 0 0 0 0 0 0 1 -1 1 -1];
matlabbatch.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{7}.tcon.name = 'both-ses II';
matlabbatch.spm.stats.con.consess{7}.tcon.weights = [1 -1 -1 -1 0 0 0 0 0 0 1 -1 -1 -1];
matlabbatch.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{8}.tcon.name = 'both-ses CI';
matlabbatch.spm.stats.con.consess{8}.tcon.weights = [-1 1 -1 -1 0 0 0 0 0 0 -1 1 -1 -1];
matlabbatch.spm.stats.con.consess{8}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{9}.tcon.name = 'both-ses IC';
matlabbatch.spm.stats.con.consess{9}.tcon.weights = [-1 -1 1 -1 0 0 0 0 0 0 -1 -1 1 -1];
matlabbatch.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{10}.tcon.name = 'both-ses CC';
matlabbatch.spm.stats.con.consess{10}.tcon.weights = [-1 -1 -1 1 0 0 0 0 0 0 -1 -1 -1 1];
matlabbatch.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{11}.tcon.name = 'both-ses CI > II';
matlabbatch.spm.stats.con.consess{11}.tcon.weights = [-1 1 0 0 0 0 0 0 0 0 -1 1 0 0];
matlabbatch.spm.stats.con.consess{11}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{12}.tcon.name = 'both-ses CI,CC > II,IC';
matlabbatch.spm.stats.con.consess{12}.tcon.weights = [-1 1 -1 1 0 0 0 0 0 0 -1 1 -1 1];
matlabbatch.spm.stats.con.consess{12}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{13}.tcon.name = 'ses3 II > ses4 II';
matlabbatch.spm.stats.con.consess{13}.tcon.weights = [1 0 0 0 0 0 0 0 0 0 -1 0 0 0];
matlabbatch.spm.stats.con.consess{13}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{14}.tcon.name = 'ses3 II < ses4 II';
matlabbatch.spm.stats.con.consess{14}.tcon.weights = [-1 0 0 0 0 0 0 0 0 0 1 0 0 0];
matlabbatch.spm.stats.con.consess{14}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{15}.tcon.name = 'ses3 II,IC > ses4 II,IC';
matlabbatch.spm.stats.con.consess{15}.tcon.weights = [1 0 1 0 0 0 0 0 0 0 -1 0 -1 0];
matlabbatch.spm.stats.con.consess{15}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{16}.tcon.name = 'ses3 II,IC < ses4 II,IC';
matlabbatch.spm.stats.con.consess{16}.tcon.weights = [-1 0 -1 0 0 0 0 0 0 0 1 0 1 0];
matlabbatch.spm.stats.con.consess{16}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{17}.tcon.name = 'ses3 CI > ses4 CI';
matlabbatch.spm.stats.con.consess{17}.tcon.weights = [0 1 0 0 0 0 0 0 0 0 0 -1 0 0];
matlabbatch.spm.stats.con.consess{17}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{18}.tcon.name = 'ses3 CI < ses4 CI';
matlabbatch.spm.stats.con.consess{18}.tcon.weights = [0 -1 0 0 0 0 0 0 0 0 0 1 0 0];
matlabbatch.spm.stats.con.consess{18}.tcon.sessrep = 'none';

matlabbatch.spm.stats.con.consess{17}.tcon.name = 'ses3 CI,CC > ses4 CI,CC';
matlabbatch.spm.stats.con.consess{17}.tcon.weights = [0 1 0 1 0 0 0 0 0 0 0 -1 0 -1];
matlabbatch.spm.stats.con.consess{17}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.consess{18}.tcon.name = 'ses3 CI,CC < ses4 CI,CC';
matlabbatch.spm.stats.con.consess{18}.tcon.weights = [0 -1 0-1 0 0 0 0 0 0 0 1 0 1];
matlabbatch.spm.stats.con.consess{18}.tcon.sessrep = 'none';
matlabbatch.spm.stats.con.delete = 0;  
    end
end