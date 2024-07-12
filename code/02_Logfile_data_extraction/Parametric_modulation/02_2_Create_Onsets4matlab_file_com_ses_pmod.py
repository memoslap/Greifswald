# use conda environmen MRI_Analysis (on user niemannf, for neuro-forschung?)

# %%
from openpyxl import load_workbook
import pandas as pd
import numpy as np
import scipy.io
import csv
import os
#import pickle
import nibabel as nib
import glob
from importlib import reload
#import Experiment as Exp
import re
import Experiment_2_pmod as Exp
reload(Exp)

#import Experiment_motion as Exp_mot
#reload(Exp_mot)


# %%
##################### change here ##################################

pmod_col_names='Response_num'

################### end  #####################################
# %%
#path_files=os.getcwd()
#äpath='/media/MeMoSLAP_Subjects/derivatives/preprocessing/motion_corrected' #MeMoSLAP Linux

path_BIDS_source='/media/Data03/Studies/MeMoSLAP'

path='/media/Data03/Studies/MeMoSLAP/derivatives/preprocessing/SPM-standard' # Meinzer_Linux



path_out = f'/media/Data03/Studies/MeMoSLAP/derivatives/analysis/first_level/task_pmod/{pmod_col_names}'



list_sub=[a for a in os.listdir(path) if os.path.isdir(os.path.join(path, a))]

list_ses=[]
for sub in list_sub:
    list_ses.append([a for a in os.listdir(os.path.join(path,sub)) if os.path.isdir(os.path.join(path,sub, a))])

# bids for behavioral recordings
#https://bids-specification.readthedocs.io/en/latest/modality-specific-files/behavioral-experiments.html

# %% only neede for motion correction

def search_files(directory, search_string):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if search_string in open(file_path).read():
                file_list.append(file_path)
    return file_list

def get_index_of_bold(sub, ses):
    #path_func = f'/media/MeMoSLAP_Subjects/derivatives/preprocessing/motion_corrected/{sub}/{ses}/func'  # MeMoSLAP Linux
    path_func =  os.path.join(path,sub,ses,'func')    #MeinzerMRI_Linux
    
    bold=glob.glob(f'{path_func}/{sub}*{ses}*run-01_bold.nii')
    
    if len(bold) !=1:
        print(f'no or to much bold images in directory {sub} {ses}')
        
    elif os.path.isfile(bold[0]): 
        bold_img=nib.load(bold[0])
        Vol=bold_img.shape[3]
        return Vol, bold[0]



# read files for parametric modeling
    
    
   

# %% only neede for motion correction
#path_FD = '/media/MeMoSLAP_Subjects/malinowskir/OUT/work/mriqc_wf/funcMRIQC/fMRI_HMC'
#threshold_FD = 1.5
#dict_outlier={}

# %%
def read_in_df(path,sub, ses, project,*args):
    if os.path.isfile(os.path.join(path,sub,ses,'beh',f'{sub}_{ses}_task-{project}_acq-1_desc-data.tsv')) or os.path.isfile(os.path.join(path,sub,ses,'beh',f'{sub}_{ses}_task-{project}_acq-2_desc-data.tsv')):

        path_files=os.path.join(path, sub, ses,'beh')
        
        if os.path.isfile(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-1_desc-pmod.tsv')) and os.path.isfile(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-1_desc-pulses.tsv')):
            df = pd.read_csv(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-1_desc-data.tsv'), delimiter="\t")
            df_pmod = pd.read_csv(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-1_desc-pmod.tsv'), delimiter="\t")
            df_pulse = pd.read_csv(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-1_desc-pulses.tsv'), delimiter="\t")
            acq='acq-1'
        elif os.path.isfile(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-2_desc-pmod.tsv')) and os.path.isfile(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-2_desc-pulses.tsv')):
            df = pd.read_csv(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-2_desc-data.tsv'), delimiter="\t")
            df_pmod = pd.read_csv(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-2_desc-pmod.tsv'), delimiter="\t")
            df_pulse = pd.read_csv(os.path.join(path_files,f'{sub}_{ses}_task-{project}_acq-2_desc-pulses.tsv'), delimiter="\t")
            acq='acq-2'
        else:
            df = None
            
        # presentation gives time in 10000 units so devide every time containing columne by 10000
        # in logfile prsentation is starting and waiting for trigger from scanner
        # first Pulse 30 gives time for first image so this should be point zero for all other onsets
##################### change here ##################################

        try: 
        #elif len(df)>0 and len(df_pulse) > 0:
                        
            for i in ["stimulus_onset", "response_onset"]:
            
                df[i]=df[i]-df_pulse['Time'][0] 
                df[i] = df[i]/10000 # add additional 2 seconds for the skipped image at the beginning which showed wrong timing 
                print('first pulse',df_pulse['Time'][0],'diff pulse',df[i][0] )

    #get's the Volume of images

                if (ses == 'ses-4') and (Vol):
                    #df[i]=df[i]+args[0]+1 # Add number of Volumes to ses-4 Onsets
                    df[i]=df[i] # Add number of Volumes to ses-4 Onsets
                elif (ses == 'ses-4') and (not Vol):
                    print(f'Vol problem for subject {sub} session 4')
    ################### end #####################################
            df=df.fillna(0)
            return df, acq, df_pmod
            #print(df.stimulus_onset[0])
            
        except:
            print(f'do data in df or df_pulse {sub} {ses} but file exists')
        
    else:
        print(f'for {sub} {ses} no data.tsv')
        return None

#def read_in_pmod(sub,ses):       
#    path_file="/media/AG-Share/01_Studien/00_FOR5429_MeMoSLAP/RU_RegularMeetings/JF_Mohamed_P1/00_LOCATO/behavioral_analyses/data_extraction"
#    file='p1_merged.csv'
#    df=pd.read_csv(os.path.join(path_file,file))
#    df['ID']='sub-' + df['ID'].astype(str).str.zfill(3)
#    df['Session']='ses-' + df['Session'].astype(str)
#    df['Task_Type']=df['Task_Type'].astype(str).str.replace('1','Learning').str.replace('2','Control')
#    df=df[(df['ID'] == sub) & (df['Session'] == ses)]
#    return df
    
for sub, sess in zip(list_sub,list_ses):
    for ses in sess:

        get_project = glob.glob(os.path.join(path, sub, ses,'beh','sub-*_ses-*_task-*_acq-*_desc-data.tsv'))

        if get_project:
            project = re.search(r'task-(.{4})', get_project[0]).group(1)

        else:
            print(f'something wrong with {sub} {ses}, no data tsv file ')
        
        Vol=None
        bold=None

        path_func = os.path.join(path,sub,'ses-3','func')    #MeinzerMRI_Linux
    
        bold=glob.glob(f'{path_func}/{sub}*ses-3*run-01_bold.nii')
    
        if len(bold) !=1:
            print(f'no or to much bold images in directory {sub} ses-3')
            
        elif os.path.isfile(bold[0]):
            bold=bold[0]
        #sub='sub-039'
        #ses='ses-3'
        if (ses=='ses-4') and (bold) and (os.path.isfile(bold)):
            [Vol, bold]=get_index_of_bold(sub,'ses-3')  # to add index


        dummy=read_in_df(path,sub,ses,project,Vol)
        if dummy:
            df=dummy[0]
            acq=dummy[1]
            df_pmod = dummy[2]
        #try: 
        if (df is None) or (df_pmod is None) or (not bold) or (not os.path.isfile(bold)) or (not os.path.isfile(os.path.join(path,sub,ses,'beh',f'{sub}_{ses}_task-{project}_{acq}_desc-data.tsv'))):
            print(f'Something went wrong at subtracting the onsets {sub} {ses}')
            continue
    #except:
        print(f'should work for {sub} {ses}')
    # choose function for each experiment
        
        if not os.path.isdir(os.path.join(path_out,sub,ses)):
            try:
                os.makedirs(os.path.join(path_out,sub,ses))
            except:
                print(f'could not run mkdir on  {os.path.join(path_out,sub,ses)}')     
        
        if glob.glob(os.path.join(path, sub, ses,'beh','*p1*log')):
            

            #print('P1 not usable at the moment...')

            # Regressor is neede to build normal mat file for onsets, durations, names
            [Regressor, df_fmriprep] = Exp.Experiment.P1_onset_generating(df)

            Exp.Experiment.save_Regressor(df,sub, ses,path_out, Regressor,None)

            if os.path.isdir(os.path.join(path_BIDS_source, sub, ses,'func')):
                df_fmriprep.to_csv(os.path.join(path_BIDS_source, sub, ses,'func',f'task-{project}_events.tsv'), index = False)
            else: 
                print(f'no BIDS path for {sub} {ses} {project}')

            # Regressor_2 contains yet wrongly format of parametric modulation parameter, will be put together in first level analysis script
            #[Regressor_2, df_fmriprep] = Exp.Experiment.P1_pmod_generating(df_pmod, Regressor)
            Regressor_2 = Exp.Experiment.P1_pmod_generating(df_pmod, Regressor, pmod_col_names)

            Exp.Experiment.save_Regressor(df_pmod, sub, ses,path_out, Regressor_2,'pmod')
            
            
        
        elif glob.glob(os.path.join(path, sub, ses,'beh','*p3*log')):
            print('P1 not usable at the moment...')
            #[Regressor, df_fmriprep] = Exp.Experiment.P3_onset_generating(df)
            #Exp.Experiment.save_Regressor(df,sub, ses,path_out, Regressor,None)

            
            #if os.path.isdir(os.path.join(path_BIDS_source, sub, ses,'func')):
            #    df_fmriprep.to_csv(os.path.join(path_BIDS_source, sub, ses,'func',f'task-{project}_events.tsv'), index = False)
            #else: 
            #    print(f'no BIDS path for {sub} {ses} {project}')

            #[Regressor_2, df_fmriprep] = Exp.Experiment.P3_pmod_generating(df, Regressor_2)
            #[Regressor_2] = Exp.Experiment.P3_pmod_generating(df, Regressor_2)

            #Exp.Experiment.save_Regressor(df,sub, ses,path_out, Regressor_2,'pmod')
            
            #df_fmriprep.to_csv(os.path.join(path_BIDS_source, sub, ses,'func',f'task-{project}_events.tsv'))

        elif glob.glob(os.path.join(path, sub, ses,'beh','*p7*log')):
            #print('Tell Leo after you run this. This is not the correct script for P7!!!!')
            print('no function for P7 at the moment')
            ###Exp.Experiment.P7_onset_generating(df,sub,ses,path,path_out)
        

        else:
            print(f'somethings odd for {sub} {ses}')



