import os
import numpy as np
import scipy
import pandas as pd

# parametric modulation; help in spm -> first level batch -> multiple conditions
"""
For  parametric  modulation  include  a structure array, which is up to 1 x n in size, called
pmod. n must be less than or equal to the number of cells in the names/onsets/durations
cell  arrays. The structure array pmod must have the fields: name, param and poly.  Each
of  these  fields  is  in  turn  a  cell  array  to allow the inclusion of one or more parametric
effects  per  column of the design. The field name must be a cell array containing strings.
The  field  param  is  a  cell  array  containing  a  vector  of  parameters.  Remember each
parameter  must be the same length as its corresponding onsets vector. The field poly is
a  cell  array  (for  consistency)  with  each  cell containing a single number specifying the
order of the polynomial expansion from 1 to 6.

Note  that  each  condition  is  assigned  its  corresponding  entry  in  the structure array
(condition  1  parametric  modulators  are in pmod(1), condition 2 parametric modulators
are  in  pmod(2),  etc. Within a condition multiple parametric modulators are accessed via
each  fields  cell  arrays.  So  for condition 1, parametric modulator 1 would be defined in 
pmod(1).name{1},   pmod(1).param{1},  and  pmod(1).poly{1}.  A  second  parametric
modulator  for  condition  1  would be defined as pmod(1).name{2}, pmod(1).param{2}
and  pmod(1).poly{2}.  If  there  was  also  a  parametric modulator for condition 2, then
remember  the  first  modulator  for  that  condition  is in cell array 1: pmod(2).name{1},
pmod(2).param{1},   and   pmod(2).poly{1}.   If   some,   but   not   all   conditions   are
parametrically  modulated,  then the non-modulated indices in the pmod structure can be
left  blank.  For  example,  if  conditions  1  and 3 but not condition 2 are modulated, then
specify  pmod(1)  and  pmod(3).  Similarly,  if conditions 1 and 2 are modulated but there
are 3 conditions overall, it is only necessary for pmod to be a 1 x 2 structure array.

EXAMPLE:
Make an empty pmod structure:
  pmod = struct('name',{''},'param',{},'poly',{});
Specify one parametric regressor for the first condition:
  pmod(1).name{1}  = 'regressor1';    
  pmod(1).param{1} = [1 2 4 5 6];
  pmod(1).poly{1}  = 1;
Specify 2 parametric regressors for the second condition: 
  pmod(2).name{1}  = 'regressor2-1';
  pmod(2).param{1} = [1 3 5 7]; 
  pmod(2).poly{1}  = 1;
  pmod(2).name{2}  = 'regressor2-2';
  pmod(2).param{2} = [2 4 6 8 10];
  pmod(2).poly{2}  = 1;

The  parametric  modulator  should  be  mean corrected if appropriate. Unused structure  
entries should have all fields left empty.
"""

class Experiment: 

    def matlab_onsets2fmriprep_onsets(df,onsets, durations, names_full):
        flattend_onsets=[item for sublist in onsets for item in sublist]
        flattend_durations=[item for sublist in durations for item in sublist]
        flattend_names=[item for sublist in names_full for item in sublist]

        df = pd.DataFrame({'onsets': flattend_onsets})
        df['durations'] = flattend_durations
        df['trial_type'] = flattend_names
        return df
    
    def P1_onset_generating(df):
    # 6 blocks and control only learning blocks and control

        onsets = [[],[],[],[],[],[],[],[]]
        durations = [[],[],[],[],[],[],[],[]]
        names= [[],[],[],[],[],[],[],[]]
        names_full= [[],[],[],[],[],[],[],[]]

        print(df.stimulus_onset.to_list)

        dict_learn={}
        motion_outlier=[]
        for i,numb,ascend in zip(["ls1","ls2","ls3","ls4","c1","c2","c3","c4"],[1,2,3,4,1,2,3,4],[0,1,2,3,4,5,6,7]):
                

            if i in ["ls1","ls2","ls3","ls4"]:
                dict_learn[i] =df.loc[((df['block'].str.contains('block')) & (df['learning_stage']==numb)) 
                                    #& df['response_type']=='correct'))
                            ,'stimulus_onset'].to_list()   

            elif i in ["c1","c2","c3","c4"]:
                dict_learn[i] = df.loc[((df['block'].str.contains('r')) & (df['learning_stage']==numb))
                                ,"stimulus_onset"].to_list()
                

            #dict_learn[i]=[ele for ele in dict_learn[i] if ele != []]

                
            onsets[ascend] = np.array(dict_learn[i], dtype = np.double) 
            durations[ascend]=np.array(np.ones((len(onsets[ascend]),), dtype = np.double)*4.5)
            names[ascend]=np.array(i,dtype=object)
            names_full[ascend]=np.array(list(np.repeat(i, len(onsets[ascend]))),dtype=object)  

            
       
        Regressor ={"onsets":onsets, "durations":durations, "names":names}

        df_fmriprep = Experiment.matlab_onsets2fmriprep_onsets(df,onsets, durations, names_full)
       
        return Regressor, df_fmriprep
        
    def P1_pmod_generating(df, Regressor,pmod_col_names):

        param=[[],[],[],[],[],[],[],[]]
        poly = [[],[],[],[],[],[],[],[]]
        name= [[],[],[],[],[],[],[],[]]
        names_full= [[],[],[],[],[],[],[],[]]


        #len(Regressor["onsets"][0:3])=60 -> 6 leaning stage , len(Regressor["onsets"][4:7])=20'
        # []
        dict_param={}
        for stage,num_stage in enumerate(Regressor['onsets']):  # #8 4 learning stage, 4 -> task, 4-> control
            len_R=len(Regressor['onsets'][stage])/10
            
            #for block in range(1,int(len_R)+1):
            for num_block, block in enumerate(df['Block'].unique()):          #  #6 6 blocks for each learning stage, beware control have only 2 block and can be named differently in df
            # idx is number of Stage, len is number of block
                df['RT']=df['RT'].fillna(0)
            # reminder: condition in contrast for learning stage 1 combined all 6 blocks in which learningstage is one
                if (stage < 4) and (block in ["1","2","3","4","5","6"]):    # first 4 onsets are learning, last 4 onsets are control   
                    Task_Type='learning'
                    #dict_param[f'{stage+1}_{block}']=np.ones(10)*df[(df['Block']==block) & (df['Stage']==stage+1) & (df['Task_Type']==Task_Type)]['Response_num'].values[0]
                    #dict_param[f'{stage+1}_{block}']=df[(df['Block']==block) & (df['Stage']==stage+1) & (df['Task_Type']==Task_Type)]['Response_num']
                    dict_param[f'{stage+1}_{block}']=df[(df['Block']==block) & (df['Stage']==stage+1) & (df['Task_Type']==Task_Type)][pmod_col_names]
                elif (stage >=4) and (block in ["c1","c2"]):
                    #only 2 control blocks with each 4 learning stages 
                        Task_Type='control'
                        block_num = 1 if block == "c1" else (2 if block == "c2" else None)
                        #dict_param[f'{stage+1}_{block_num}']=np.ones(10)*df[(df['Block']==block) & (df['Stage']==stage%4+1) & (df['Task_Type']==Task_Type)]['Response_num'].values[0]
                        #dict_param[f'{stage+1}_{block_num}']=df[(df['Block']==block) & (df['Stage']==stage%4+1) & (df['Task_Type']==Task_Type)]['Response_num']
                        dict_param[f'{stage+1}_{block_num}']=df[(df['Block']==block) & (df['Stage']==stage%4+1) & (df['Task_Type']==Task_Type)][pmod_col_names]

                else:
                    continue
        param_conc={}
        for stage,num in enumerate(Regressor['onsets']):
            if stage < 4:
                param_conc[f'{stage+1}']=np.concatenate((dict_param[f'{stage+1}_1'].round(3),
                                                         dict_param[f'{stage+1}_2'].round(3),
                                                         dict_param[f'{stage+1}_3'].round(3),
                                                         dict_param[f'{stage+1}_4'].round(3),
                                                         dict_param[f'{stage+1}_5'].round(3),
                                                         dict_param[f'{stage+1}_6'].round(3)),dtype = np.double)
            else:
                param_conc[f'{stage+1}']=np.concatenate((dict_param[f'{stage+1}_1'].round(3),
                                                         dict_param[f'{stage+1}_2'].round(3)),dtype = np.double)
 
            
        for ascend, name_l in enumerate(["ls1","ls2","ls3","ls4","c1","c2","c3","c4"]):

            """
            Make an empty pmod structure:
            pmod = struct('name',{''},'param',{},'poly',{});
            Specify one parametric regressor for the first condition:
            pmod(1).name{1}  = 'regressor1';    
            pmod(1).param{1} = [1 2 4 5 6];
            pmod(1).poly{1}  = 1;
            """
               
            param[ascend] = np.array(param_conc[str(ascend+1)], dtype = np.double)
            poly[ascend] = np.array(1, dtype = np.double)
            name[ascend] = np.array(name_l,dtype=object)
            names_full[ascend]=np.array(list(np.repeat(name_l, len(param[ascend]))),dtype=object) 

        Regressor_2 ={'param':param, 'poly':poly, 'name':name}
        
        Regressor_3={'pmod_wrong':Regressor_2}

        #df_fmriprep = Experiment.matlab_onsets2fmriprep_onsets(df,param, poly, names_full)
        #return Regressor_3, df_fmriprep
        return Regressor_3

    def save_Regressor(df,sub,ses,path_out,Regressor, mod):

        if ~os.path.isdir(os.path.join(path_out,sub, ses)):
            try: 
                os.mkdirs(os.path.join(path_out,sub, ses))
            except:
                print('path already exists')
        if mod==None:
            path_final=os.path.join(path_out,sub, ses,f'Onsets_Durations_Names.mat')        
        elif mod == 'pmod':
            path_final=os.path.join(path_out,sub, ses,f'Onsets_Durations_Names_pmod.mat')
        
        scipy.io.savemat(path_final, Regressor)




    def P3_onset_generating(df):
        # 6 blocks and control only learning blocks and control
        
        onsets = [[],[],[],[],[],[],[],[]]
        durations = [[],[],[],[],[],[],[],[]]
        names= [[],[],[],[],[],[],[],[]]
        names_full= [[],[],[],[],[],[],[],[]]
        print(df.stimulus_onset.to_list)

        dict_learn={}
        for i,numb,ascend in zip(["ls1","ls2","ls3","ls4","c1","c2","c3","c4"],[1,2,3,4,1,2,3,4],[0,1,2,3,4,5,6,7]):
                

            if i in ["ls1","ls2","ls3","ls4"]:
                dict_learn[i] =df.loc[((df['block'].str.contains('block')) & (df['learning_stage']==numb)) 
                                    #& df['response_type']=='correct'))
                            ,'stimulus_onset'].to_list()   

            elif i in ["c1","c2","c3","c4"]:
                dict_learn[i] = df.loc[((df['block'].str.contains('r')) & (df['learning_stage']==numb))
                                ,"stimulus_onset"].to_list()
                

            #dict_learn[i]=[ele for ele in dict_learn[i] if ele != []]

                
            onsets[ascend] = np.array(dict_learn[i], dtype = np.double) 
            durations[ascend]=np.array(np.ones((len(onsets[ascend]),), dtype = np.double)*4.5)
            names[ascend]=np.array(i,dtype=object)
            names_full[ascend]=np.array(list(np.repeat(i, len(onsets[ascend]))),dtype=object)  
   

        #Regressor = {names:"names","onsets": onsets,'durations':durations}
        Regressor ={"onsets":onsets, "durations":durations, "names":names}

        df_fmriprep = Experiment.matlab_onsets2fmriprep_onsets(df,onsets, durations, names_full)
        return Regressor, df_fmriprep
        

    def P7_onset_generating(df,sub, ses,path, path_out):


        df['congruency_s']=df['congruency'].shift(1)
        df['correct_s']=df['correct'].shift(1)
        df=df.reset_index()
        df['index_s']=df['index'].shift(1)


        dict={}
        for i in ['CC','II','CI','IC']:
            #to shift columns entries above use
            shift = 1
            if i=="II":
                df_II=df.reset_index().loc[(((df["congruency"]=='I') & (df["congruency_s"]=='I')) & ((df['correct']=="1") & (df['correct_s']=="1"))) | (((df["congruency"]=='I') & (df["congruency_s"]=='I')) &((df['correct']=="1")&(df['correct_s']=="p")))]    
                dict[i]=df_II['stimulus_onset']
                #df_II['index_s']=df_II['index'].shift(shift)
                #dict[i]=df_II.loc[df_II['index']==df_II['index_s']+shift  ,'stimulus_onset']
            elif i=='CC':
                df_CC=df.reset_index().loc[(((df["congruency"]=='C') & (df["congruency_s"]=='C')) & ((df['correct']=="1") & (df['correct_s']=="1"))) | (((df["congruency"]=='C') & (df["congruency_s"]=='C')) &((df['correct']=="1")&(df['correct_s']=="p")))] 
                #df_CC['index_s']=df_CC['index'].shift(shift)
                #dict[i]=df_CC.loc[df_CC['index']==df_CC['index_s']+shift  ,'stimulus_onset']
                dict[i]=df_CC['stimulus_onset']
            elif i=='CI':
                df_CI=df.reset_index().loc[(((df["congruency"]=='I') & (df["congruency_s"]=='C')) & ((df['correct']=="1") & (df['correct_s']=="1"))) | (((df["congruency"]=='I') & (df["congruency_s"]=='C')) &((df['correct']=="1")&(df['correct_s']=="p")))] 
                dict[i]=df_CI['stimulus_onset']
                #df_CI['index_s']=df_CI['index'].shift(shift)
                #dict[i]=df_CI.loc[df_CI['index']==df_CI['index_s']+shift  ,'stimulus_onset']
            elif i=='IC':
                df_IC=df.reset_index().loc[(((df["congruency"]=='C') & (df["congruency_s"]=='I')) & ((df['correct']=="1") & (df['correct_s']=="1"))) | (((df["congruency"]=='C') & (df["congruency_s"]=='I')) &((df['correct']=="1")&(df['correct_s']=="p")))] 
                dict[i]=df_IC['stimulus_onset']
                #df_IC['index_s']=df_IC['index'].shift(shift)
                #dict[i]=df_IC.loc[df_IC['index']==df_IC['index_s']+shift  ,'stimulus_onset']

        onsets = [[],[],[],[]]
        durations = [[],[],[],[]]
        names = [[],[],[],[]]
        names_full = [[],[],[],[]]

        for numb,i in enumerate(dict):

            onsets[numb - 1] = np.array(dict[i], dtype = np.double) 
            durations[numb - 1]=np.array(np.ones((len(onsets[numb - 1]),), dtype = np.double))
            names[numb - 1]=np.array(i,dtype=object)
            names_full[numb - 1]=np.array(list(np.repeat(i, len(onsets[numb - 1]))),dtype=object)  

        

        Regressor ={"onsets":onsets, "durations":durations, "names":names}

        #df_fmriprep = Experiment.matlab_onsets2fmriprep_onsets(df,onsets, durations, names_full)
        #return Regressor, df_fmriprep
        return Regressor
        
        #if ~os.path.isdir(os.path.join(path_out,sub, ses,)):
        #    try:
        #        os.mkdirs(os.path.join(path_out,sub, ses))
        #    except:
        #        print('path already exists')
        
        #path_final=os.path.join(path_out,sub, ses,f'Onsets_Durations_Names.mat')

        #scipy.io.savemat(path_final, Regressor)
