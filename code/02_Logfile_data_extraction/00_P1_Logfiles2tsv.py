#how to use, open virtual environment and run script in there
# conda activate preproc
# python3 00_P1_Logfiles2tsv.py

#if error No module named 'polars' or similar, wrong environment, type it in correctly

import os
import re
import matplotlib.pyplot as plt
import polars as pl
import sbj_level_pl as sl
import shutil
import plots

base_path="/media/Data03/Studies/MeMoSLAP" # do not change!
output_path= "/media/Data03/Studies/MeMoSLAP/derivatives/preprocessing/SPM-standard" #change depending on your preprocessing folder

list_sub=[a for a in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, a))]

list_ses=[]
for sub in list_sub:
    list_ses.append([a for a in os.listdir(os.path.join(base_path,sub)) if os.path.isdir(os.path.join(base_path,sub, a))])


def tweak_data(data: pl.DataFrame) -> pl.DataFrame:
    response_time = get_response_time(data)
    responses = data.pipe(filter_data).pipe(get_responses)
    stimuli = data.pipe(filter_data).pipe(get_stimuli)

    return (
        stimuli.join(response_time, on="Trial", how="left")
        .join(responses, on="Trial", how="left")
        .sort("Trial")
    )


def filter_data(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data.filter(pl.col("Event Type") == "Picture")
        .drop(["Event Type"])
        .with_columns(
            pl.col("Code").str.split(by=";").list.to_struct(n_field_strategy="max_width")
        )
        .unnest("Code")
        .rename(
            {
                "field_1": "block",
                "field_2": "learning_stage",
            }
        )
    )


def get_responses(data: pl.DataFrame) -> pl.DataFrame:
    return (
        (
            data.filter(pl.col("field_0").is_in(["correct", "incorrect", "to late"]))
            .drop(["learning_stage"])
            .rename(
                {
                    "field_0": "response_type",
                    "Time": "feedback_onset",
                }
            )
        )
        .with_columns(
            pl.when(pl.col("response_type") == "to late")
            .then(pl.col("Trial") - 1)
            .otherwise(pl.col("Trial") - 2)
            .alias("Trial"),
            pl.col("response_type").cast(pl.Categorical),
        )
        .drop(["Subject", "block", "TTime"])
    )


def get_response_time(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data.filter(pl.col("Event Type") == "Response")
        .drop(["Subject", "Event Type"])
        .rename(
            {
                "Code": "response_button",
                "Time": "response_onset",
                "TTime": "response_latency",
            }
        )
        .unique(subset="Trial", keep="first")
        .with_columns(
            pl.col("response_button").cast(pl.Categorical),
        )
    )


def get_stimuli(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data.filter(
            (pl.col("block").str.starts_with("block"))
            | (pl.col("block").str.starts_with("c"))
        )
        .rename(
            {
                "field_0": "stimulus",
                "Time": "stimulus_onset",
                "Subject": "subject",
            }
        )
        .drop("TTime")
        .with_columns(
            pl.col("learning_stage").str.strip("LS").cast(pl.Int16),
            pl.col("block").cast(pl.Categorical),
            pl.col("stimulus").cast(pl.Categorical),
        )
    )

for sub, sess in zip(list_sub,list_ses):
    for ses in sess:
        curr_path=base_path+"/"+sub+"/"+ses
        main_dir= os.listdir(curr_path)
        if "beh" in main_dir:
            curr_path=base_path+"/"+sub+"/"+ses+"/"+"beh"
            main_dir= os.listdir(curr_path)
            if main_dir:
                if __name__ == "__main__":
                    for file in main_dir:
                        try:    
                            project = re.search(r'task-(.{4})', file).group(1)                         
                            version = re.search(r'acq-(.{1})', file).group(1)
                            if version =='OLMM':                            
                                bids= f'{sub}_{ses}_task-{project}_acq-{version}_beh.log'
                                print(f'{sub}_{ses}_task-{project}_acq-{version}_desc-data.tsv')
                                if  not os.path.isfile(os.path.join(output_path,sub,ses,'beh', f"{sub}_{ses}_task-{project}_acq-{version}_desc-data.tsv")): #put a 2 between t and s to overwrite  
                                    #os.path.isfile(os.path.join(base_path,sub,ses,'beh',bids)) 

                                    data = sl.read_data(
                                        os.path.join(
                                            base_path,sub,ses,'beh',
                                            bids
                                        )
                                    )

                                    bids2= f'{sub}_{ses}_task-{project}_acq-{version}'
                                    data = sl.data_from_pandas(data)
                                    d = data.pipe(sl.calibrate_time_to_mri).pipe(tweak_data)
                                    pulse = data.pipe(sl.calibrate_time_to_mri).pipe(sl.get_pulse)

                                    if not os.path.isdir(os.path.join(output_path,sub,ses,'beh')):
                                        os.makedirs(os.path.join(output_path,sub,ses,'beh'))

                                    data_bids= bids2 +"_desc-data.tsv"
                                    d.write_csv(os.path.join(output_path,sub,ses,'beh', data_bids), separator="\t")
                                    pulses_bids= bids2 +"_desc-pulses.tsv"
                                    pulse.write_csv(os.path.join(output_path,sub,ses,'beh', pulses_bids), separator="\t")
                        except:
                            print("This file could not be used:"+file)
                            #accuracy = plots.accuracy_plot_over_learning_stage(d)
                            #plt.show()

                            #responses = {"hit": 2, "miss": 2}
                            #rt = plots.rt_plot(d, responses)
                            #plt.show()
