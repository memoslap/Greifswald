#how to use, open virtual environment and run script in there
# conda activate preproc
# python3 00_P1_Logfiles2tsv.py

#if error import polars could not be resolve or similar, wrong environment

import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import polars as pl
import shutil



base_path="/media/Data03/Studies/MeMoSLAP" # do not change!
output_path= "/media/Data03/Studies/MeMoSLAP/derivatives/preprocessing/SPM-standard" #change depending on your preprocessing folder

list_sub=[a for a in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, a))]

list_ses=[]
for sub in list_sub:
    list_ses.append([a for a in os.listdir(os.path.join(base_path,sub)) if os.path.isdir(os.path.join(base_path,sub, a))])



def read_data(file: str) -> pd.DataFrame:
    cols = ["Subject", "Trial", "Event Type", "Code", "TTime", "Time"]
    return pd.read_csv(
        file,
        sep="\t",
        skiprows=3,
        usecols=cols,
        encoding="latin",
    )


def data_from_pandas(data: pd.DataFrame) -> pl.DataFrame:
    data = pl.from_pandas(data)

    return data.filter(pl.col("Subject") == data[0, "Subject"]).with_columns(
        pl.col("Subject").cast(pl.Categorical),
        pl.col("Trial").cast(pl.Int16),
        pl.col("Event Type").cast(pl.Categorical),
        pl.col("Code").cast(pl.Utf8),
        pl.col("TTime").cast(pl.Int32),
        pl.col("Time").cast(pl.Int32),
    )


def get_stimuli(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data
        # only include presented stimuli
        .filter(pl.col("Event Type") == "Picture")
        # drop unnecessary columns
        .drop(["Event Type", "TTime"])
        # split the 'Code' column
        .with_columns(
            pl.col("Code").str.split(by=";")
            # create a new column for each instance of split
            .list.to_struct(n_field_strategy="max_width")
        )
        .unnest("Code")
        # give the created split instances new names
        .rename(
            {
                "field_0": "condition",
                "field_1": "word",
                "field_2": "picture",
                "field_3": "coherence",
                "field_4": "block",
                "field_5": "learning_stage",
                "Time": "stimulus_onset",
            }
        )
        .drop_nulls()  # drops the feedback trials
        .with_columns(
            # set the data type of the new columns
            pl.col("word").cast(pl.Categorical),
            pl.col("picture").cast(pl.Categorical),
            pl.col("coherence").cast(pl.Categorical),
            pl.col("block").str.strip("block").cast(pl.Int16),
            pl.col("learning_stage").str.strip("LS").cast(pl.Int16),
            pl.col("condition")
            .str.split("/")
            .list.to_struct(n_field_strategy="max_width"),
        )
        .unnest("condition")
        .rename(
            {
                "field_0": "condition",
                "field_1": "picture_path",
            }
        )
        .with_columns(
            pl.col("condition").cast(pl.Categorical),
            pl.col("picture_path").cast(pl.Categorical),
        )
    )


def get_responses(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data
        # only include responses
        .filter(pl.col("Event Type") == "Response")
        # drop unnecessary columns
        .drop(["Subject", "Event Type"])
        .unique(subset="Trial", keep="first")
        .with_columns(pl.col("Code").cast(pl.Categorical))
        .rename(
            {
                "Code": "response",
                "TTime": "reaction_time",
                "Time": "response_onset",
            }
        )
    )


def get_pulse(data: pl.DataFrame) -> pl.DataFrame:
    pulses = data.filter(pl.col("Event Type") == "Pulse")
    return pulses


def calibrate_time_to_mri(data: pl.DataFrame) -> pl.DataFrame:
    return data.with_columns(
        pl.col("Time")
        - data.filter(pl.col("Event Type") == "Pulse")["Time"][0]
    )


def tweak_data(data: pl.DataFrame) -> pl.DataFrame:
    _data = data.filter(pl.col("Code") != "fix")

    stimuli = get_stimuli(_data)
    response = get_responses(_data)

    _data = (
        stimuli.join(response, on="Trial", how="left")
        .with_columns(
            # Signal detection theory classifies responses as hit,
            # miss, false alarm, or correct rejection
            signal_detection_theory=pl.when(
                (pl.col("coherence") == "c") & (pl.col("response") == "yes")
            )
            .then(pl.lit("hit"))
            .when((pl.col("coherence") == "c") & (pl.col("response") == "no"))
            .then(pl.lit("miss"))
            .when((pl.col("coherence") == "i") & (pl.col("response") == "yes"))
            .then(pl.lit("false alarm"))
            .when((pl.col("coherence") == "i") & (pl.col("response") == "no"))
            .then(pl.lit("correct rejection"))
        )
        .with_columns(
            # create a column 'hit' with 0 for false responses and
            # 1 for correct responses
            hit=pl.when(
                pl.col("signal_detection_theory").is_in(
                    ["hit", "correct rejection"]
                )
            )
            .then(pl.lit(1))
            .when(
                pl.col("signal_detection_theory").is_in(
                    ["miss", "false alarm"]
                )
            )
            .then(pl.lit(0))
        )
        .drop("Trial")
    )

    return _data

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
                        
                            bids= f'{sub}_{ses}_task-{project}_acq-{version}_beh.log'
                            #print(bids)
                            if  not os.path.isfile(os.path.join(output_path,sub,ses,'beh', "*.tsv")): #put a 2 between t and s to overwrite  
                                #os.path.isfile(os.path.join(base_path,sub,ses,'beh',bids)) 
                                #if os.path.isfile(os.path.join(base_path,sub,ses,'beh',f'{sub}p3-{ses}_task-learning.log')) and (not os.path.isfile(os.path.join(base_path,sub,ses,'beh', "data.tsv"))):
                                
                                data = read_data(os.path.join(base_path,sub,ses,'beh',
                                            bids))
                                
                                bids2= f'{sub}_{ses}_task-{project}_acq-{version}'
                                data = data_from_pandas(data)
                                d = data.pipe(calibrate_time_to_mri).pipe(tweak_data)
                                pulses = data.pipe(calibrate_time_to_mri).pipe(get_pulse)

                                if not os.path.isdir(os.path.join(output_path,sub,ses,'beh')):
                                    os.makedirs(os.path.join(output_path,sub,ses,'beh'))

                                data_bids= bids2 +"_desc-data.tsv"
                                d.write_csv(os.path.join(output_path,sub,ses,'beh',data_bids), separator="\t")
                                pulses_bids= bids2 +"_desc-pulses.tsv"
                                pulses.write_csv(os.path.join(output_path,sub,ses,'beh',pulses_bids), separator="\t")
                        except:
                            print("This file could not be used:"+file)
                        
