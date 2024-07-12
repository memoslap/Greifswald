#how to use, open virtual environment and run script in there
# conda activate preproc
# python3 00_P1_Logfiles2tsv.py

#if error import polars could not be resolve or similar, wrong environment

import os
import re
import shutil
import polars as pl
import sbj_level_pl as sl

base_path="/media/Data03/Studies/MeMoSLAP" # do not change!
output_path= "/media/Data03/Studies/MeMoSLAP/derivatives/preprocessing/SPM-HGW-P7SP-clean" #change depending on your preprocessing folder
list_sub=[a for a in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, a))]

list_ses=[]
for sub in list_sub:
    list_ses.append([a for a in os.listdir(os.path.join(base_path,sub)) if os.path.isdir(os.path.join(base_path,sub, a))])

pl.enable_string_cache(True)


def get_stimuli(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data.filter(pl.col("Code").str.starts_with("Trialnum"))
        .with_columns(
            pl.col("Code").str.split(by=",").list.to_struct(n_field_strategy="max_width") # in older version arr.to_struct
        )
        .unnest("Code")
        .with_columns(
            pl.col("^field.*$").str.split(" = ").list.get(1).cast(pl.Categorical)    # in older version arr.get
        )
        .rename(
            {
                "field_0": "trialnum",
                "field_1": "face",
                "field_2": "name",
                "field_3": "congruency",
                "field_4": "correct_response",
                "Time": "stimulus_onset",
            }
        )
        .drop(["Event Type", "TTime"])
    )

def get_isi(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data.filter(pl.col("Code").str.starts_with("isi"))
        .with_columns(
            pl.col("Code").str.split("=").list.get(1).cast(pl.Float32)/1000,
            pl.col("Trial")+2   # in older version arr.get
        )
        .rename(
            {
                "Code": "isi",
            }
        )
        .drop(["Event Type", "TTime", "Time","Subject"])
    )

def tweak_data(data: pl.DataFrame) -> pl.DataFrame:
    responses = data.pipe(sl.get_responses)
    stimuli = data.pipe(get_stimuli)

    return (
        stimuli.join(responses, on="Trial", how="left")
        .rename({"Subject": "subject", "Trial": "trial"})
        .with_columns(
            pl.when(pl.col("trialnum") == "97")
            .then("p")
            .when(pl.col("trialnum") == "193")
            .then("p")
            .when(pl.col("trialnum") == "289")
            .then("p")
            .when(pl.col("correct_response") == pl.col("response"))
            .then(pl.lit(1))
            .when(pl.col("response").is_null())
            .then("m")
            .otherwise(pl.lit(0))
            .alias("correct")
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
                        
                            bids= f'{sub}_{ses}_task-{project}_acq-{version}_beh.log'
                            #print(bids)
                            if  not os.path.isfile(os.path.join(output_path,sub,ses,'beh', "*.tsv")) and os.path.isdir(os.path.join(output_path,sub,ses)): #put a 2 between t and s to overwrite  
                                #os.path.isfile(os.path.join(base_path,sub,ses,'beh',bids))  
                                data = sl.read_data(
                                    os.path.join(
                                        base_path,sub,ses,'beh',
                                        bids
                                    )
                                )
                                bids2= f'{sub}_{ses}_task-{project}_acq-{version}'
                                data = sl.data_from_pandas(data)
                                isi = data.pipe(get_isi)

                                d = data.pipe(sl.calibrate_time_to_mri).pipe(tweak_data).join(isi,right_on="Trial",how="left",left_on="trial")
                                pulses = data.pipe(sl.calibrate_time_to_mri).pipe(sl.get_pulse)
                                if not os.path.isdir(os.path.join(output_path,sub,ses,'beh')):
                                    os.makedirs(os.path.join(output_path,sub,ses,'beh'))
                                data_bids= bids2 +"_desc-data.tsv"
                                d.write_csv(os.path.join(output_path,sub,ses,'beh', data_bids), separator="\t")
                                pulses_bids= bids2 +"_desc-pulses.tsv"
                                pulses.write_csv(os.path.join(output_path,sub,ses,'beh', pulses_bids), separator="\t")
                                isi_bids= bids2 +"_desc-isi.tsv"
                                isi.write_csv(os.path.join(output_path,sub,ses,'beh', isi_bids), separator="\t")
                                print("files created for "+sub+ses)

                        except:
                            print("This file could not be used:"+file)
