import os

import matplotlib.pyplot as plt
import pandas as pd
import polars as pl

import plots

print(pl.__version__)


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
            {	 "field_0": "condition",
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


if __name__ == "__main__":
	path=os.getcwd()
	data = read_data(os.path.join(path,"pilot13_P7_Task.log"))
	#data = read_data(os.path.join(path,"pilot06_P7_Task.log"))
	
	data = data_from_pandas(data)
	d = data.pipe(calibrate_time_to_mri).pipe(tweak_data)
	pulses = data.pipe(calibrate_time_to_mri).pipe(get_pulse)

	d.write_csv("data.tsv", separator="\t")
	pulses.write_csv("pulses.tsv", separator="\t")

	accuracy = plots.accuracy_plot_over_learning_stage(d)
	plt.show()

	responses = {"correct": 1, "miss": 2}
	rt = plots.rt_plot(d, responses)
	plt.show()
