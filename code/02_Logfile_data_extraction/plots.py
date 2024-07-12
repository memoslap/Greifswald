import numpy as np
import polars as pl
import seaborn as sns


def accuracy_plot_over_learning_stage(data: pl.DataFrame):
    accuracy = sns.barplot(
        data=data
        # ignore trials where no answer was given
        .drop_nulls().to_pandas(),
        x="learning_stage",
        y="hit",
        estimator=np.mean,
        hue="condition",
        alpha=0.5,
        edgecolor="black",
        capsize=0.2,
    )

    sns.move_legend(
        accuracy,
        "upper left",
        bbox_to_anchor=(1, 1),
    )
    accuracy.set(
        xlabel="Learning stage",
        ylabel="proportion correct",
    )

    return accuracy


def rt_plot(data, response):
    rt = sns.histplot(
        data=data.filter(pl.col("hit") == response["hit"]).to_pandas(),
        x="reaction_time",
        edgecolor="black",
        hue="condition",
    )

    sns.move_legend(
        rt,
        "upper left",
        bbox_to_anchor=(1, 1),
    )
    rt.set(
        xlabel="Reaction time (0.1 * ms)",
    )

    return rt


def rt_plot_over_learning_stage(data, response):
    rt_over_learning_stage = sns.barplot(
        data=data.filter(pl.col("hit") == response["hit"]).to_pandas(),
        x="learning_stage",
        y="reaction_time",
        estimator=np.mean,
        hue="condition",
        alpha=0.5,
        edgecolor="black",
        capsize=0.2,
    )

    sns.move_legend(
        rt_over_learning_stage,
        "upper left",
        bbox_to_anchor=(1, 1),
    )
    rt_over_learning_stage.set(
        xlabel="Learning Stage",
        ylabel="Reaction time (0.1 * ms)",
    )

    return rt_over_learning_stage
