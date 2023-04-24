"""
This file is solely concerned with reading out directory of plate files with well data
and then merging them with a layout file.
The produced artifact is a csv file that contains all well data in long format
with the metadata from the layout file merged into it.
"""
import os
import math
import numpy as np
import pandas as pd

# the working directory is wherever the script is saved
# Script configuration
LAYOUT_FILE = "./MM20221209_summary_ROS_DC3000_ef_library.xlsx"
PLATES_DIR = "." # enter "." if it is in the same directory as the script
OUTPUT_PREFIX = "MM20230106_2_out"

# In the strings input a way to identify the plate id from a filename
# if you have a filename like MM_ROS_plate01_processed.xslx
# then the prefix would be 'MM_ROS_' and the suffix would be '_processed'
# To select a prefix and suffix, identify the common and constant parts in filenames 

PLATE_ID_PREFIX = "ROS_"
PLATE_ID_SUFFIX = "_processed"
PLATE_ROW_COUNT = 8

# Code
def get_plate_files(plates_dir: str) -> list[str]:
    """
    Scan the given directory and select all files that containing the prefixes for plate ID identification.

    TODO: simplify
    """
    return [
        os.path.join(plates_dir, f) for f in os.listdir(plates_dir) if PLATE_ID_PREFIX in f and PLATE_ID_SUFFIX in f
    ]

def load_plate_df(files: list[str]) -> list[pd.DataFrame]:
    """
    Load all plate well files into a list of dataframes
    Additionally:
    - Add a column for plate_id derived from the filename
    - Drop all empty rows
    - Add a column for timepoint derived from the row index

    TODO: document partitioning
    """
    frames = []
    for f in files:
        df = pd.read_excel(f)
        _, _, suffix = f.partition(PLATE_ID_PREFIX)
        plate_id, _, _ = suffix.partition(PLATE_ID_SUFFIX)
        df["plate_id"] = plate_id
        df["Col1"].replace("", np.nan, inplace=True)
        df.dropna(subset=["Col1"], inplace=True)
        relative_times = (row_index / PLATE_ROW_COUNT for row_index in df.index)
        absolute_times = (math.ceil(value) if not value.is_integer() else math.ceil(value) + 1 for value in relative_times)
        timepoints = [abs_time * 2.5 for abs_time in absolute_times]
        df["time_point"] = timepoints
        frames.append(df)
    return frames


def melt_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """
    First, combine all the individual frames into a single dataframe
    Then change the orientation of well plate data from wide to long
    Additionally:
    - Add a column indicating the well position
    - Rename the measured value to 'rlu' (Relative Light Unit, an arbitrary number depending on machine and sample)
    """
    merged = pd.concat(frames, ignore_index=True)
    melted = pd.melt(
        merged,
        id_vars=["Row", "plate_id", "time_point"],
        value_vars=['Col1', 'Col2', 'Col3', 'Col4', 'Col5', 'Col6', 'Col7', 'Col8', 'Col9', 'Col10', 'Col11', 'Col12']
    )
    melted["variable"] = melted["variable"].str.replace("Col", "")
    melted["Well"] = melted["Row"].str.cat(melted["variable"], sep=":")
    melted.drop("variable", axis=1, inplace=True)
    melted.rename(columns={"value": "rlu"}, inplace=True)
    return melted

def attach_layout(plate_df: pd.DataFrame, layout_file: str) -> pd.DataFrame:
    """
    Read the layout file and merge it with the plate dataframe
    """
    layout_df = pd.read_excel(layout_file)
    layout_df.rename(columns={"PlateID": "plate_id"}, inplace=True)
    return pd.merge(layout_df, plate_df, on=["Well", "plate_id"])

def rename_cols(df: pd.DataFrame) -> pd.DataFrame:
    df["well_row"] = df["Row"]
    df["well_column"] = df["Well"].str[2:].astype(int)
    df["assay_id"] = df["plate_id"]
    df["treatment_id"] = df["Treatment"].str.cat(df["Mock"], sep="_")
    df["rlu_value"] = df["rlu"]
    df["sample_id"] = df["assay_id"].str.cat(df["well_row"], sep="_").str.cat(df["well_column"].astype(str), sep="_")
    df["timepoint"] = df["time_point"]
    return df[[
        "assay_id",
        "treatment_id",
        "sample_id",
        "well_row",
        "well_column",
        "rlu_value",
        "timepoint"
    ]]

plate_files = get_plate_files(PLATES_DIR)
plate_frames = load_plate_df(plate_files)
plate_df = melt_frames(plate_frames)
df = attach_layout(plate_df, LAYOUT_FILE)
df = rename_cols(df)
df.to_csv(f"{OUTPUT_PREFIX}_base_data.csv")