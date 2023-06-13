import pandas as pd
import numpy as np
import os
import math

layout = pd.read_excel("C:/Users/Mende012/Documents/Bioinformatics/ROS assays/XZ2022_assay_optimisation/XZ20211117_plate006_layout.xlsx")
plate = pd.read_excel("C:/Users/Mende012/Documents/Bioinformatics/ROS assays/XZ2022_assay_optimisation/XZ20211117_plate006_RLU_processed.xlsx")

plate["Col1"].replace("", np.nan, inplace=True)
plate.dropna(subset=["Col1"], inplace=True)
relative_times = (row_index / 8 for row_index in plate.index)
absolute_times = (math.ceil(value) if not value.is_integer() else math.ceil(value) + 1 for value in relative_times)
plate["timepoint"] = [abs_time * 2.5 for abs_time in absolute_times]
plate["assay_id"] = "plate006"
platemelt = pd.melt(
    plate,
    id_vars=["Row", "assay_id", "timepoint"],
    value_vars=['Col1', 'Col2', 'Col3', 'Col4', 'Col5', 'Col6', 'Col7', 'Col8', 'Col9', 'Col10', 'Col11', 'Col12']
)

layoutmelt = pd.melt(
    layout,
    id_vars=["Row"],
    value_vars=['Col1', 'Col2', 'Col3', 'Col4', 'Col5', 'Col6', 'Col7', 'Col8', 'Col9', 'Col10', 'Col11', 'Col12']
)
layoutmelt.rename(columns={"value": "treatment_id"}, inplace=True)
merged = platemelt.merge(layoutmelt, on=["Row", "variable"])
merged.rename(columns={"Row": "well_row", "value": "rlu_value"}, inplace=True)
merged["well_column"] = merged["variable"].str[3:].astype(int)
merged["sample_id"] = merged["assay_id"].str.cat(merged["well_row"], sep="_").str.cat(merged["well_column"].astype(str), sep="_")
merged.drop(columns=["variable"], inplace=True)
merged["rlu_value"] = merged["rlu_value"].astype(int)
merged.to_csv("MM20230116_plate006_long.csv")