# Import relevant libraries
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math


df = pd.read_csv("./Step_1_MM20230105_out_base_data.csv")

def group_dx(group: pd.DataFrame) -> pd.DataFrame:
    group.sort_values("timepoint", inplace=True)
    first_timepoint = group.loc[group["timepoint"] == group["timepoint"].min()]
    first_rlu = first_timepoint["rlu_value"].values[0]
    group["total_fold_change"] = [x / first_rlu for x in group["rlu_value"]]
    group["inter_timepoint_fold_change"] = [1 + x for x in group["rlu_value"].pct_change().fillna(0)]
    return group

df = df.groupby(["sample_id"], as_index=False, group_keys=False).apply(group_dx)
df.head()


# ASSAY_ID = "plate035"
# EXCLUDE = [
#     "plate035_A_1"
# ]
# filtered_df = df.loc[df["assay_id"] == ASSAY_ID]
# filtered_df = filtered_df[~filtered_df["sample_id"].isin(EXCLUDE)]

# Toggle comments for below lines to run visualization of the plate grid

#grid = sns.FacetGrid(data=filtered_df, col="well_column", row="well_row")
#grid.map_dataframe(sns.lineplot, x="timepoint", y="inter_timepoint_fold_change")

clean_df = df[~df["sample_id"].isin(EXCLUDE)]
print(clean_df.head())


CONTROL_ID = "D36E_MQ"

ctrl_df = clean_df.loc[df["treatment_id"] == CONTROL_ID]
ctrl_group = ctrl_df.groupby(["assay_id", "timepoint"], as_index=False, group_keys=False)\
    .agg(
        inter_timepoint_fold_change_ctrl=("inter_timepoint_fold_change", "mean"),
    )
normalised_df = clean_df.merge(
    ctrl_group[[
        "assay_id",
        "timepoint",
        "inter_timepoint_fold_change_ctrl",
    ]],
    on=["assay_id", "timepoint"]
)


normalised_df["inter_timepoint_fold_change"] = normalised_df["inter_timepoint_fold_change"] / normalised_df["inter_timepoint_fold_change_ctrl"]
intertimepoint = [x for x in normalised_df["inter_timepoint_fold_change"] if x <=0]

print(normalised_df["inter_timepoint_fold_change"])

normalised_df["inter_timepoint_fold_change"] = [math.log2(x) for x in normalised_df["inter_timepoint_fold_change"]]
normalised_df.head(100)

PLOT_VALUE = "inter_timepoint_fold_change"
PLOT_ASSAY = "plate018"

fig, axes = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(15, 5))
assay_df = normalised_df[normalised_df["assay_id"] == PLOT_ASSAY]
for idx, (id, treatment_df) in enumerate(assay_df.groupby("treatment_id")):
    sns.lineplot(
        data=treatment_df,
        x="timepoint",
        y=PLOT_VALUE,
        hue="treatment_id",
        ax=axes[math.floor(idx / 2), idx % 2]
    )


TIMEPOINT = 50

def auc_calc(group):
    group.sort_values("timepoint", inplace=True)
    auc = np.trapz(group["inter_timepoint_fold_change"], dx=2.5)
    return auc


auc_df = normalised_df[normalised_df["timepoint"] <= TIMEPOINT]
auc_df = auc_df.groupby(["treatment_id", "assay_id", "sample_id"], as_index=False, group_keys=False).apply(auc_calc)
auc_df.rename({None: "area_under_curve"}, inplace=True, axis="columns")
auc_df.to_csv("Step_2_MM20230106_2_auc_data.csv") # change name/update date
auc_df




sns.set(rc={"figure.figsize":(100, 30)})

boxplot = sns.boxplot(
    data=auc_df,
    x="treatment_id",
    y="area_under_curve",
    hue="treatment_id",
    width=10,
)
boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation=20)
boxplot
