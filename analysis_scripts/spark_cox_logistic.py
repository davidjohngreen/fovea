# ----------------------------
# OCT Cox Regression Pipeline
# ----------------------------
"""
Run survival analyses testing associations between foveal imaging traits and incident
ophthalmic disease outcomes in UK Biobank.

This script sets up and runs a Cox proportional hazards pipeline in the DNAnexus /
UK Biobank RAP environment to assess whether OCT-derived foveal parameters are
associated with future risk of outcomes such as age-related macular degeneration
(AMD) or glaucoma.

Main steps:
1. Installs required Python packages and downloads the input phenotype / outcome files.
2. Connects to the UK Biobank dataset using dxdata and retrieves relevant fields,
   including dates of diagnosis, month and year of birth, year of scan, genetic sex,
   and death records.
3. Cleans and reformats the UK Biobank data, including calculation of age at scan.
4. Loads the external file containing foveal parameters and merges it with the
   UK Biobank participant data by participant ID.
5. Prepares covariates for modelling, including sex, age at scan, spherical equivalent,
   and ancestry.
6. For each foveal parameter and each disease date field, constructs a time-to-event
   dataset in which:
   - follow-up begins at the imaging date,
   - events are defined by a qualifying diagnosis after baseline,
   - censoring occurs at death or the study censor date.
7. Fits Cox proportional hazards models adjusted for the selected covariates.
8. Converts model coefficients into hazard ratios, including hazard ratios expressed
   per 1 standard deviation increase in the imaging trait.
9. Saves the final results table to CSV for downstream interpretation.

Outputs:
- A CSV file containing, for each imaging trait and outcome:
  - number of subjects
  - number of events
  - hazard ratio
  - 95% confidence intervals
  - p-value
  - standard deviation used
  - hazard ratio per 1 SD increase

Notes:
- The script uses imaging date as baseline for survival analysis.
- Diagnosis dates before baseline are excluded from the incident event definition.
- Death dates are used where available to define participant-specific censoring.
- The active section runs univariable trait-by-trait Cox models with covariate adjustment,
  while additional commented sections include alternative analyses such as multivariable
  Cox models, correlation-pruned models, FDR correction, and logistic regression.

This is useful for testing whether structural variation in foveal morphology is
associated with later ophthalmic disease risk in a large population cohort.
"""


# ----------------------------
# Install Dependencies & Download Input Files
# ----------------------------
# Install Python packages
!pip install lifelines


# Download input files from DNAnexus
!dx download "/Users/David/DaveGreen_temp/cox_regression/ICD10_codes_subset.csv"
!dx download "/Users/David/DaveGreen_temp/cox_regression/for_cox.csv"

!dx download "/Users/David/right_vs_left_checks/cox/for_cox_left.csv"



!dx download "/Users/David/right_vs_left_checks/cox/for_cox_smart_avg.csv"



# Imports
import os
import gc
import math
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import fdrcorrection
import dxpy
import dxdata
import pyspark
from pandas.api.types import is_numeric_dtype
from distutils.version import LooseVersion
import numpy as np
from math import exp

# ----------------------------
# Spark & Dataset Setup
# ----------------------------
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

# Locate the dataset programmatically
dispensed_dataset = dxpy.find_one_data_object(
    typename="Dataset",
    name="app*dataset",
    folder="/",
    name_mode="glob")

# Load dataset
dataset = dxdata.load_dataset(id=dispensed_dataset["id"])
participant = dataset['participant']

# ----------------------------
# Prepare Field List
# ----------------------------
df_outcomes = pd.read_csv("ICD10_codes_subset.csv")



df_outcomes['Field ID'] = df_outcomes['Field ID'].astype(str)
field_ids = df_outcomes['Field ID'].dropna().unique().tolist()
field_ids += ["52", "21836", "33", "34", "22001", "21836_i1", "40000"]
field_ids = [fid for fid in field_ids if fid not in ['131034', '131035', '130690', '130691', '130754', '130755', '130834', '130835', '33']]

def fields_for_id(field_id):
    matches = participant.find_fields(name_regex=rf'^p{field_id}(_i\d+)?(_a\d+)?$')
    return sorted(matches, key=lambda f: LooseVersion(f.name))

timestamp_fields = []
missing = []
for fid in field_ids:
    fs = fields_for_id(fid)
    timestamp_fields.append(fs[0]) if fs else missing.append(fid)
if missing:
    raise ValueError(f"No UKB fields matched IDs: {missing}")
timestamp_fields.append(participant.find_field(name='eid'))

engine = dxdata.connect(dialect='hive+pyspark', connect_args={'config': {'spark.kryoserializer.buffer.max': '256m','spark.sql.autoBroadcastJoinThreshold': '-1'}})
df = participant.retrieve_fields(engine=engine, fields=timestamp_fields, coding_values="replace").toPandas()

# ----------------------------
# Clean & Prepare Phenotype Data
# ----------------------------
mapping = {f"p{fld}": name for fld, name in zip(df_outcomes['Field ID'], df_outcomes['Description'])}
df.rename(columns=mapping, inplace=True)
df.rename(columns={"p52": "MOB", "p21836_i0": "YOS", "p21836_i1": "YOS2", "p34": "YOB", "p22001": "Genetic_sex"}, inplace=True)
df['YOS'] = df['YOS'].fillna(df['YOS2'])
df = df[df['YOS'].notna()]

month_dict = {month: i+1 for i, month in enumerate(['January','February','March','April','May','June','July','August','September','October','November','December'])}
df['year_scan_0'] = df['YOS'].astype(str).str.split('-').str[0]
df['month_scan_0'] = df['YOS'].astype(str).str.split('-').str[1]
df['MOB_numeric'] = df['MOB'].map(month_dict)
df['Age_at_scan'] = df['year_scan_0'].astype(int) - df['YOB'].astype(int) - (df['month_scan_0'].astype(int) < df['MOB_numeric'].astype(int))



# ----------------------------
# Merge with Input Embeddings
# ----------------------------
CFP = pd.read_csv("for_cox_left.csv")

# Trim SE outside ±6D by setting them to NaN
CFP.loc[(CFP['SE'] <= -1) | (CFP['SE'] >= 1), 'SE'] = np.nan

CFP['IID'] = CFP['patient_id'].astype(int)
df['eid'] = df['eid'].astype(int)

merged = pd.merge(df, CFP, left_on="eid", right_on="IID", how="inner")
merged['sex_binary'] = (merged['Genetic_sex'] == 'Male').astype(int)

merged['ancestry'] = pd.Categorical(
    merged['ancestry'],
    categories=['EUR', 'AFR', 'EAS', 'SAS'],
    ordered=True
)

# ----------------------------
# Prepare Columns for Cox Models
# ----------------------------
emb_cols = [c for c in merged.columns if c.startswith('embedding_')]
date_cols = [c for c in merged.columns if c.startswith('Date')]
cutoff_lower = pd.to_datetime('2006-01-01')
censor_date = pd.to_datetime('2022-05-31')
merged['YOS_date'] = pd.to_datetime(merged['YOS'], errors='coerce')
merged['date_of_death'] = pd.to_datetime(merged.get('p40000_i0'), errors='coerce')
for col in date_cols:
    merged[col] = pd.to_datetime(merged[col], errors='coerce')







# ----------------------------
# Cox with SD scaling
# ----------------------------
Z = 1.96
cox_models = []
cph = CoxPHFitter()

for emb in emb_cols:
    for date_col in date_cols:
        tmp = merged.copy()
        tmp['subject_censor'] = tmp['date_of_death'].fillna(censor_date)
        tmp = tmp[tmp['YOS_date'].notna()]
        mask = tmp[date_col].notna() & (tmp[date_col] > cutoff_lower) & (tmp[date_col] <= tmp['subject_censor'])
        tmp['event'] = mask.astype(int)
        tmp['duration'] = np.where(
            tmp['event'] == 1,
            (tmp[date_col] - tmp['YOS_date']).dt.days,
            (tmp['subject_censor'] - tmp['YOS_date']).dt.days
        )
        tmp = tmp[tmp['duration'].notna() & (tmp['duration'] >= 0)]

        ancestry_dummies = pd.get_dummies(tmp['ancestry'], drop_first=True, prefix='ancestry')

        df_model = pd.concat([
            tmp[['duration', 'event', 'sex_binary', 'Age_at_scan', 'SE', emb]],
            ancestry_dummies
        ], axis=1).dropna()

        n_subjects = len(df_model)
        n_events = int(df_model['event'].sum())
        if n_events < 5 or n_subjects < 20:
            print(f"⚠️ Skipping emb={emb}, date={date_col}: too few events/subjects ({n_events}/{n_subjects})")
            continue

        try:
            cph.fit(df_model, duration_col='duration', event_col='event')
        except ValueError as e:
            print(f"⚠️ Skipping emb={emb}, date={date_col}: {e}")
            continue

        # pull coef & SE for this predictor
        row = cph.summary.loc[emb]
        beta = float(row['coef'])
        se   = float(row['se(coef)'])
        hr   = float(row['exp(coef)'])
        lci  = float(row['exp(coef) lower 95%'])
        uci  = float(row['exp(coef) upper 95%'])
        pval = float(row['p'])

        # --- per z-score (per 1 SD) ---
        sd = float(df_model[emb].std(ddof=1))          # SD on the exact analysis set
        HR_perZ  = float(np.exp(beta * sd))
        LCI_perZ = float(np.exp((beta - Z*se) * sd))
        UCI_perZ = float(np.exp((beta + Z*se) * sd))

        # clean event name
        try:
            event_name = date_col.split('(', 1)[1].split(')', 1)[0]
        except Exception:
            event_name = str(date_col)

        cox_models.append({
            'emb': emb,
            'event_col': event_name,
            'n_subjects': n_subjects,
            'n_events': n_events,
            # native units
            'HR': hr, 'LCI_95': lci, 'UCI_95': uci, 'p_value': pval, 'SE': se, 
            # per 1 SD (z)
            'SD_used': sd,
            'HR_perZ': HR_perZ, 'LCI_95_perZ': LCI_perZ, 'UCI_95_perZ': UCI_perZ,
        })

        print(f"✅ {emb} × {event_name}: N={n_subjects}, events={n_events}, SD={sd:.4g}, HR_perZ={HR_perZ:.3g}")

# Optionally turn into a DataFrame/CSV
cox_df = pd.DataFrame(cox_models)
cox_df.to_csv("cox_results_perZ_1D.csv", index=False)












# # ----------------------------
# # Run Logistic Regression Loop
# # ----------------------------
import statsmodels.api as sm

logit_models = []

for emb in emb_cols:
    for date_col in date_cols:
        tmp = merged.copy()
        tmp['subject_censor'] = tmp['date_of_death'].fillna(censor_date)
        tmp = tmp[tmp['YOS_date'].notna()]
        mask = tmp[date_col].notna() & (tmp[date_col] > cutoff_lower) & (tmp[date_col] <= tmp['subject_censor'])
        tmp['event'] = mask.astype(int)

        ancestry_dummies = pd.get_dummies(tmp['ancestry'], drop_first=True, prefix='ancestry')

        df_model = pd.concat([
            tmp[['event', 'sex_binary', 'Age_at_scan', 'SE', emb]],
            ancestry_dummies
        ], axis=1).dropna()

        # 🔧 Convert bools to ints
        for col in df_model.columns:
            if df_model[col].dtype == bool:
                df_model[col] = df_model[col].astype(int)

        try:
            X = df_model.drop(columns='event')
            X = sm.add_constant(X)
            y = df_model['event']
            model = sm.Logit(y, X).fit(disp=0)
            smry = model.summary2().tables[1]
            sm_emb = smry.loc[emb]
        except Exception as e:
            print(f"⚠️ Logit failed for emb={emb}, date={date_col}: {e}")
            continue

        event_name = date_col.split('(', 1)[1].split(')', 1)[0]
        total_events = y.sum()

        logit_models.append({
            'emb': emb,
            'event_col': event_name,
            'total_cases': total_events,
            'OR': np.exp(sm_emb['Coef.']),
            'LCI_95': np.exp(sm_emb['[0.025']),
            'UCI_95': np.exp(sm_emb['0.975]']),
            'p_value': sm_emb['P>|z|']
        })

        print(f"✅ Fitted Logit for emb={emb}, date={date_col} ({len(df_model)} subjects, {total_events} events)")

# ----------------------------
# Save Logistic Results
# ----------------------------
df_logit_results = pd.DataFrame(logit_models)
df_logit_results.to_csv("david_logit_results.csv", index=False)


