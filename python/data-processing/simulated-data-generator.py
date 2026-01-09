# ---------------------------------------------
# Title: Simulated Clinical/Genomic Data Generator
# Description: Generate simulated patient cohorts with mutations
# Input: Configuration parameters
# Output: Pandas DataFrames for each data type
# Source: NBCT_Gene_Analysis
# ---------------------------------------------

import pandas as pd
import numpy as np
import random
import uuid
from datetime import datetime, timedelta


def random_date(start_date, end_date):
    """Generate random date between start and end."""
    time_between = end_date - start_date
    random_days = random.randrange(time_between.days)
    return (start_date + timedelta(days=random_days)).strftime('%Y%m%d')


def generate_patient_cohorts(n_patients, mutation_rates):
    """
    Generate patient IDs and stratify into mutation cohorts.
    
    Args:
        n_patients: Total number of patients
        mutation_rates: Dict with keys 'cna_rate', 'snv_rate', 'other_rate', 'death_rate'
    
    Returns:
        Dict with patient ID lists for each cohort
    """
    patient_ids = [str(uuid.uuid4()) for _ in range(n_patients)]
    random.shuffle(patient_ids)
    
    # Mutation cohorts (mutually exclusive)
    n_cna = int(n_patients * mutation_rates['cna_rate'])
    n_snv = int(n_patients * mutation_rates['snv_rate'])
    n_other = int(n_patients * mutation_rates['other_rate'])
    n_death = int(n_patients * mutation_rates['death_rate'])
    
    unassigned = list(patient_ids)
    
    cna_ids = random.sample(unassigned, k=n_cna)
    unassigned = [p for p in unassigned if p not in cna_ids]
    
    snv_ids = random.sample(unassigned, k=n_snv)
    
    # Other mutations can overlap
    other_ids = random.sample(patient_ids, k=n_other)
    
    # Death cohort
    deceased_ids = random.sample(patient_ids, k=n_death)
    
    return {
        'all': patient_ids,
        'cna': cna_ids,
        'snv': snv_ids,
        'other': other_ids,
        'deceased': deceased_ids
    }


def generate_demographics(patient_ids, columns):
    """Generate patient demographic data."""
    data = []
    for pid in patient_ids:
        data.append({
            columns['key']: pid,
            columns['year_of_birth']: str(random.randint(1940, 1980)),
            columns['gender']: str(random.choices([1, 2], weights=[0.55, 0.45], k=1)[0]),
            columns['height']: str(random.randint(150, 185)),
            columns['weight']: str(random.randint(50, 95))
        })
    return pd.DataFrame(data)


def generate_survival_data(patient_ids, deceased_ids, columns, max_days=1095):
    """Generate death and follow-up data."""
    death_data = []
    followup_data = []
    
    death_dict = {}
    for pid in deceased_ids:
        death_day = random.randint(30, max_days)
        death_dict[pid] = death_day
        death_data.append({columns['key']: pid, columns['death_date']: str(death_day)})
    
    for pid in patient_ids:
        if pid in death_dict:
            followup_day = death_dict[pid]
        else:
            followup_day = random.randint(max_days - 100, max_days)
        followup_data.append({columns['key']: pid, columns['followup_date']: str(followup_day)})
    
    return pd.DataFrame(death_data), pd.DataFrame(followup_data)
