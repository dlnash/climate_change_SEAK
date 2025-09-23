"""
Filename:    time_helpers.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: scripts dealing with time

"""
import re
import pandas as pd

def find_date_based_on_filename(fname):
    numbers = re.findall(r'\d+', fname)
    df = pd.DataFrame([numbers[2:]], columns=['Year', 'Month', 'Day'])
    date = pd.to_datetime(df)
    date = date.values

    return date

