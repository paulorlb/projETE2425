
import os
import sys
import importlib
from pathlib import Path
import pandas as pd
import numpy as np

def fill_missing_with_neighborhood_avg(row, w, shape, colname):
    if pd.isnull(row[colname]):  # Check if the value is missing
        neighbors = w.neighbors[row.name]  # Get neighbors of the current row
        neighbor_values = shape.loc[neighbors, colname]  # Get values of neighbors
        avg_neighbor_value = np.nanmean(neighbor_values)  # Compute average of neighbors
        return avg_neighbor_value
    else:
        return row[colname]  # Return the original value if not missing



