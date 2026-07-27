#!/usr/bin/env python3
"""Simple data analysis demo using pandas and numpy"""

from pathlib import Path
import pandas as pd
import numpy as np

# Setup paths
script_dir = Path(__file__).resolve().parent
out_dir = script_dir.parent / "code_out" / "demo_analysis"
out_dir.mkdir(parents=True, exist_ok=True)

# Generate sample experimental data
np.random.seed(42)
n_samples = 20

data = {
    'sample_id': [f'Sample_{i:02d}' for i in range(1, n_samples + 1)],
    'group': ['Control'] * 10 + ['Treatment'] * 10,
    'measurement_1': np.random.normal(100, 15, n_samples),
    'measurement_2': np.random.normal(50, 10, n_samples)
}

# Create DataFrame
df = pd.DataFrame(data)

# Calculate summary statistics
df['mean_measurement'] = df[['measurement_1', 'measurement_2']].mean(axis=1)
df['ratio'] = df['measurement_1'] / df['measurement_2']

# Group statistics
summary = df.groupby('group').agg({
    'measurement_1': ['mean', 'std'],
    'measurement_2': ['mean', 'std'],
    'mean_measurement': ['mean', 'std']
}).round(2)

# Save results
df.to_csv(out_dir / "sample_measurements.csv", index=False)
summary.to_csv(out_dir / "group_summary.csv")

print(f"Done! Results in: {out_dir}")
print("\nSample data (first 5 rows):")
print(df.head().to_string(index=False))
print("\nGroup summary statistics:")
print(summary)
