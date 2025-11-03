"""
Quality Control Analysis Example
=================================

This script demonstrates how to analyze QC results, generate statistics,
and create filtered datasets based on quality flags.

Author: Benedetta Torelli
Repository: https://github.com/BennyTorelli/Seaexplorer_pyglider
Date: October 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================================
# STEP 1: Load QC Data
# ============================================================================

print("=" * 70)
print("QUALITY CONTROL ANALYSIS")
print("=" * 70)

qc_file = 'output/analysis/seaexplorer_qc_variables.csv'

if not os.path.exists(qc_file):
    print(f"ERROR: QC file not found: {qc_file}")
    print("Please run: python scripts/qc_variables.py --latest")
    exit(1)

print(f"\nLoading: {qc_file}")
qc = pd.read_csv(qc_file, parse_dates=['TIME'])
print(f"Loaded {len(qc):,} measurements")

# ============================================================================
# STEP 2: QC Summary Statistics
# ============================================================================

print("\n" + "=" * 70)
print("QC SUMMARY STATISTICS")
print("=" * 70)

# Get all QC columns
qc_cols = [c for c in qc.columns if '_QC' in c]

print(f"\nTotal QC tests: {len(qc_cols)}")
print(f"Total measurements: {len(qc):,}\n")

# Count flags for each QC column
for col in qc_cols:
    counts = qc[col].value_counts().sort_index()
    total = len(qc)
    
    print(f"{col}:")
    for flag, count in counts.items():
        pct = (count / total) * 100
        flag_name = {1: 'GOOD', 4: 'BAD', 9: 'MISSING'}.get(flag, 'UNKNOWN')
        print(f"  {flag_name:8} ({flag}): {count:>10,} ({pct:>6.2f}%)")
    print()

# ============================================================================
# STEP 3: Identify Problem Variables
# ============================================================================

print("=" * 70)
print("PROBLEM VARIABLE IDENTIFICATION")
print("=" * 70)

print("\nVariables with >5% BAD flags:\n")

range_qc_cols = [c for c in qc_cols if 'Range_QC' in c]

for col in range_qc_cols:
    bad_count = (qc[col] == 4).sum()
    bad_pct = (bad_count / len(qc)) * 100
    
    if bad_pct > 5.0:
        var_name = col.replace('_Range_QC', '')
        print(f"  {var_name:6} : {bad_count:>8,} BAD flags ({bad_pct:>5.1f}%)")
        print(f"           Recommendation: Expert review required")
        print()

# ============================================================================
# STEP 4: Create Filtered Datasets
# ============================================================================

print("=" * 70)
print("CREATING FILTERED DATASETS")
print("=" * 70)

# Filter 1: STRICT (all QC flags = 1)
strict_filter = (qc[qc_cols] == 1).all(axis=1)
qc_strict = qc[strict_filter]
strict_file = 'output/analysis/seaexplorer_data_QC_strict.csv'
qc_strict.to_csv(strict_file, index=False)
print(f"\n1. STRICT filter:")
print(f"   Criteria: ALL QC flags = GOOD (1)")
print(f"   Retained: {len(qc_strict):,} / {len(qc):,} ({len(qc_strict)/len(qc)*100:.1f}%)")
print(f"   Saved to: {strict_file}")

# Filter 2: MODERATE (Range_QC and Sensor_QC = 1, allow MISSING)
range_qc_cols = [c for c in qc_cols if 'Range_QC' in c]
sensor_qc_cols = [c for c in qc_cols if 'Sensor_QC' in c]

moderate_filter = (
    (qc[range_qc_cols] == 1).all(axis=1) &
    (qc[sensor_qc_cols].isin([1, 9])).all(axis=1)
)
qc_moderate = qc[moderate_filter]
moderate_file = 'output/analysis/seaexplorer_data_QC_moderate.csv'
qc_moderate.to_csv(moderate_file, index=False)
print(f"\n2. MODERATE filter (RECOMMENDED):")
print(f"   Criteria: Range_QC and Sensor_QC = GOOD, allow MISSING")
print(f"   Retained: {len(qc_moderate):,} / {len(qc):,} ({len(qc_moderate)/len(qc)*100:.1f}%)")
print(f"   Saved to: {moderate_file}")

# Filter 3: LIBERAL (only Sensor_QC = 1 or 9)
liberal_filter = (qc[sensor_qc_cols].isin([1, 9])).all(axis=1)
qc_liberal = qc[liberal_filter]
liberal_file = 'output/analysis/seaexplorer_data_QC_liberal.csv'
qc_liberal.to_csv(liberal_file, index=False)
print(f"\n3. LIBERAL filter:")
print(f"   Criteria: Sensor_QC = GOOD or MISSING (excludes instrument failures)")
print(f"   Retained: {len(qc_liberal):,} / {len(qc):,} ({len(qc_liberal)/len(qc)*100:.1f}%)")
print(f"   Saved to: {liberal_file}")

# ============================================================================
# STEP 5: Temporal QC Pattern Analysis
# ============================================================================

print("\n" + "=" * 70)
print("TEMPORAL QC PATTERN ANALYSIS")
print("=" * 70)

# Resample by day and count BAD flags
qc['DATE'] = qc['TIME'].dt.date

daily_stats = []
for date in sorted(qc['DATE'].unique()):
    day_data = qc[qc['DATE'] == date]
    
    bad_temp = (day_data['TEMP_Range_QC'] == 4).sum()
    bad_doxy = (day_data['DOXY_Range_QC'] == 4).sum()
    bad_chla = (day_data['CHLA_Range_QC'] == 4).sum()
    
    daily_stats.append({
        'Date': date,
        'Total': len(day_data),
        'TEMP_BAD': bad_temp,
        'DOXY_BAD': bad_doxy,
        'CHLA_BAD': bad_chla
    })

daily_df = pd.DataFrame(daily_stats)

print("\nDaily BAD flag counts (first 5 days):")
print(daily_df.head(5).to_string(index=False))

# Check for increasing trends (biofouling indicator)
print("\n\nChecking for biofouling patterns...")
if len(daily_df) >= 7:
    early_chla = daily_df.head(7)['CHLA_BAD'].mean()
    late_chla = daily_df.tail(7)['CHLA_BAD'].mean()
    
    if late_chla > early_chla * 2:
        print(f"⚠️  WARNING: Chlorophyll BAD flags increased from {early_chla:.1f} to {late_chla:.1f} per day")
        print("   This may indicate biofouling on bio-optical sensors")
    else:
        print(f"✓ No significant biofouling detected (early={early_chla:.1f}, late={late_chla:.1f} BAD/day)")

# ============================================================================
# STEP 6: Generate QC Visualization
# ============================================================================

print("\n" + "=" * 70)
print("GENERATING QC VISUALIZATIONS")
print("=" * 70)

# Create figure with 3 subplots
fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

# Subsample for faster plotting (every 100th point)
qc_plot = qc.iloc[::100]

# Temperature with QC flags
scatter1 = axes[0].scatter(qc_plot['TIME'], qc_plot['TEMP'], 
                          c=qc_plot['TEMP_Range_QC'], 
                          cmap='RdYlGn_r', s=1, vmin=1, vmax=9)
axes[0].set_ylabel('Temperature (°C)', fontsize=11)
axes[0].set_title('Temperature QC (Green=Good, Red=Bad, Gray=Missing)', fontsize=12)
axes[0].grid(True, alpha=0.3)

# Oxygen with QC flags
scatter2 = axes[1].scatter(qc_plot['TIME'], qc_plot['DOXY'], 
                          c=qc_plot['DOXY_Range_QC'], 
                          cmap='RdYlGn_r', s=1, vmin=1, vmax=9)
axes[1].set_ylabel('Oxygen (µmol/kg)', fontsize=11)
axes[1].set_title('Dissolved Oxygen QC', fontsize=12)
axes[1].grid(True, alpha=0.3)

# Chlorophyll with QC flags
scatter3 = axes[2].scatter(qc_plot['TIME'], qc_plot['CHLA'], 
                          c=qc_plot['CHLA_Range_QC'], 
                          cmap='RdYlGn_r', s=1, vmin=1, vmax=9)
axes[2].set_ylabel('Chlorophyll (mg/m³)', fontsize=11)
axes[2].set_xlabel('Time', fontsize=11)
axes[2].set_title('Chlorophyll QC', fontsize=12)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()

output_file = 'output/analysis/qc_timeseries_overview.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ Visualization saved: {output_file}")

plt.close()

# ============================================================================
# STEP 7: Create QC Report
# ============================================================================

print("\n" + "=" * 70)
print("GENERATING QC REPORT")
print("=" * 70)

report_file = 'output/analysis/QC_Report.txt'

with open(report_file, 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("QUALITY CONTROL REPORT\n")
    f.write("=" * 70 + "\n\n")
    
    f.write(f"Generated: {pd.Timestamp.now()}\n")
    f.write(f"Dataset: {qc_file}\n")
    f.write(f"Total measurements: {len(qc):,}\n\n")
    
    f.write("SUMMARY BY QC TEST\n")
    f.write("-" * 70 + "\n\n")
    
    for col in qc_cols:
        counts = qc[col].value_counts().sort_index()
        f.write(f"{col}:\n")
        for flag, count in counts.items():
            pct = (count / len(qc)) * 100
            flag_name = {1: 'GOOD', 4: 'BAD', 9: 'MISSING'}.get(flag, 'UNKNOWN')
            f.write(f"  {flag_name:8} ({flag}): {count:>10,} ({pct:>6.2f}%)\n")
        f.write("\n")
    
    f.write("\n" + "=" * 70 + "\n")
    f.write("FILTERED DATASETS CREATED\n")
    f.write("=" * 70 + "\n\n")
    
    f.write(f"1. STRICT:    {len(qc_strict):>10,} measurements ({len(qc_strict)/len(qc)*100:.1f}%)\n")
    f.write(f"   File: {strict_file}\n\n")
    
    f.write(f"2. MODERATE:  {len(qc_moderate):>10,} measurements ({len(qc_moderate)/len(qc)*100:.1f}%)\n")
    f.write(f"   File: {moderate_file}\n\n")
    
    f.write(f"3. LIBERAL:   {len(qc_liberal):>10,} measurements ({len(qc_liberal)/len(qc)*100:.1f}%)\n")
    f.write(f"   File: {liberal_file}\n\n")
    
    f.write("\n" + "=" * 70 + "\n")
    f.write("RECOMMENDATIONS\n")
    f.write("=" * 70 + "\n\n")
    
    f.write("For most scientific applications, use MODERATE filtered dataset.\n")
    f.write("For high-precision studies, use STRICT filtered dataset.\n")
    f.write("For visualization and exploratory analysis, LIBERAL filter is acceptable.\n\n")
    
    # Check for issues
    issues = []
    for col in range_qc_cols:
        bad_pct = (qc[col] == 4).sum() / len(qc) * 100
        if bad_pct > 5.0:
            var_name = col.replace('_Range_QC', '')
            issues.append(f"  • {var_name}: {bad_pct:.1f}% BAD flags - Expert review recommended")
    
    if issues:
        f.write("⚠️  ISSUES DETECTED:\n\n")
        for issue in issues:
            f.write(issue + "\n")
    else:
        f.write("✓ No significant quality issues detected.\n")
    
    f.write("\n" + "=" * 70 + "\n")

print(f"✓ Report saved: {report_file}")

# ============================================================================
# COMPLETE
# ============================================================================

print("\n" + "=" * 70)
print("QC ANALYSIS COMPLETE!")
print("=" * 70)

print("\nGenerated files:")
print(f"  • {strict_file}")
print(f"  • {moderate_file}")
print(f"  • {liberal_file}")
print(f"  • {output_file}")
print(f"  • {report_file}")

print("\nNext steps:")
print("  1. Review QC_Report.txt for detailed statistics")
print("  2. Examine qc_timeseries_overview.png for visual patterns")
print("  3. Use filtered CSV files for your analysis")
print("\nDocumentation: https://github.com/BennyTorelli/Seaexplorer_pyglider/tree/main/documentation")
