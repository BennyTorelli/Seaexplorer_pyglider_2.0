#!/usr/bin/env python3
"""
Visualize test points and La Palma polygon on a map.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import numpy as np

# Load La Palma polygon
polygon_file = 'config/shapefiles/LaPalmaDissolve_coords.json'
with open(polygon_file, 'r') as f:
    coords_list = json.load(f)

polygon_coords = coords_list[0]
lapalma_polygon = Polygon(polygon_coords)

# Extract polygon boundary for plotting
poly_lons = [coord[0] for coord in polygon_coords]
poly_lats = [coord[1] for coord in polygon_coords]

# Load test results
test_results = pd.read_csv('tests/data/test_land_qc_results.csv')
expected = pd.read_csv('tests/data/expected_land_qc.csv')

# Merge to get descriptions
test_results = test_results.merge(
    expected[['latitude', 'longitude', 'description', 'expected_land_qc']], 
    left_on=['LATITUDE', 'LONGITUDE'], 
    right_on=['latitude', 'longitude'],
    how='left'
)

# Separate points by result
correct_ocean = test_results[(test_results['LAND_QC'] == 1) & (test_results['expected_land_qc'] == 1)]
correct_land = test_results[(test_results['LAND_QC'] == 4) & (test_results['expected_land_qc'] == 4)]
incorrect = test_results[test_results['LAND_QC'] != test_results['expected_land_qc']]

# Create figure
fig, ax = plt.subplots(figsize=(12, 14))

# Plot La Palma polygon
ax.plot(poly_lons, poly_lats, 'k-', linewidth=2, label='La Palma coastline', zorder=1)
ax.fill(poly_lons, poly_lats, color='lightgray', alpha=0.5, zorder=0)

# Plot test points
if len(correct_ocean) > 0:
    ax.scatter(correct_ocean['LONGITUDE'], correct_ocean['LATITUDE'], 
              c='blue', s=200, marker='o', edgecolors='black', linewidth=2,
              label='✓ Correct - Ocean (flag 1)', zorder=3)

if len(correct_land) > 0:
    ax.scatter(correct_land['LONGITUDE'], correct_land['LATITUDE'], 
              c='red', s=200, marker='X', edgecolors='black', linewidth=2,
              label='✓ Correct - Land (flag 4)', zorder=3)

if len(incorrect) > 0:
    ax.scatter(incorrect['LONGITUDE'], incorrect['LATITUDE'], 
              c='yellow', s=300, marker='*', edgecolors='red', linewidth=3,
              label='✗ Incorrect', zorder=4)

# Add labels for each point
for idx, row in test_results.iterrows():
    # Offset text slightly to avoid overlap
    offset_x = 0.01
    offset_y = 0.01
    
    # Create label with point number and flag
    point_num = idx + 1
    actual_flag = int(row['LAND_QC']) if not pd.isna(row['LAND_QC']) else '?'
    expected_flag = int(row['expected_land_qc'])
    
    status = '✓' if actual_flag == expected_flag else '✗'
    label = f"{point_num}. {status}"
    
    ax.annotate(label, 
               xy=(row['LONGITUDE'], row['LATITUDE']),
               xytext=(row['LONGITUDE'] + offset_x, row['LATITUDE'] + offset_y),
               fontsize=10, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.8),
               zorder=5)

# Add description text box
textstr = "Test Points:\n" + "-" * 40 + "\n"
for idx, row in test_results.iterrows():
    point_num = idx + 1
    desc = row['description']
    actual = int(row['LAND_QC']) if not pd.isna(row['LAND_QC']) else '?'
    expected = int(row['expected_land_qc'])
    status = '✓' if actual == expected else '✗'
    textstr += f"{point_num}. {status} {desc}\n"
    textstr += f"   Expected: {expected} | Actual: {actual}\n"

props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', bbox=props, family='monospace')

# Labels and formatting
ax.set_xlabel('Longitude (°E)', fontsize=12, fontweight='bold')
ax.set_ylabel('Latitude (°N)', fontsize=12, fontweight='bold')
ax.set_title('LAND_QC Test Results - La Palma', fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='lower right', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_aspect('equal')

# Add coordinate bounds text
bounds_text = f"La Palma Bounds:\nLon: {lapalma_polygon.bounds[0]:.4f} to {lapalma_polygon.bounds[2]:.4f}\nLat: {lapalma_polygon.bounds[1]:.4f} to {lapalma_polygon.bounds[3]:.4f}"
ax.text(0.98, 0.02, bounds_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

plt.tight_layout()

# Save figure
output_file = 'tests/data/land_qc_test_map.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ Map saved to: {output_file}")

# Also create a summary statistics
print(f"\n{'='*60}")
print(f"LAND_QC TEST SUMMARY")
print(f"{'='*60}")
print(f"Total test points: {len(test_results)}")
print(f"Correct ocean (flag 1): {len(correct_ocean)}")
print(f"Correct land (flag 4): {len(correct_land)}")
print(f"Incorrect: {len(incorrect)}")
print(f"Accuracy: {(len(correct_ocean) + len(correct_land)) / len(test_results) * 100:.1f}%")
print(f"{'='*60}\n")

if len(incorrect) > 0:
    print("Failed points:")
    for idx, row in incorrect.iterrows():
        print(f"  - {row['description']}")
        print(f"    LAT: {row['LATITUDE']:.6f}, LON: {row['LONGITUDE']:.6f}")
        print(f"    Expected: {int(row['expected_land_qc'])}, Got: {int(row['LAND_QC'])}")

plt.show()
