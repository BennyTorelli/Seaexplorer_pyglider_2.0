#!/usr/bin/env python3
"""
Visualize LAND_QC test results on a map.
Shows the La Palma polygon and test points colored by their LAND_QC flag.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPLPolygon
import numpy as np

def plot_test_results():
    """Create a map visualization of test results."""
    
    # Load La Palma polygon
    with open('config/shapefiles/LaPalmaDissolve_coords.json', 'r') as f:
        coords_list = json.load(f)
    polygon_coords = coords_list[0]
    
    # Load test results
    results = pd.read_csv('tests/data/test_land_qc_results.csv')
    expected = pd.read_csv('tests/data/expected_land_qc.csv')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 14))
    
    # Plot La Palma island polygon
    lons = [coord[0] for coord in polygon_coords]
    lats = [coord[1] for coord in polygon_coords]
    island_polygon = MPLPolygon(list(zip(lons, lats)), 
                                fill=True, 
                                facecolor='lightgray', 
                                edgecolor='black', 
                                linewidth=2,
                                alpha=0.5,
                                label='La Palma Island')
    ax.add_patch(island_polygon)
    
    # Plot test points
    for i in range(len(results)):
        lat = results['LATITUDE'].iloc[i]
        lon = results['LONGITUDE'].iloc[i]
        land_qc = results['LAND_QC'].iloc[i]
        expected_qc = expected['expected_land_qc'].iloc[i]
        desc = expected['description'].iloc[i]
        
        # Determine color and marker
        if land_qc == expected_qc:
            # Correct prediction
            if land_qc == 1:
                color = 'blue'
                marker = 'o'
                label_text = 'Ocean (Correct)'
            else:
                color = 'red'
                marker = 's'
                label_text = 'Land (Correct)'
        else:
            # Wrong prediction
            color = 'orange'
            marker = 'X'
            label_text = 'Mismatch'
        
        ax.scatter(lon, lat, c=color, marker=marker, s=200, 
                  edgecolors='black', linewidths=2, zorder=10)
        
        # Add point number
        ax.text(lon, lat, str(i+1), fontsize=10, fontweight='bold',
               ha='center', va='center', color='white', zorder=11)
    
    # Create custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Ocean (Flag 1 - Correct)',
               markerfacecolor='blue', markersize=12, markeredgecolor='black', markeredgewidth=2),
        Line2D([0], [0], marker='s', color='w', label='Land (Flag 4 - Correct)',
               markerfacecolor='red', markersize=12, markeredgecolor='black', markeredgewidth=2),
        Line2D([0], [0], marker='X', color='w', label='Mismatch',
               markerfacecolor='orange', markersize=12, markeredgecolor='black', markeredgewidth=2),
        MPLPolygon([(0,0)], facecolor='lightgray', edgecolor='black', 
                   linewidth=2, alpha=0.5, label='La Palma Island')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)
    
    # Set axis labels and title
    ax.set_xlabel('Longitude (°E)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Latitude (°N)', fontsize=12, fontweight='bold')
    ax.set_title('LAND_QC Test Results - La Palma Island\nTest Points with Expected vs Actual Flags', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Set axis limits with some padding
    ax.set_xlim(-18.15, -17.65)
    ax.set_ylim(28.35, 28.95)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_aspect('equal')
    
    # Add text box with statistics
    n_correct = sum(results['LAND_QC'] == expected['expected_land_qc'])
    n_total = len(results)
    textstr = f'Results: {n_correct}/{n_total} correct ({n_correct/n_total*100:.0f}%)'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props, fontweight='bold')
    
    # Save figure
    output_file = 'tests/data/land_qc_map.png'
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Map saved to: {output_file}")
    
    # Also create a detailed legend with point descriptions
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    ax2.axis('off')
    
    # Create table with results
    table_data = []
    for i in range(len(results)):
        lat = results['LATITUDE'].iloc[i]
        lon = results['LONGITUDE'].iloc[i]
        land_qc = results['LAND_QC'].iloc[i]
        expected_qc = expected['expected_land_qc'].iloc[i]
        desc = expected['description'].iloc[i]
        status = '✓' if land_qc == expected_qc else '✗'
        
        table_data.append([
            f"{i+1}",
            desc,
            f"{lat:.4f}",
            f"{lon:.4f}",
            str(expected_qc),
            str(land_qc),
            status
        ])
    
    table = ax2.table(cellText=table_data,
                     colLabels=['#', 'Description', 'Lat', 'Lon', 'Expected', 'Actual', 'Status'],
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.05, 0.4, 0.1, 0.1, 0.1, 0.1, 0.08])
    
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # Color cells based on status
    for i in range(len(table_data)):
        if table_data[i][6] == '✓':
            table[(i+1, 6)].set_facecolor('#90EE90')  # Light green
        else:
            table[(i+1, 6)].set_facecolor('#FFB6C1')  # Light red
    
    # Header styling
    for j in range(7):
        table[(0, j)].set_facecolor('#4472C4')
        table[(0, j)].set_text_props(weight='bold', color='white')
    
    ax2.set_title('LAND_QC Test Results - Detailed Table', 
                 fontsize=14, fontweight='bold', pad=20)
    
    output_file2 = 'tests/data/land_qc_table.png'
    plt.tight_layout()
    plt.savefig(output_file2, dpi=300, bbox_inches='tight')
    print(f"✓ Table saved to: {output_file2}")
    
    print("\nVisualization complete!")
    print(f"  - Map: {output_file}")
    print(f"  - Table: {output_file2}")


if __name__ == '__main__':
    plot_test_results()
