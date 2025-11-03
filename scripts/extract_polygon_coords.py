#!/usr/bin/env python3
"""
Extract polygon coordinates from La Palma shapefile.
This script reads the shapefile and extracts the coordinates of the polygon vertices.
"""

import geopandas as gpd
import json

def extract_polygon_coordinates(shapefile_path):
    """
    Extract coordinates from a shapefile polygon.
    
    Parameters:
    -----------
    shapefile_path : str
        Path to the shapefile
    
    Returns:
    --------
    list : List of (lon, lat) tuples representing polygon vertices
    """
    print(f"\n{'='*60}")
    print(f"Reading shapefile: {shapefile_path}")
    print(f"{'='*60}\n")
    
    # Read shapefile
    gdf = gpd.read_file(shapefile_path)
    
    # Print basic info
    print(f"Number of features: {len(gdf)}")
    print(f"CRS: {gdf.crs}")
    print(f"Geometry types: {gdf.geometry.type.unique()}")
    print(f"\nBounds:")
    print(f"  Min Lon: {gdf.total_bounds[0]:.6f}")
    print(f"  Min Lat: {gdf.total_bounds[1]:.6f}")
    print(f"  Max Lon: {gdf.total_bounds[2]:.6f}")
    print(f"  Max Lat: {gdf.total_bounds[3]:.6f}")
    
    # Extract coordinates from first geometry
    geom = gdf.geometry.iloc[0]
    
    # Handle different geometry types
    coords_list = []
    
    if geom.geom_type == 'Polygon':
        # Simple polygon - get exterior coordinates
        coords = list(geom.exterior.coords)
        coords_list.append(coords)
        print(f"\nGeometry type: Polygon")
        print(f"Number of vertices: {len(coords)}")
        
    elif geom.geom_type == 'MultiPolygon':
        # Multiple polygons - get all exterior coordinates
        print(f"\nGeometry type: MultiPolygon")
        print(f"Number of polygons: {len(geom.geoms)}")
        for i, poly in enumerate(geom.geoms):
            coords = list(poly.exterior.coords)
            coords_list.append(coords)
            print(f"  Polygon {i+1}: {len(coords)} vertices")
    
    # Print first few coordinates as example
    print(f"\nFirst 5 coordinates (lon, lat):")
    for i, coord in enumerate(coords_list[0][:5]):
        print(f"  {i+1}. ({coord[0]:.6f}, {coord[1]:.6f})")
    
    if len(coords_list[0]) > 5:
        print(f"  ... ({len(coords_list[0])-5} more vertices)")
    
    return coords_list


if __name__ == "__main__":
    import sys
    
    # Try both shapefiles
    shapefiles = [
        "config/shapefiles/LaPalmaDissolve.shp",
        "config/shapefiles/La_Palma.shp"
    ]
    
    for shp_file in shapefiles:
        try:
            coords_list = extract_polygon_coordinates(shp_file)
            
            # Save to JSON file
            output_file = shp_file.replace('.shp', '_coords.json')
            with open(output_file, 'w') as f:
                json.dump(coords_list, f, indent=2)
            print(f"\n✓ Coordinates saved to: {output_file}")
            
        except Exception as e:
            print(f"\n✗ Error reading {shp_file}: {e}")
        
        print("\n")
