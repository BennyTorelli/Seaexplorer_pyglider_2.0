# Technical Manual: SeaExplorer Data Processing Pipeline

**Repository:** https://github.com/BennyTorelli/Seaexplorer_pyglider  
**Version:** 1.0  
**Last Updated:** October 2025

---

## Table of Contents

1. [System Requirements](#1-system-requirements)
2. [Installation](#2-installation)
3. [Configuration](#3-configuration)
4. [Pipeline Execution](#4-pipeline-execution)
5. [Quality Control](#5-quality-control)
6. [Output Products](#6-output-products)
7. [Troubleshooting](#7-troubleshooting)
8. [Advanced Usage](#8-advanced-usage)

---

## 1. System Requirements

### 1.1 Hardware
- **Minimum:**
  - CPU: Dual-core processor
  - RAM: 8 GB
  - Storage: 5 GB free space per mission (raw + processed data)
  
- **Recommended:**
  - CPU: Quad-core processor (M1/M2 or Intel i5/i7)
  - RAM: 16 GB or more
  - Storage: SSD with 10+ GB free space

### 1.2 Software
- **Operating System:**
  - macOS 13.0+ (Ventura or later)
  - Linux (Ubuntu 20.04+, CentOS 8+)
  - Windows 10/11 (with WSL2 recommended)

- **Python:** Version 3.11+ (3.13.5 recommended)

- **Git:** For repository cloning and version control

### 1.3 Dependencies
See `requirements.txt` for complete list. Key packages:
```
pyglider==0.0.7
xarray==2025.9.0
numpy==2.3.3
pandas==2.3.2
gsw==3.6.20
shapely==2.1.2
geopandas==1.1.1
netCDF4==1.7.2
matplotlib==3.9.0
```

---

## 2. Installation

### 2.1 Clone Repository
```bash
git clone https://github.com/BennyTorelli/Seaexplorer_pyglider.git
cd Seaexplorer_pyglider
```

### 2.2 Create Virtual Environment
```bash
# Create environment
python3 -m venv .venv

# Activate environment
# macOS/Linux:
source .venv/bin/activate
# Windows:
.venv\Scripts\activate
```

### 2.3 Install Dependencies
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### 2.4 Verify Installation
```bash
# Test PyGlider import
python -c "import pyglider; print(f'PyGlider version: {pyglider.__version__}')"

# Test GSW import
python -c "import gsw; print(f'GSW version: {gsw.__version__}')"

# Test Shapely import
python -c "import shapely; print(f'Shapely version: {shapely.__version__}')"
```

Expected output:
```
PyGlider version: 0.0.7
GSW version: 3.6.20
Shapely version: 2.1.2
```

---

## 3. Configuration

### 3.1 Directory Structure
The repository expects the following structure:
```
Seaexplorer_pyglider/
├── config/
│   ├── seaexplorer_0067.yml       # Deployment configuration
│   └── shapefiles/                 # Geographic boundaries
│       ├── LaPalmaDissolve.shp
│       ├── LaPalmaDissolve_coords.json
│       └── ...
├── input/
│   └── raw/                        # Raw glider files (.pld1, .gli, .gz)
├── output/
│   ├── l0_data/
│   │   ├── rawnc/                 # Intermediate parquet files
│   │   └── timeseries/            # L0 netCDF files
│   └── analysis/                   # Final products (CSV, standardized NC)
├── scripts/
│   ├── 01_setup_pyglider.py
│   ├── qc_variables.py
│   └── extract_polygon_coords.py
├── tests/                          # Test suite
├── MASTER_pyglider_pipeline.py     # Main processing script
└── requirements.txt
```

### 3.2 Deployment Configuration File

Edit `config/seaexplorer_0067.yml` for your mission:

```yaml
# Basic mission metadata
metadata:
  acknowledgement: 'Funding source and deployment team'
  comment: 'Mission-specific notes'
  contributor_name: 'Your Name'
  contributor_role: 'Data processor'
  creator_email: 'your.email@institution.org'
  creator_name: 'Your Name'
  creator_url: 'https://your-institution.org'
  institution: 'Your Institution'
  project: 'Project Name'
  sea_name: 'Eastern North Atlantic'
  
# Glider identifica information
glider_model: seaexplorer
glider_name: sea074
glider_serial: '074'
wmo_id: '1234567'  # If assigned

# Mission timing
deployment_start: '2025-02-12'
deployment_end: '2025-03-05'

# Geographic boundaries (for initial checks)
lon1: -18.5  # West
lon2: -17.5  # East
lat1: 28.0   # South
lat2: 29.0   # North

# Sensor configuration
# Temperature
temperature:
  sensor_name: 'SBE 41CP CTD'
  long_name: 'Temperature'
  standard_name: 'sea_water_temperature'
  units: 'Celsius'
  comment: 'CTD temperature, ITS-90 scale'

# Conductivity
conductivity:
  sensor_name: 'SBE 41CP CTD'
  long_name: 'Conductivity'
  standard_name: 'sea_water_electrical_conductivity'
  units: 'S m-1'
  comment: 'CTD conductivity, calibrated'

# Pressure
pressure:
  sensor_name: 'SBE 41CP CTD'
  long_name: 'Pressure'
  standard_name: 'sea_water_pressure'
  units: 'dbar'
  comment: 'CTD pressure, absolute'

# Chlorophyll
chlorophyll:
  sensor_name: 'WET Labs ECO Puck'
  long_name: 'Chlorophyll-a concentration'
  standard_name: 'concentration_of_chlorophyll_in_sea_water'
  units: 'mg m-3'
  comment: 'Fluorescence-derived chlorophyll-a'
  dark_count: 50  # Factory calibration
  scale_factor: 0.0121  # Factory calibration

# CDOM
cdom:
  sensor_name: 'WET Labs ECO Puck'
  long_name: 'Colored Dissolved Organic Matter'
  units: 'ppb'
  comment: 'CDOM fluorescence, 370ex/460em'
  dark_count: 45
  scale_factor: 0.091

# Backscatter (turbidity)
backscatter_700:
  sensor_name: 'WET Labs ECO Puck'
  long_name: 'Optical backscatter at 700nm'
  standard_name: 'volume_backwards_scattering_coefficient_of_radiative_flux_in_sea_water'
  units: 'NTU'
  comment: 'Converted from backscatter using scale factor 1/0.002727'
  dark_count: 48
  scale_factor: 0.002727  # Used in conversion (NTU = β / 0.002727)

# Dissolved oxygen
oxygen_concentration:
  sensor_name: 'Aanderaa Optode 4330'
  long_name: 'Dissolved oxygen concentration'
  standard_name: 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
  units: 'umol kg-1'
  comment: 'Mass-based units, converted from volume-based using density'
  coarsen: 8  # Averaging window
  correct_oxygen: 1  # Apply lag correction

# Salinity (derived)
salinity:
  long_name: 'Practical Salinity'
  standard_name: 'sea_water_practical_salinity'
  units: '1'
  comment: 'Calculated using GSW TEOS-10, recalculated after conductivity conversion'

# Density (derived)
density:
  long_name: 'Sea Water Density'
  standard_name: 'sea_water_density'
  units: 'kg m-3'
  comment: 'Calculated using GSW TEOS-10 equation of state'
```

### 3.3 QC Configuration

QC parameters are defined in `scripts/qc_variables.py` (lines 12-22):

```python
QC_RANGES_LAPALMA = {
    'TEMP': {'min': 5.0, 'max': 30.0, 'unit': '°C'},
    'PSAL': {'min': 33.0, 'max': 38.0, 'unit': 'PSU'},
    'CNDC': {'min': 0.0, 'max': 7.0, 'unit': 'S/m'},
    'DOXY': {'min': 110.0, 'max': 250.0, 'unit': 'µmol/L'},
    'CHLA': {'min': 0.0, 'max': 42.0, 'unit': 'mg/m³'},
    'TURB': {'min': 0.0, 'max': 5.0, 'unit': 'NTU'},
    'PRES': {'min': 0.0, 'max': 2000.0, 'unit': 'dbar'},
}
```

**To adapt for your region:**
1. Consult climatological atlases (World Ocean Atlas, regional climatologies)
2. Review historical data from the deployment area
3. Adjust min/max values conservatively (capture 99% of expected values)

### 3.4 Geographic Boundary Configuration

For LAND_QC, you need a coastline shapefile:

**Option 1: Use existing (La Palma example)**
- File: `config/shapefiles/LaPalmaDissolve.shp`
- Coordinates: Already extracted to JSON
- No action needed if deploying near La Palma

**Option 2: Create new for different region**
```bash
# 1. Obtain coastline shapefile for your region
#    Sources: Natural Earth, GSHHG, national hydrographic offices

# 2. Place shapefile in config/shapefiles/
#    Required files: .shp, .shx, .dbf, .prj

# 3. Extract coordinates to JSON
python scripts/extract_polygon_coords.py

# 4. Update qc_variables.py line 595 with new JSON filename
polygon_file = 'config/shapefiles/YOUR_REGION_coords.json'
```

---

## 4. Pipeline Execution

### 4.1 Prepare Raw Data

```bash
# 1. Copy raw files to input/raw/
cp /path/to/mission/data/*.pld1 input/raw/
cp /path/to/mission/data/*.gli input/raw/
cp /path/to/mission/data/*.pld0 input/raw/  # If available

# 2. If files are compressed, they will be auto-decompressed by pipeline
# Supported formats: .gz

# 3. Verify files are present
ls -lh input/raw/
```

Expected output:
```
sea074.0000.pld1.gz
sea074.0001.pld1.gz
...
sea074.0000.gli.gz
sea074.0001.gli.gz
...
```

### 4.2 Run Complete Pipeline

**Method 1: All-in-one execution**
```bash
python MASTER_pyglider_pipeline.py
```

Expected console output:
```
Inizializzazione processore SeaExplorer...
PYGLIDER SEAEXPLORER - PIPELINE COMPLETA
===============================================

STEP 0: Verifica e decompressione file .gz
------------------------------------------
  Decompressione di 555 file .gz...
  Completato: 555 file decompressi

STEP 1: Conversione raw → parquet
-----------------------------------
   Processando 555 file totali...
Step 1 completato: 549 file parquet creati

STEP 2: Merge file parquet
---------------------------
Step 2 completato: file parquet uniti

STEP 3: Creazione L0 timeseries
---------------------------------
Step 3 completato: sea074-2025.nc (306.5 MB)

STEP 3b: Applicazione conversioni personalizzate
-----------------------------------------------
   Conversione TURB: backscatter → NTU
   Conversione CNDC: mS/cm → S/m
   Ricalcolo SALINITÀ: median=35.97 PSU
   Ricalcolo DENSITÀ: median=1028.7 kg/m³
   Conversione DOXY: → µmol/kg
   NetCDF aggiornato con conversioni

STEP 4: Conversione NetCDF → CSV
---------------------------------
Step 4 completato: CSV creato (565.0 MB)

STEP 5: Analisi qualità dati
------------------------------
Analisi L0:
   Campioni: 1,825,662
   Salinita: 0.00-37.22 PSU
   Temperatura: 7.53-28.14 °C

STEP 6: Rinomina variabili con nomi standard
---------------------------------------------
   L0 rinominato: sea074-2025_standard_names.nc
   CSV rinominato: seaexplorer_data_complete_standard_names.csv

PIPELINE COMPLETATA CON SUCCESSO!
```

**Method 2: Run individual steps** (for debugging or selective processing)
```python
from MASTER_pyglider_pipeline import SeaExplorerProcessor

processor = SeaExplorerProcessor()

# Step 1 only
processor.step1_raw_to_parquet()

# Step 2 only
processor.step2_merge_parquet()

# ... and so on
```

### 4.3 Run Quality Control

```bash
# Process latest L0 netCDF file
python scripts/qc_variables.py --latest

# Or specify file
python scripts/qc_variables.py output/analysis/sea074-2025_standard_names.nc
```

Expected output:
```
Wrote QC CSV to output/analysis/seaexplorer_qc_variables.csv
```

### 4.4 Verify Outputs

```bash
# Check file sizes
ls -lh output/analysis/

# Expected files:
# sea074-2025_standard_names.nc          ~300-350 MB
# seaexplorer_data_complete_standard_names.csv  ~500-600 MB
# seaexplorer_qc_variables.csv           ~300-400 MB

# Preview NetCDF metadata
ncdump -h output/analysis/sea074-2025_standard_names.nc | head -50

# Preview CSV headers
head -1 output/analysis/seaexplorer_data_complete_standard_names.csv

# Check QC column count
head -1 output/analysis/seaexplorer_qc_variables.csv | tr ',' '\n' | wc -l
# Expected: 30 columns
```

---

## 5. Quality Control

### 5.1 Interpreting QC Flags

**QC Flag Values:**
- `1` = GOOD: Value passes all quality checks
- `4` = BAD: Value fails quality check (out of range, suspect)
- `9` = MISSING: Value not available (NaN, null)

**QC Column Naming Convention:**
```
{VARIABLE}_{QC_TYPE}_QC

Examples:
- TEMP_Range_QC: Temperature oceanographic range check
- TEMP_Sensor_QC: Temperature sensor physical limit check
- TEMP_Na_QC: Temperature missing value check
```

### 5.2 QC CSV Structure

The QC CSV (`seaexplorer_qc_variables.csv`) contains 30 columns:

**Columns 1-11: Data**
1. TIME
2. DEPTH
3. LATITUDE
4. LONGITUDE
5. TEMP
6. CNDC
7. PRES
8. CHLA
9. TURB
10. DOXY
11. PSAL

**Columns 12-19: Range QC (oceanographic)**
12. Date_QC
13. Location_QC
14. TEMP_Range_QC
15. CNDC_Range_QC
16. CHLA_Range_QC
17. TURB_Range_QC
18. DOXY_Range_QC
19. PSAL_Range_QC

**Columns 20-24: Missing Value QC**
20. TEMP_Na_QC
21. CNDC_Na_QC
22. DOXY_Na_QC
23. CHLA_Na_QC
24. TURB_Na_QC

**Columns 25-29: Sensor QC (physical limits)**
25. TEMP_Sensor_QC
26. CNDC_Sensor_QC
27. DOXY_Sensor_QC
28. CHLA_Sensor_QC
29. TURB_Sensor_QC

**Column 30: Geographic QC**
30. LAND_QC

### 5.3 Common QC Patterns

**All flags = 1 (perfect data):**
```csv
TIME,DEPTH,...,TEMP_Range_QC,TEMP_Sensor_QC,TEMP_Na_QC,...
2025-02-12 10:43:37,0.5,...,1,1,1,...
```
✓ No action needed

**Missing value (Na_QC = 9):**
```csv
TIME,DEPTH,TEMP,...,TEMP_Na_QC
2025-02-12 10:43:37,0.5,NaN,...,9
```
⚠️ Sensor malfunction or data gap

**Out of range (Range_QC = 4):**
```csv
TIME,DEPTH,TEMP,TEMP_Range_QC
2025-02-12 10:43:37,0.5,32.5,4
```
⚠️ Unexpected value (check for calibration drift, biofouling, or real anomaly)

**Sensor limit exceeded (Sensor_QC = 4):**
```csv
TIME,DEPTH,TEMP,TEMP_Sensor_QC
2025-02-12 10:43:37,0.5,45.0,4
```
🚨 Sensor failure (physically impossible value)

**On land (LAND_QC = 4):**
```csv
TIME,LATITUDE,LONGITUDE,LAND_QC
2025-02-12 10:43:37,28.65,-17.85,4
```
🚨 GPS error or navigation failure (glider cannot be on land!)

### 5.4 QC Workflow Recommendations

**Step 1: Overview Statistics**
```bash
# Count QC flags by type
python -c "
import pandas as pd
qc = pd.read_csv('output/analysis/seaexplorer_qc_variables.csv')
qc_cols = [c for c in qc.columns if '_QC' in c]
for col in qc_cols:
    print(f'{col}:')
    print(qc[col].value_counts().sort_index())
    print()
"
```

**Step 2: Identify Problem Periods**
```python
import pandas as pd
import matplotlib.pyplot as plt

qc = pd.read_csv('output/analysis/seaexplorer_qc_variables.csv', parse_dates=['TIME'])

# Plot QC flags over time
fig, axes = plt.subplots(5, 1, figsize=(12, 10), sharex=True)
variables = ['TEMP', 'CNDC', 'DOXY', 'CHLA', 'TURB']

for ax, var in zip(axes, variables):
    qc_col = f'{var}_Range_QC'
    bad = qc[qc[qc_col] == 4]
    ax.scatter(bad['TIME'], bad['DEPTH'], c='red', s=1, label='Flag 4 (Bad)')
    ax.set_ylabel(f'{var} Depth (m)')
    ax.invert_yaxis()
    ax.legend()

axes[-1].set_xlabel('Time')
plt.tight_layout()
plt.savefig('output/analysis/qc_overview.png', dpi=300)
```

**Step 3: Expert Review**
For data with significant flagging (>5%), manual inspection is recommended:
1. Plot flagged data in context (time series, T-S diagrams)
2. Compare with climatology or nearby observations
3. Check sensor calibration history
4. Consult with deployment team about known issues

**Step 4: Data Filtering** 
```python
# Create "good data only" subset
qc = pd.read_csv('output/analysis/seaexplorer_qc_variables.csv')

# Keep only data where all Range_QC flags = 1
range_qc_cols = [c for c in qc.columns if 'Range_QC' in c]
good_data = qc[(qc[range_qc_cols] == 1).all(axis=1)]

# Save filtered dataset
good_data.to_csv('output/analysis/seaexplorer_data_qc_passed.csv', index=False)

print(f'Original: {len(qc)} points')
print(f'QC passed: {len(good_data)} points ({len(good_data)/len(qc)*100:.1f}%)')
```

---

## 6. Output Products

### 6.1 NetCDF Files

**File:** `output/analysis/sea074-2025_standard_names.nc`

**Structure:**
```
dimensions:
    time = 1825662  # Number of measurements
    
variables:
    TIME(time): datetime64[ns]
    DEPTH(time): float64 [m]
    LATITUDE(time): float64 [degrees_north]
    LONGITUDE(time): float64 [degrees_east]
    TEMP(time): float32 [Celsius]
    CNDC(time): float32 [S/m]
    PRES(time): float32 [dbar]
    CHLA(time): float32 [mg/m³]
    CDOM(time): float32 [ppb]
    TURB(time): float32 [NTU]
    DOXY(time): float32 [µmol/kg]
    PSAL(time): float32 [PSU]
    density(time): float32 [kg/m³]
    ... (auxiliary variables)

global attributes:
    Conventions: CF-1.8
    title: SeaExplorer glider data - Mission sea074
    institution: [Your Institution]
    source: SeaExplorer autonomous underwater glider
    history: [Processing history timestamp]
    references: https://github.com/BennyTorelli/Seaexplorer_pyglider
    comment: Processed with PyGlider, TEOS-10 conversions applied
```

**Reading with Python:**
```python
import xarray as xr

ds = xr.open_dataset('output/analysis/sea074-2025_standard_names.nc')

# Access variables
temperature = ds['TEMP'].values
time = ds['TIME'].values
lat = ds['LATITUDE'].values

# Subset by time
ds_feb = ds.sel(time=slice('2025-02-01', '2025-03-01'))

# Subset by depth
shallow = ds.where(ds['DEPTH'] < 100, drop=True)
```

**Reading with MATLAB:**
```matlab
filename = 'output/analysis/sea074-2025_standard_names.nc';

% Read variables
temp = ncread(filename, 'TEMP');
time = ncread(filename, 'TIME');
depth = ncread(filename, 'DEPTH');

% Read attributes
units = ncreadatt(filename, 'TEMP', 'units');
```

### 6.2 CSV Files

**File 1:** `seaexplorer_data_complete_standard_names.csv`

- **Purpose:** Full dataset with all variables
- **Size:** ~565 MB (1.8M rows × 22 columns)
- **Use cases:** 
  - GIS software (QGIS, ArcGIS)
  - Statistical analysis (R, SPSS)
  - Spreadsheet import (Excel - note: may truncate at ~1M rows)

**File 2:** `seaexplorer_qc_variables.csv`

- **Purpose:** QC flags and selected variables
- **Size:** ~348 MB (1.8M rows × 30 columns)
- **Use cases:**
  - Quality control review
  - Data filtering
  - Repository submission

**CSV Import Best Practices:**

```python
# Python/Pandas (efficient for large files)
import pandas as pd

# Read in chunks to manage memory
chunksize = 100000
chunks = []
for chunk in pd.read_csv('output/analysis/seaexplorer_data_complete_standard_names.csv', 
                         chunksize=chunksize, 
                         parse_dates=['TIME']):
    chunks.append(chunk)
df = pd.concat(chunks, ignore_index=True)

# Or read specific columns only
df = pd.read_csv('output/analysis/seaexplorer_data_complete_standard_names.csv',
                 usecols=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE', 'TEMP', 'PSAL'])
```

```R
# R (readr package for large files)
library(readr)

data <- read_csv("output/analysis/seaexplorer_data_complete_standard_names.csv",
                 col_types = cols(TIME = col_datetime()))

# Or use data.table for very large files
library(data.table)
data <- fread("output/analysis/seaexplorer_data_complete_standard_names.csv")
```

### 6.3 Visualization Examples

**Quick Profile Plot:**
```python
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('output/analysis/seaexplorer_data_complete_standard_names.csv',
                   usecols=['DEPTH', 'TEMP', 'PSAL', 'DOXY'])

# Select first profile (depths increasing then decreasing)
profile = data.iloc[:500]  # Adjust based on your data

fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharey=True)

axes[0].plot(profile['TEMP'], profile['DEPTH'], 'b-')
axes[0].set_xlabel('Temperature (°C)')
axes[0].set_ylabel('Depth (m)')
axes[0].invert_yaxis()
axes[0].grid(True)

axes[1].plot(profile['PSAL'], profile['DEPTH'], 'g-')
axes[1].set_xlabel('Salinity (PSU)')
axes[1].grid(True)

axes[2].plot(profile['DOXY'], profile['DEPTH'], 'r-')
axes[2].set_xlabel('Oxygen (µmol/kg)')
axes[2].grid(True)

plt.tight_layout()
plt.savefig('output/analysis/example_profile.png', dpi=300)
```

**Trajectory Map:**
```python
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

data = pd.read_csv('output/analysis/seaexplorer_data_complete_standard_names.csv',
                   usecols=['LATITUDE', 'LONGITUDE'])

# Subsample for faster plotting
data_sub = data.iloc[::100]  # Every 100th point

plt.figure(figsize=(10, 8))
plt.scatter(data_sub['LONGITUDE'], data_sub['LATITUDE'], 
           c=range(len(data_sub)), cmap='viridis', s=1)
plt.colorbar(label='Time (relative)')
plt.xlabel('Longitude (°E)')
plt.ylabel('Latitude (°N)')
plt.title('Glider Trajectory')
plt.grid(True)
plt.axis('equal')
plt.savefig('output/analysis/trajectory.png', dpi=300)
```

---

## 7. Troubleshooting

### 7.1 Common Errors

#### Error: "ModuleNotFoundError: No module named 'pyglider'"
**Cause:** PyGlider not installed or virtual environment not activated

**Solution:**
```bash
# Ensure virtual environment is activated
source .venv/bin/activate  # macOS/Linux
# or
.venv\Scripts\activate  # Windows

# Install PyGlider
pip install pyglider
```

#### Error: "FileNotFoundError: [Errno 2] No such file or directory: 'input/raw/'"
**Cause:** Directory structure not set up

**Solution:**
```bash
# Create required directories
mkdir -p input/raw
mkdir -p output/l0_data/rawnc
mkdir -p output/l0_data/timeseries
mkdir -p output/analysis
```

#### Error: "ValueError: No objects to concatenate" (Step 2)
**Cause:** No parquet files found from Step 1

**Solution:**
```bash
# Check if Step 1 completed successfully
ls -l output/l0_data/rawnc/*.parquet

# If empty, re-run Step 1
python MASTER_pyglider_pipeline.py  # It will auto-detect and run from Step 1
```

#### Error: "KeyError: 'temperature'" (Step 3b)
**Cause:** Variable name mismatch in configuration

**Solution:**
Check `config/seaexplorer_0067.yml` for correct sensor names. Ensure:
```yaml
temperature:  # Lowercase, matches PyGlider convention
  sensor_name: 'SBE 41CP'
  ...
```

#### Warning: "LAND_QC: La Palma polygon file not found"
**Cause:** Shapefile coordinates not extracted

**Solution:**
```bash
# Extract coordinates from shapefile
python scripts/extract_polygon_coords.py

# Or copy existing JSON
cp config/shapefiles/LaPalmaDissolve_coords.json config/shapefiles/
```

### 7.2 Performance Issues

#### Slow Step 1 (>15 minutes for 555 files)
**Possible causes:**
- Slow disk I/O (e.g., network drive)
- Insufficient RAM (swapping to disk)
- Large number of small files

**Solutions:**
1. Copy raw files to local SSD before processing
2. Increase RAM allocation if running in VM
3. Close other applications to free memory

#### Out of Memory Error (Step 3 or 4)
**Symptoms:**
```
MemoryError: Unable to allocate array
```

**Solutions:**
1. **Increase system RAM** (16 GB minimum recommended)
2. **Process in chunks:**
   ```python
   # Modify MASTER script to process profiles in batches
   # Contact repository maintainer for chunked processing variant
   ```
3. **Use compression:**
   ```python
   ds.to_netcdf('output.nc', encoding={'TEMP': {'zlib': True, 'complevel': 4}})
   ```

### 7.3 Data Quality Issues

#### All CHLA values negative
**Cause:** Dark count not subtracted, or incorrect factory calibration

**Solution:**
1. Check dark count in deployment YAML:
   ```yaml
   chlorophyll:
     dark_count: 50  # Verify with manufacturer calibration sheet
   ```
2. If persists, apply post-processing correction:
   ```python
   import xarray as xr
   ds = xr.open_dataset('output/analysis/sea074-2025_standard_names.nc')
   ds['CHLA'] = ds['CHLA'] - ds['CHLA'].quantile(0.05)  # Subtract 5th percentile
   ds.to_netcdf('output/analysis/sea074-2025_corrected.nc')
   ```

#### DOXY values unrealistic (>500 µmol/kg)
**Cause:** Mixed units not fully corrected

**Solution:**
1. Check Step 3b output messages for conversion statistics
2. Verify all values were multiplied by 1000:
   ```python
   import xarray as xr
   ds = xr.open_dataset('output/l0_data/timeseries/sea074-2025.nc')
   print(f"DOXY range: {ds['oxygen_concentration'].min().values} - {ds['oxygen_concentration'].max().values}")
   # Should be 160-6500 µmol/kg, not 0.16-6.5 mmol/L
   ```

#### Geographic QC all flags = 4 (all out of bounds)
**Cause:** Incorrect lat/lon bounds in configuration

**Solution:**
1. Check actual glider coordinates:
   ```python
   import pandas as pd
   data = pd.read_csv('output/analysis/seaexplorer_data_complete.csv')
   print(f"Lat range: {data['latitude'].min()} to {data['latitude'].max()}")
   print(f"Lon range: {data['longitude'].min()} to {data['longitude'].max()}")
   ```
2. Update `QC_RANGES_LAPALMA` in `scripts/qc_variables.py`:
   ```python
   # Expand bounds to encompass actual mission area
   'LATITUDE': {'min': 28.0, 'max': 29.0},  # Adjust as needed
   'LONGITUDE': {'min': -18.5, 'max': -17.5},
   ```

### 7.4 Getting Help

**Before asking for help, collect:**
1. Error message (full traceback)
2. Python version: `python --version`
3. Package versions: `pip list | grep -E "pyglider|xarray|numpy|pandas"`
4. File structure: `tree -L 2` (or `ls -R`)
5. Configuration file: `cat config/seaexplorer_0067.yml`

**Where to get help:**
- **GitHub Issues:** https://github.com/BennyTorelli/Seaexplorer_pyglider/issues
- **PyGlider Documentation:** https://pyglider.readthedocs.io/
- **Email:** [Your contact email]

**When reporting bugs, include:**
```bash
# Generate diagnostic report
python -c "
import sys
import platform
print('Python version:', sys.version)
print('Platform:', platform.platform())

import pyglider, xarray, numpy, pandas, gsw
print('PyGlider:', pyglider.__version__)
print('xarray:', xarray.__version__)
print('numpy:', numpy.__version__)
print('pandas:', pandas.__version__)
print('gsw:', gsw.__version__)
" > diagnostic_report.txt

# Attach diagnostic_report.txt to GitHub issue
```

---

## 8. Advanced Usage

### 8.1 Batch Processing Multiple Missions

Create a batch script to process multiple deployments:

```python
# batch_process.py
from MASTER_pyglider_pipeline import SeaExplorerProcessor
import os
import glob

# List of missions to process
missions = [
    {'deployment': 'sea074_mission1', 'config': 'config/sea074_feb2025.yml'},
    {'deployment': 'sea074_mission2', 'config': 'config/sea074_mar2025.yml'},
    {'deployment': 'sea075_mission1', 'config': 'config/sea075_feb2025.yml'},
]

for mission in missions:
    print(f"\n{'='*60}")
    print(f"Processing: {mission['deployment']}")
    print(f"{'='*60}\n")
    
    # Set up processor with mission-specific config
    processor = SeaExplorerProcessor()
    processor.deploymentyaml = mission['config']
    
    # Run complete pipeline
    try:
        processor.run_complete_pipeline()
        print(f"✓ {mission['deployment']} completed successfully")
    except Exception as e:
        print(f"✗ {mission['deployment']} failed: {e}")
        continue
    
    # Run QC
    os.system(f"python scripts/qc_variables.py --latest")
    
    # Move outputs to mission-specific folder
    output_dir = f"output/missions/{mission['deployment']}/"
    os.makedirs(output_dir, exist_ok=True)
    os.system(f"mv output/analysis/*.nc {output_dir}")
    os.system(f"mv output/analysis/*.csv {output_dir}")
    
    print(f"✓ Outputs saved to {output_dir}\n")
```

Run batch processing:
```bash
python batch_process.py > batch_log.txt 2>&1
```

### 8.2 Custom QC Thresholds

Create mission-specific QC configuration:

```python
# qc_custom.py
from scripts.qc_variables import *

# Define custom ranges for Arctic deployment
QC_RANGES_ARCTIC = {
    'TEMP': {'min': -2.0, 'max': 15.0, 'unit': '°C'},  # Colder waters
    'PSAL': {'min': 30.0, 'max': 35.0, 'unit': 'PSU'},  # Lower salinity
    'CNDC': {'min': 0.0, 'max': 5.0, 'unit': 'S/m'},
    'DOXY': {'min': 200.0, 'max': 400.0, 'unit': 'µmol/L'},  # Higher oxygen
    'CHLA': {'min': 0.0, 'max': 20.0, 'unit': 'mg/m³'},
    'TURB': {'min': 0.0, 'max': 2.0, 'unit': 'NTU'},
    'PRES': {'min': 0.0, 'max': 1000.0, 'unit': 'dbar'},  # Shallower deployment
}

# Use custom ranges in QC
# ... (modify qc_variables.py to accept custom ranges)
```

### 8.3 Integration with Data Visualization Tools

**Ocean Data View (ODV):**
```bash
# Export to ODV-compatible format
python -c "
import pandas as pd
import numpy as np

data = pd.read_csv('output/analysis/seaexplorer_data_complete_standard_names.csv')

# ODV required columns: Cruise, Station, Type, yyyy-mm-ddThh:mm:ss, Longitude, Latitude, Depth
odv = pd.DataFrame({
    'Cruise': 'sea074',
    'Station': np.arange(len(data)),
    'Type': 'G',  # Glider
    'yyyy-mm-ddThh:mm:ss': data['TIME'],
    'Longitude [degrees_east]': data['LONGITUDE'],
    'Latitude [degrees_north]': data['LATITUDE'],
    'Depth [m]': data['DEPTH'],
    'Temperature [degC]': data['TEMP'],
    'Salinity [PSU]': data['PSAL'],
    'Oxygen [umol/kg]': data['DOXY'],
    'Chlorophyll [mg/m^3]': data['CHLA'],
})

odv.to_csv('output/analysis/sea074_ODV.txt', sep='\t', index=False)
print('ODV file created: output/analysis/sea074_ODV.txt')
"
```

**QGIS (Geographic Information System):**
1. Open QGIS
2. Layer → Add Layer → Add Delimited Text Layer
3. Select `seaexplorer_data_complete_standard_names.csv`
4. Set X field: `LONGITUDE`, Y field: `LATITUDE`
5. CRS: EPSG:4326 (WGS84)
6. Style by depth or variable values

**Jupyter Notebook Interactive Analysis:**
```python
# analysis.ipynb
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean  # Oceanographic colormaps

# Load data
ds = xr.open_dataset('output/analysis/sea074-2025_standard_names.nc')

# Interactive plot with hover tooltips
import holoviews as hv
from holoviews import opts
hv.extension('bokeh')

points = hv.Points(ds, ['LONGITUDE', 'LATITUDE'], ['TEMP', 'DEPTH'])
points.opts(color='TEMP', cmap='thermal', size=3, tools=['hover'], 
           width=800, height=600, colorbar=True)
```

### 8.4 Exporting to Data Repositories

**ERDDAP Preparation:**
```python
# Ensure CF-1.8 compliance
import xarray as xr

ds = xr.open_dataset('output/analysis/sea074-2025_standard_names.nc')

# Add required global attributes
ds.attrs['cdm_data_type'] = 'Trajectory'
ds.attrs['featureType'] = 'trajectory'
ds.attrs['Conventions'] = 'CF-1.8, ACDD-1.3'
ds.attrs['keywords'] = 'ocean glider, SeaExplorer, temperature, salinity, oxygen, chlorophyll'
ds.attrs['keywords_vocabulary'] = 'GCMD Science Keywords'
ds.attrs['standard_name_vocabulary'] = 'CF Standard Name Table v79'
ds.attrs['license'] = 'These data may be used freely for research and educational purposes'

# Add trajectory coordinate
ds['trajectory'] = ([], 'sea074')
ds['trajectory'].attrs['cf_role'] = 'trajectory_id'
ds['trajectory'].attrs['long_name'] = 'trajectory identifier'

# Save ERDDAP-ready file
ds.to_netcdf('output/analysis/sea074-2025_ERDDAP.nc')
```

**SeaDataNet CDI/CSR Export:**
```python
# Generate CDI (Common Data Index) metadata
cdi_metadata = {
    'edmo_code': '1234',  # Your institution's EDMO code
    'local_cdi_id': 'sea074-2025',
    'cdi_version': '1.0',
    'platform_code': 'SEA074',
    'platform_name': 'SeaExplorer glider sea074',
    'cruise_name': 'DELTA Mission Feb-Mar 2025',
    'chief_scientist': 'Your Name',
    'originator': 'Your Institution',
    'data_type': 'Ocean glider profile',
    'geospatial_lat_min': float(ds['LATITUDE'].min()),
    'geospatial_lat_max': float(ds['LATITUDE'].max()),
    'geospatial_lon_min': float(ds['LONGITUDE'].min()),
    'geospatial_lon_max': float(ds['LONGITUDE'].max()),
    'depth_min': 0.0,
    'depth_max': float(ds['DEPTH'].max()),
    'time_coverage_start': str(ds['TIME'].values[0]),
    'time_coverage_end': str(ds['TIME'].values[-1]),
    'parameters': 'TEMP, PSAL, DOXY, CHLA, TURB',
    'doi': ''  # To be assigned
}

import json
with open('output/analysis/sea074_CDI_metadata.json', 'w') as f:
    json.dump(cdi_metadata, f, indent=2)
```

---

## Appendix A: File Format Specifications

### A.1 Raw File Formats

**Payload files (.pld1, .pld0):**
- Binary format, proprietary to ALSEAMAR
- Contains: sensor data, timestamps, navigation
- Processed by PyGlider raw_to_rawnc()

**Navigation files (.gli):**
- ASCII text format
- Contains: GPS fixes, surface data
- Merged with payload data during processing

### A.2 Intermediate Formats

**Parquet files:**
- Columnar storage format (Apache Parquet)
- Fast read/write, efficient compression
- Used internally by PyGlider

### A.3 Output Formats

**NetCDF-4:**
- Self-describing binary format
- CF-1.8 compliant
- Contains data + metadata
- Tools: ncdump, ncview, Python/xarray, MATLAB

**CSV:**
- Plain text, comma-separated
- Universal compatibility
- Large file sizes (no compression)
- Tools: Excel, R, Python/pandas, QGIS

---

## Appendix B: Coordinate Systems

**Temporal:**
- Reference: UTC (Coordinated Universal Time)
- Format: ISO 8601 (YYYY-MM-DDTHH:MM:SS)
- Encoding: datetime64[ns] (nanosecond precision)

**Spatial:**
- Horizontal: WGS84 (EPSG:4326)
  - Latitude: degrees_north, -90 to +90
  - Longitude: degrees_east, -180 to +180
- Vertical: Depth below sea surface (positive down)
  - Units: meters
  - Reference: local sea surface (not geoid)

**Pressure:**
- Units: decibars (dbar)
- Conversion: 1 dbar ≈ 1 meter depth in ocean
- Absolute pressure (includes atmospheric)

---

## Appendix C: Version History

**v1.0 (October 2025):**
- Initial release
- Complete pipeline automation
- 30-column QC system
- LAND_QC implementation
- TEOS-10 conversions
- Variable standardization

**Planned for v1.1:**
- Real-time processing mode
- Automated dark count correction
- Multi-glider merging tools
- Enhanced visualization gallery

---

**Document Version:** 1.0  
**Last Updated:** October 2025  
**Maintainer:** Benedetta Torelli  
**Repository:** https://github.com/BennyTorelli/Seaexplorer_pyglider
