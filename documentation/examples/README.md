# Documentation Examples

This directory contains practical Python examples demonstrating how to use the SeaExplorer data processing pipeline.

## Available Examples

### 1. Basic Pipeline Usage (`basic_pipeline_usage.py`)

**Purpose:** Complete end-to-end processing from raw glider files to quality-controlled outputs.

**What it does:**
- Initializes the SeaExplorerProcessor
- Runs all 7 processing steps (raw → L0 → standardization)
- Executes quality control analysis
- Verifies output files
- Displays data summary statistics

**Usage:**
```bash
cd /Users/benedettatorelli/Desktop/PyGlider_SeaExplorer_Project
python documentation/examples/basic_pipeline_usage.py
```

**Expected runtime:** 20-30 minutes (for ~555 raw files)

**Output files:**
- `output/analysis/sea074-2025_standard_names.nc` (NetCDF)
- `output/analysis/seaexplorer_data_complete_standard_names.csv` (full dataset)
- `output/analysis/seaexplorer_qc_variables.csv` (QC flags)

---

### 2. Quality Control Analysis (`qc_analysis_example.py`)

**Purpose:** Comprehensive QC results analysis and data filtering.

**What it does:**
- Loads QC CSV and generates detailed statistics
- Identifies variables with quality issues (>5% BAD flags)
- Creates 3 filtered datasets (strict, moderate, liberal)
- Performs temporal pattern analysis (biofouling detection)
- Generates QC visualization plots
- Produces QC summary report

**Usage:**
```bash
# First ensure QC has been run
python scripts/qc_variables.py --latest

# Then run the analysis
python documentation/examples/qc_analysis_example.py
```

**Expected runtime:** 1-2 minutes

**Output files:**
- `seaexplorer_data_QC_strict.csv` (most conservative, ~80-95% data retained)
- `seaexplorer_data_QC_moderate.csv` (recommended, ~90-98% retained)
- `seaexplorer_data_QC_liberal.csv` (most permissive, ~95-99% retained)
- `qc_timeseries_overview.png` (visualization)
- `QC_Report.txt` (detailed text report)

---

## Requirements

All examples require:
- Completed pipeline execution (raw data processed)
- Python 3.11+ with installed dependencies (see `requirements.txt`)
- Sufficient disk space (~2-3 GB for outputs)

**Install dependencies:**
```bash
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Example Workflow

**Complete processing from scratch:**

```bash
# 1. Activate virtual environment
cd /Users/benedettatorelli/Desktop/PyGlider_SeaExplorer_Project
source .venv/bin/activate

# 2. Ensure raw data is in input/raw/
ls -l input/raw/*.pld1
ls -l input/raw/*.gli

# 3. Run complete pipeline
python documentation/examples/basic_pipeline_usage.py

# 4. Analyze QC results
python documentation/examples/qc_analysis_example.py

# 5. Use filtered data for your analysis
# Choose appropriate filtered CSV based on your needs:
#   - High-precision studies: seaexplorer_data_QC_strict.csv
#   - General analysis: seaexplorer_data_QC_moderate.csv (RECOMMENDED)
#   - Visualization: seaexplorer_data_QC_liberal.csv
```

---

## Customization

### Adapting for Your Deployment

**1. Update configuration file:**
```bash
cp config/seaexplorer_0067.yml config/seaexplorer_YOUR_GLIDER.yml
# Edit with your mission details (dates, location, sensor calibrations)
```

**2. Modify QC thresholds:**
```python
# In scripts/qc_variables.py, lines 12-22
QC_RANGES_YOUR_REGION = {
    'TEMP': {'min': YOUR_MIN, 'max': YOUR_MAX, 'unit': '°C'},
    'PSAL': {'min': YOUR_MIN, 'max': YOUR_MAX, 'unit': 'PSU'},
    # ... adjust for your oceanographic region
}
```

**3. Update geographic boundary:**
```python
# In scripts/qc_variables.py or examples
# For LAND_QC: provide shapefile for your coastline
# For Location_QC: adjust lat/lon bounds
```

---

## Additional Examples (Coming Soon)

- **Data Visualization:** Creating publication-quality plots
- **T-S Diagram Analysis:** Water mass identification
- **Oxygen Saturation Profiles:** Calculating and plotting O₂%
- **Batch Processing:** Multiple missions in one script
- **ERDDAP Export:** Preparing data for online repositories
- **Integration with R:** Exporting for statistical analysis

---

## Troubleshooting

### Example fails with "No module named 'MASTER_pyglider_pipeline'"

**Solution:**
```python
# Add project root to Python path at top of script
import sys
sys.path.append('/Users/benedettatorelli/Desktop/PyGlider_SeaExplorer_Project')
```

### Example fails with "No raw files found"

**Solution:**
```bash
# Ensure raw data is present
ls input/raw/*.pld1

# If empty, copy raw files:
cp /path/to/your/raw/data/*.pld1 input/raw/
cp /path/to/your/raw/data/*.gli input/raw/
```

### Example fails with "FileNotFoundError: qc_variables.csv"

**Solution:**
```bash
# Run QC script first
python scripts/qc_variables.py --latest
```

---

## Contact & Support

- **Repository:** https://github.com/BennyTorelli/Seaexplorer_pyglider
- **Issues:** https://github.com/BennyTorelli/Seaexplorer_pyglider/issues
- **Documentation:** See `documentation/` directory for comprehensive guides

---

**License:** MIT  
**Author:** Benedetta Torelli  
**Last Updated:** October 2025
