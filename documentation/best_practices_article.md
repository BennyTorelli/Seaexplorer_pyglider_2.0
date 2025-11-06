# Best Practices for Processing and Quality Control of SeaExplorer Glider Data Using an Open-Source Python Pipeline

**Article Type:** Best Practices in Ocean Observing  
**Target Journal:** Frontiers in Marine Science

---

## Abstract

Autonomous underwater gliders have become essential tools for sustained ocean observations, but the lack of standardized, end-to-end data processing and quality control workflows remains a significant challenge for the scientific community. This article presents a comprehensive best-practice framework for transforming raw SeaExplorer glider data into analysis-ready, high-quality datasets. We detail a complete pipeline, built upon the open-source `PyGlider` toolbox, which automates the entire workflow from initial data ingestion to the generation of science-ready products. The process includes critical custom enhancements for unit standardization according to TEOS-10 oceanographic conventions and a multi-tiered Quality Control (QC) system. This QC framework integrates established methods with novel procedures specifically designed for glider data, such as automated coastal proximity detection (LAND_QC), advanced spike detection for physical and bio-optical sensors, and density inversion checks. We demonstrate the workflow using data from a 21-day SeaExplorer mission near La Palma, Canary Islands, processing over 1.8 million data points. The results highlight the pipeline's capability to correct unit inconsistencies, identify sensor anomalies, and produce fully-documented, reproducible datasets. All code, configuration files, and methodologies are openly available, providing a robust and reusable template for the glider community to enhance data quality, ensure interoperability, and adhere to FAIR data principles.

**Keywords:** autonomous underwater glider, SeaExplorer, PyGlider, data processing, quality control, oceanographic data, best practices, TEOS-10, open science, reproducibility

---

## 1. Introduction

### 1.1 The Rise of Glider Oceanography and the Data Challenge

Autonomous underwater gliders have revolutionized oceanographic observation over the past two decades (Rudnick, 2016; Testor et al., 2019). These platforms provide sustained, high-resolution measurements of physical and biogeochemical properties across vast ocean regions, complementing traditional ship-based surveys and fixed moorings. The SeaExplorer glider (manufactured by ALSEAMAR), in particular, has proven to be a versatile and robust platform for extended missions in differents ocean environments, thanks to its high payload capacity and energy efficiency.

However, the proliferation of glider missions has created a new challenge: the management and processing of the vast and complex datasets they generate. Raw glider data often suffer from a range of issues that hinder direct scientific analysis:
1.  **Proprietary and Complex File Formats:** Data are often stored in multiple, instrument-specific binary files.
2.  **Inconsistent Naming and Units:** Variables may lack standardized names, and units can be mixed or non-standard (e.g., mS/cm vs. S/m).
3.  **Data Quality Issues:** Raw data streams are susceptible to errors from sensor noise, biofouling, environmental contamination, or transmission failures.
4.  **Lack of Standardization:** The absence of a unified processing workflow makes it difficult to compare data across different missions or institutions, hindering large-scale oceanographic synthesis.

### 1.2 Bridging the Gap with Open-Source Tools

The open-source Python package `PyGlider` (Gregor et al., 2019) has made significant strides in addressing some of these issues by providing a foundational framework for converting raw data from multiple glider types (including Seaglider, Slocum, and SeaExplorer) into a more usable, standardized netCDF format. However, `PyGlider` is a toolbox, not a complete, end-to-end solution. The production of truly "science-ready" data requires significant additional steps, particularly in the areas of unit standardization and rigorous, multi-faceted quality control.

### 1.3 Objectives of This Best-Practice Guide

This article details a comprehensive, best-practice workflow that extends the capabilities of `PyGlider` to create a fully automated and reproducible data processing pipeline. Our primary objectives are to:
1.  **Provide an End-to-End Solution:** Document a complete data processing chain, from raw file decompression to the generation of analysis-ready products.
2.  **Implement Rigorous Quality Control:** Describe a multi-tiered QC framework that applies both established and novel tests to ensure data integrity, following international standards (e.g., QARTOD, Argo).
3.  **Ensure Oceanographic Accuracy:** Detail the critical steps for standardizing units and recalculating key variables (salinity, density) according to the Thermodynamic Equation of Seawater 2010 (TEOS-10).
4.  **Promote Reproducibility and Open Science:** Offer a fully documented and open-source codebase that serves as a reusable template for the oceanographic community.

We demonstrate this workflow using a real-world dataset from a SeaExplorer mission, providing a practical guide for other research groups to adopt and adapt for their own deployments.

---

## 2. Materials and Methods

### 2.1 The La Palma Mission: A Case Study

To illustrate our pipeline, we use data from a representative SeaExplorer mission conducted off the coast of La Palma, Canary Islands, Spain.

*   **Location:** Waters surrounding La Palma, Canary Islands (approx. 28.45°N - 28.86°N, 18.01°W - 17.72°W)
*   **Mission Duration:** February 12 - March 5, 2025 (21 days)
*   **Platform:** SeaExplorer glider (ID: sea074)
*   **Key Statistics:** The mission yielded 358 vertical profiles (179 down-casts, 179 up-casts), reaching a maximum depth of 1,030 m and collecting a total of 1,825,662 multi-parameter measurements.

### 2.2 Sensor Suite

The SeaExplorer was equipped with a standard suite of physical and biogeochemical sensors, whose characteristics are crucial for defining QC parameters.

*   **Physical:** SBE 41CP CTD (Conductivity, Temperature, Depth).
*   **Biogeochemical:** Aanderaa Optode 4330 (Dissolved Oxygen) and a WET Labs ECO Puck measuring Chlorophyll-a fluorescence, CDOM fluorescence, and Optical Backscatter (700 nm).

### 2.3 The Software Ecosystem

Our pipeline is built entirely on open-source software, ensuring accessibility and reproducibility.

*   **Python Version:** 3.13.5
*   **Core Libraries:**
    *   `PyGlider` (0.0.7): For core data conversion.
    *   `xarray` (2025.9.0): For handling multi-dimensional labeled data (netCDF).
    *   `pandas` (2.3.2): For data manipulation and CSV export.
    *   `gsw` (3.6.20): For TEOS-10 thermodynamic calculations of seawater properties.
    *   `shapely` (2.1.2): For computational geometry, used in our coastal proximity QC.
*   **Repository:** All scripts and configuration files are available at: https://github.com/BennyTorelli/Seaexplorer_pyglider

### 2.4 The Data Processing Pipeline: A Discursive Walkthrough

The complete processing workflow is orchestrated by a master script, `MASTER_pyglider_pipeline.py`, which executes a sequence of steps designed to be both modular and robust. Here, we describe the rationale and function of each step.

#### Step 0: File Preparation and Sanitization

**Rationale:** Raw data from the glider is often delivered in a compressed format and may contain formatting inconsistencies that can cause processing tools to fail. This initial step prepares the data for robust parsing.

**Implementation:**
1.  **Decompression:** The pipeline first scans the raw data directory for any `.gz` files and decompresses them. This is a necessary prerequisite as `PyGlider` operates on the uncompressed text-based log files.
2.  **Sanitization:** A crucial custom function, `_sanitize_raw_logs`, is executed. It reads the raw text files and specifically targets the timestamp columns (e.g., `PLD_REALTIMECLOCK`). It searches for empty or invalid timestamp entries (e.g., 'NULL', '-') and replaces them with a placeholder date ('01/01/1971 00:00:00.000'). This prevents fatal errors during the date-parsing stage in `PyGlider`, which expects a valid datetime format in every row.

#### Step 1 & 2: Core PyGlider Conversion and Merging

**Rationale:** This is the core function of `PyGlider`, converting the numerous, fragmented raw text files into a single, consolidated dataset.

**Implementation:**
1.  **Raw to Parquet (`raw_to_rawnc`):** The `seaexplorer.raw_to_rawnc` function is called. It reads all the sanitized `.pld1`, `.gli`, and `.pld0` files, parses them according to the sensor definitions in the `deployment.yml` file, and converts them into an efficient intermediate `parquet` format.
    ```python
    # From MASTER_pyglider_pipeline.py
    seaexplorer.raw_to_rawnc(
        rawdir='input/sanitized/',
        rawncdir='output/l0_data/rawnc/',
        deploymentyaml='config/seaexplorer_0067.yml'
    )
    ```
2.  **Merge Parquet (`merge_parquet`):** The resulting collection of `parquet` files is then merged into a single, time-ordered master `parquet` file. This step consolidates all sensor data into a unified timeline.

#### Step 3: Generation of Level-0 Time-Series NetCDF

**Rationale:** The merged `parquet` file is an intermediate product. The goal is to create a standard, self-describing netCDF file that adheres to Climate and Forecast (CF) conventions. This is the first "usable" data product.

**Implementation:** The `seaexplorer.raw_to_timeseries` function is used to convert the merged `parquet` file into a Level-0 (L0) time-series netCDF file. This file contains all variables with their initial calibrations applied, but before advanced unit conversions and quality control.

#### Step 3b: Custom Unit Conversions and Recalculations (Critical Enhancement)

**Rationale:** This is one of the most critical contributions of our pipeline. Raw sensor outputs are often not in standard scientific units, which prevents accurate oceanographic calculations. This step corrects units and recalculates key derived parameters using the TEOS-10 standard, ensuring the data are physically consistent and comparable.

**Implementation:** A custom function, `_apply_custom_conversions`, opens the L0 netCDF file and performs the following operations directly on the `xarray.Dataset`:

1.  **Turbidity (NTU):** Converts optical backscatter (units of m⁻¹sr⁻¹) to Nephelometric Turbidity Units (NTU) using the manufacturer's scale factor. NTU is the community-standard unit for turbidity.
2.  **Conductivity (S/m):** Converts conductivity from milliSiemens per centimeter (mS/cm) to the SI unit of Siemens per meter (S/m) by dividing by 10. This is essential for accurate TEOS-10 calculations.
3.  **Dissolved Oxygen (µmol/L to µmol/kg):** This involves a two-part correction. First, it addresses an instrument-specific issue where data units were mixed (some values in µmol/L, others in mmol/L). The script identifies and scales the mmol/L values by 1000. Second, it converts the volume-based concentration (µmol/L) to a mass-based concentration (µmol/kg) using the newly calculated in-situ density. This makes the measurement independent of pressure and temperature effects.
4.  **Salinity (TEOS-10 Recalculation):** Practical Salinity is recalculated from scratch using the corrected conductivity (S/m), in-situ temperature, and pressure via the `gsw.SP_from_C` function. This removes any legacy calculations and ensures full compliance with TEOS-10.
5.  **Density (TEOS-10 Recalculation):** In-situ density is recalculated using the new Practical Salinity, in-situ temperature, and pressure. The calculation first derives Absolute Salinity (`SA`) and Conservative Temperature (`CT`), then uses `gsw.rho` to compute the most accurate in-situ density.

#### Step 4 & 6: Final Product Generation

**Rationale:** To maximize usability, the pipeline generates final products in standard formats with clear, community-accepted variable names.

**Implementation:**
1.  **CSV Export:** The fully processed netCDF file is exported to a CSV file (`seaexplorer_data_complete.csv`). This format is easily accessible for users who may not be familiar with netCDF or `xarray`.
2.  **Variable Standardization:** A mapping is applied to rename all variables to a consistent, uppercase convention (e.g., `temperature` becomes `TEMP`, `oxygen_concentration` becomes `DOXY`). This aligns with conventions used in major oceanographic databases (e.g., Argo). The final, standardized files are saved with a `_standard_names` suffix.

---

## 3. The Multi-Tier Quality Control Framework: Ensuring Data Integrity

The generation of a standardized, TEOS-10 compliant dataset is only the first step. The true scientific value of glider data is unlocked only after a rigorous and transparent assessment of its quality. Data from autonomous platforms are susceptible to a myriad of issues, from sensor biofouling and calibration drift to environmental contamination and electronic noise. Failure to identify and flag these issues can lead to erroneous scientific conclusions.

Therefore, a systematic Quality Control (QC) process is not an optional post-processing step but a fundamental component of the data life cycle. The philosophy of our framework is to **flag, not remove**, providing the end-user with a comprehensive set of indicators to make informed decisions. This approach preserves the original data while offering clear guidance on its reliability.

Our pipeline executes a dedicated QC script (`scripts/qc_variables.py`) that applies a battery of tests to the processed data. To ensure interoperability with global data systems like Argo and EMODnet, we use a simplified version of the **IODE/ARGO QC flag scheme**:

*   **1 (GOOD):** The data point has passed all relevant QC tests and is considered reliable.
*   **4 (BAD):** The data point has failed at least one critical QC test and is considered erroneous or highly suspect.
*   **9 (MISSING):** The data point was not recorded or is unavailable (`NaN`).
*   **0 (NOT EVALUATED):** A specific QC test was not performed on the data point, often due to missing contextual data (e.g., at the edge of a profile).

The following sections provide a detailed, test-by-test description of the framework, from foundational checks to advanced, platform-specific diagnostics.

### 3.1 Foundational Checks: Time, Position, and Data Gaps

These initial tests form the bedrock of the QC process, ensuring that the data has a valid temporal and spatial context.

*   **Granular Temporal QC (`Date_QC`):** A robust check on the integrity of the timestamp is performed. Our enhanced methodology moves beyond a simple range check to a **granular validation** of each component of the timestamp. For every data point, the pipeline verifies:
    *   **Year:** Plausible range (> 1990 and not in the future).
    *   **Month:** Valid range (1-12).
    *   **Day:** Valid for the given month and year, correctly handling leap years.
    *   **Hour, Minute, Second:** Within their valid ranges (0-23, 0-59, 0-59).
    This fine-grained approach is critical for detecting subtle data corruption or formatting errors that a coarse bounding box might miss, ensuring the temporal accuracy required for calculating derived parameters like velocity.

*   **Geographic QC (`Location_QC`):** A "bounding box" test to ensure the glider's reported latitude and longitude are within the expected operational area for the mission (e.g., the Canary Islands region). This is a first-pass check to catch gross position errors.

*   **Coastal Proximity QC (`LAND_QC`):** A more sophisticated spatial test, critical for coastal missions. This test uses a high-resolution shoreline shapefile (in this case, for La Palma) to determine if a glider's reported position is physically on land. It calculates the distance to the nearest coast and flags any point located within the land polygon as **bad (4)**. This prevents the erroneous interpretation of data collected while the glider may have been grounded or reporting inaccurate GPS fixes near the coast.

*   **Missing Value QC (`Na_QC`):** A simple but vital test that flags every `NaN` (Not a Number) value with a `9`. This explicitly distinguishes missing data from measured data and ensures that data gaps are not misinterpreted.

### 3.2 Sensor and Environmental QC: From Plausibility to Dynamics

This tier of tests scrutinizes the sensor data itself, checking for physical plausibility, sensor malfunctions, and environmental contamination.

*   **Regional Range Test (`{VAR}_Range_QC`):** This test has been specifically tailored for the mission area. Instead of a generic "global" range, we apply a **regional climatological check** based on known oceanographic conditions around the Canary Islands. For each variable (e.g., `TEMP`, `PSAL`), we define a plausible range derived from historical data (e.g., World Ocean Atlas) and expert knowledge of the region. For example, temperature values are expected to fall within a specific range characteristic of the North Atlantic Central Water. A value falling outside this regional box is flagged as **bad (4)**. This is more powerful than a global test as it can detect anomalous values that might still be globally plausible but are unrealistic for the specific study area.

*   **Stuck Value Test (`{VAR}_Stuck_QC`):** This test is a crucial diagnostic for sensor failure. It identifies situations where a sensor reports the exact same value for a prolonged period within a single profile (ascent or descent). The script segments the data into profiles based on pressure changes and then checks for constant, non-`NaN` values within each segment. If found, all points in that segment for that variable are flagged as **bad (4)**, as this indicates a "stuck" sensor that is no longer responding to environmental changes.

*   **Spike Test (3-Point Median for Physical Sensors):** A standard test for physical variables (`TEMP`, `PSAL`, `POTDEN`) that flags a point if it represents a significant, non-physical deviation from its two immediate neighbors. A point `p(i)` is flagged as a spike if `|p(i) - median(p(i-1), p(i), p(i+1))|` exceeds a pre-defined threshold. This is effective at catching transient electronic noise.

*   **Derived Velocity and Physical Plausibility (`VELOCITY_QC`):** A significant enhancement to our pipeline is the calculation and quality control of the glider's vertical speed. This serves as a powerful diagnostic for both data quality and glider performance.
    *   **Calculation:** The vertical velocity (`W`) is first calculated as the rate of change of depth over time (`dDEPTH/dTIME`). This is then converted to the glider's speed through the water (`VELOCITY`) by correcting for the glider's pitch angle (`PITCH`), which is recorded by the glider's internal flight computer. The formula is:
        $$ VELOCITY_{cm/s} = \left| \frac{dDEPTH/dTIME}{\sin(\text{PITCH}_{\text{radians}})} \right| \times 100 $$
        The absolute value is used to represent speed as a positive quantity, irrespective of ascent or descent. The `PITCH` variable itself is also included in the final data product, providing valuable context on the glider's flight characteristics.
    *   **Quality Control:** A physical plausibility test is applied to the calculated `VELOCITY`. A data point is flagged as **bad (4)** if the speed exceeds a threshold of **401 cm/s**. This threshold is derived from a conservative estimate of the maximum possible speed, considering the glider's own maximum propulsion speed (approx. 51 cm/s) and the most extreme ocean currents ever recorded (e.g., Gulf Stream, Kuroshio, at approx. 350 cm/s). Speeds exceeding this limit are considered physically unrealistic and likely indicate issues with the underlying `DEPTH`, `TIME`, or `PITCH` data points used in the calculation.

### 3.3 Bio-Optical Sensor QC: Advanced Spike and Contamination Detection

Bio-optical sensors, such as those measuring Chlorophyll-a (`CHLA`) and Turbidity (`TURB`), are inherently noisier than CTD sensors. They are prone to significant, sharp spikes caused by encounters with marine aggregates, large plankton, or internal sensor noise. This requires a more specialized set of QC tests.

*   **Negative Value Test:** A fundamental check that flags any non-physical negative readings from these sensors as **bad (4)**. This is often indicative of a calibration issue or an incorrect dark count offset.

*   **Enhanced Spike Test (5-Point Strict Median):** To handle the specific noise characteristics of these sensors, we have implemented a more robust spike detection algorithm that differs from the one used for physical sensors.
    *   **Methodology:** This test uses a 5-point median filter. A point `p(i)` is evaluated against the median of a window including its two preceding and two succeeding neighbors: `median(p(i-2), p(i-1), p(i), p(i+1), p(i+2))`. A wider window is used to better capture the baseline signal around potentially noisy spikes.
    *   **Strict Neighbor Requirement:** A crucial enhancement is our strict handling of `NaN` values within the moving window. If **any** of the four neighbors (`p(i-2)`, `p(i-1)`, `p(i+1)`, `p(i+2)`) is `NaN`, the central point `p(i)` is not evaluated and flagged as **0 (not evaluated)**. This prevents the erroneous flagging of data points at the edges of data gaps, ensuring that spikes are only identified when there is a continuous, high-quality local context. This strictness significantly reduces false positives compared to more lenient methods that might interpolate across gaps.

*   **Surface Contamination QC (for Bio-optics):** Bio-optical measurements in the top few meters of the water column are highly susceptible to contamination. This can be caused by bubbles injected during the glider's submersion, interference from surface wave action, or biological phenomena like non-photochemical quenching (NPQ) in high sunlight. To account for this, this test prophylactically flags all `CHLA` and `TURB` data collected at pressures less than 5 dbar (~5 meters depth) as **bad (4)**. While this is a conservative measure, it ensures that data used for primary productivity or water clarity studies are free from near-surface artifacts.

### 2.6 Final Output Products

The combined pipeline generates three primary, user-friendly output products:

1.  **L0 NetCDF Time-Series (`..._standard_names.nc`):** A CF-1.8 compliant netCDF file containing the fully processed, time-ordered data with all unit conversions and standardized variable names applied. This is the primary product for scientific analysis with tools that support netCDF.
2.  **Complete Data CSV (`..._standard_names.csv`):** A flattened CSV file containing all the core data, intended for easy import into spreadsheets, GIS, or statistical software.
3.  **Quality Control CSV (`seaexplorer_qc_variables.csv`):** A comprehensive CSV where each row corresponds to a measurement and each column represents a specific QC test for a specific variable. This file provides the full matrix of QC flags, allowing for detailed data screening and quality assessment.

---

## 3. Results

### 3.1 Pipeline Performance and Data Throughput

The automated pipeline demonstrates high efficiency. Running on a standard laptop (MacBook Air M2, 16GB RAM), the entire workflow—from raw file decompression to the generation of final, quality-controlled products—processed the 21-day mission dataset (1.8 million measurements) in approximately **25 minutes**. This rapid turnaround is critical for operational oceanography and near-real-time data delivery.

### 3.2 Impact of Unit Conversions and Recalculations

The custom conversion step (3b) proved essential for data integrity.
*   **Conductivity:** Initial values were in mS/cm. Conversion to S/m yielded a realistic range of 0.77 - 5.16 S/m.
*   **Salinity:** Recalculation using TEOS-10 produced a median Practical Salinity of 35.97 PSU, consistent with the known oceanography of the Eastern North Atlantic and in agreement with the World Ocean Atlas climatology.
*   **Density:** The TEOS-10 recalculation corrected numerous non-physical values in the original data, producing a physically plausible range of 996.2 - 1032.3 kg/m³.
*   **Dissolved Oxygen:** The pipeline successfully identified and corrected over 1.8 million values that were stored in inconsistent units (mmol/L vs. µmol/L), then converted them to mass-based units (µmol/kg), which is the standard for oceanographic databases.

### 3.3 Key Findings from the Quality Control Analysis

The multi-tier QC framework provided a detailed assessment of data quality.
*   **Temporal and Geographic Validity:** All 1.8 million data points passed the `Date_QC`, `Location_QC`, and the novel `LAND_QC`, indicating no gross errors in time or position.
*   **Physical Data Quality:** The core physical variables (TEMP, CNDC, PSAL) showed excellent quality, with over 99.9% of data passing the climatological range checks.
*   **Bio-Optical Sensor Issues:** The `CHLA_Sensor_QC` flagged 12.4% of chlorophyll values as bad. This was traced to a known issue of negative values resulting from an uncorrected dark count offset in the sensor, a common issue that this QC test successfully identified.
*   **Spike Detection:** The spike tests identified several hundred isolated spikes in the TEMP, PSAL, and DOXY data, successfully cleaning the dataset of transient electronic noise.
*   **Stuck Sensor:** The `Stuck_QC` test did not find any instances of stuck sensors, indicating healthy sensor performance throughout the mission.

---

## 4. Discussion

### 4.1 Advantages of an Integrated and Automated Workflow

The presented pipeline offers several key advantages over manual or fragmented processing approaches.

*   **Automation and Efficiency:** The ability to process an entire multi-week mission in under 30 minutes dramatically accelerates the data-to-science timeline.
*   **Reproducibility and Traceability:** By encapsulating all processing steps and parameters in version-controlled scripts and configuration files, the entire workflow is fully reproducible. Any scientist can rerun the pipeline on the same raw data and obtain identical results.
*   **Standardization and Interoperability:** The strict adherence to TEOS-10 standards and CF-compliant metadata ensures that the output products are immediately compatible with international data repositories (e.g., CMEMS, EMODnet, GOOS) and standard analysis software (e.g., Ocean Data View).
*   **Enhanced Quality Assurance:** The comprehensive, 30+ flag QC system provides a granular and transparent assessment of data quality. The novel `LAND_QC` provides a critical, previously missing check for coastal missions.

### 4.2 The Importance of Discursive and Granular QC

A key philosophy of this framework is that quality control should be as detailed and informative as possible. Instead of a single, generic "QC" flag, our pipeline generates specific flags for specific tests (e.g., `TEMP_Range_QC`, `TEMP_Spike_QC`, `TEMP_Sensor_QC`). This allows scientists to understand *why* a data point was flagged. For example, a value flagged by `Range_QC` might indicate an interesting oceanographic anomaly, whereas a value flagged by `Sensor_QC` is almost certainly an instrument malfunction. This level of detail is crucial for robust scientific interpretation.

### 4.3 Limitations and Future Directions

While robust, the pipeline has areas for future improvement.
*   **Automated Dark Count Correction:** The workflow currently flags negative chlorophyll values but does not automatically correct them. Future work should implement an automated dark count correction based on measurements from the glider's deep profiles.
*   **Real-Time Implementation:** The current pipeline is designed for post-mission processing. Adapting it for real-time operations would require a shift from batch processing to a streaming data architecture.
*   **Advanced QC Tests:** Additional tests, such as frequency-domain analysis for detecting sensor noise or inter-variable checks (e.g., comparing optical properties to biogeochemical parameters), could be integrated.

### 4.4 Recommendations for the Glider Community

Based on our experience developing this workflow, we offer the following recommendations:
1.  **Adopt Automated Pipelines:** Move away from manual, step-by-step processing and embrace fully automated, scripted workflows.
2.  **Prioritize Unit Standardization:** Ensure all oceanographic variables are converted to community-accepted standard units (e.g., TEOS-10) as an early step in the processing chain.
3.  **Implement Granular QC:** Use a multi-flag QC system that provides detailed information on why data are considered suspect.
4.  **Embrace Open Science:** Share processing code and configuration files alongside data to ensure full reproducibility and allow for community-driven improvement.

---

## 5. Conclusions

We have presented a comprehensive, reproducible, and open-source workflow for processing SeaExplorer glider data. The pipeline provides an end-to-end solution that transforms raw, complex instrument files into standardized, quality-controlled, and science-ready datasets. Key innovations include the integration of critical TEOS-10 unit conversions and a multi-tiered QC framework featuring novel tests like automated coastal proximity detection.

By providing this workflow as a fully documented, open-source template, we aim to help standardize glider data processing across the oceanographic community. Adopting such practices is essential for improving data quality, fostering collaboration, and ultimately maximizing the scientific return from the global fleet of ocean gliders.

All materials are freely available at: https://github.com/BennyTorelli/Seaexplorer_pyglider

---
