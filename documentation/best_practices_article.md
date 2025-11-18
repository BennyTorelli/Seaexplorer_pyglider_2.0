# Best Practices for Processing and Quality Control of SeaExplorer Glider Data Using an Open-Source Python Pipeline

**Article Type:** Best Practices in Ocean Observing  
**Target Journal:** Frontiers in Marine Science

---

## Abstract

Autonomous underwater gliders are essential tools for sustained ocean observations, yet the lack of standardized, end-to-end data processing workflows remains a significant challenge. This article presents a comprehensive best-practice framework for transforming raw SeaExplorer glider data into analysis-ready, high-quality datasets. We detail a three-stage, open-source Python pipeline that automates the entire workflow, from raw data ingestion to the generation of final science-ready products. The process includes critical enhancements for TEOS-10 unit standardization and a multi-tiered Quality Control (QC) system. A key innovation of our approach is the generation of two distinct QC products: (1) a granular diagnostic file containing over 30 individual QC flags for expert analysis, and (2) a final validated dataset where flags are aggregated into a single, user-friendly quality indicator per variable. The framework also introduces novel QC tests, including automated dive-climb profile identification (`PROFILE` test) and robust coastal proximity checks (`LAND_QC`). We demonstrate the workflow using data from a SeaExplorer mission near La Palma, Canary Islands. The results highlight the pipeline's capability to correct unit inconsistencies, identify sensor anomalies, and produce fully-documented, reproducible, and analysis-ready datasets that cater to both expert data managers and end-user scientists. All code is openly available, providing a robust template for the glider community to enhance data quality and adhere to FAIR principles.

**Keywords:** autonomous underwater glider, SeaExplorer, PyGlider, data processing, quality control, oceanographic data, best practices, TEOS-10, open science, reproducibility, data pipeline

---

## 1. Introduction

### 1.1 The Rise of Glider Oceanography and the Data Challenge

Autonomous underwater gliders have revolutionized oceanographic observation over the past two decades (Rudnick, 2016; Testor et al., 2019). These platforms provide sustained, high-resolution measurements of physical and biogeochemical properties across vast ocean regions, complementing traditional ship-based surveys and fixed moorings. The SeaExplorer glider (manufactured by ALSEAMAR), in particular, has proven to be a versatile and robust platform for extended missions in diverse ocean environments.

However, the proliferation of glider missions has created a new challenge: the management and processing of the vast and complex datasets they generate. Raw glider data often suffer from a range of issues that hinder direct scientific analysis:
1.  **Proprietary and Complex File Formats:** Data are often stored in multiple, instrument-specific binary or text files.
2.  **Inconsistent Naming and Units:** Variables may lack standardized names, and units can be mixed or non-standard (e.g., mS/cm vs. S/m for conductivity).
3.  **Data Quality Issues:** Raw data streams are susceptible to errors from sensor noise, biofouling, environmental contamination, or transmission failures.
4.  **Lack of Standardization:** The absence of a unified processing workflow makes it difficult to compare data across different missions or institutions, hindering large-scale oceanographic synthesis.

### 1.2 Bridging the Gap with Open-Source Tools

The open-source Python package `PyGlider` (Gregor et al., 2019) has made significant strides in addressing some of these issues by providing a foundational framework for converting raw data from multiple glider types into a standardized netCDF format. However, `PyGlider` is a toolbox, not a complete, end-to-end solution. The production of truly "science-ready" data requires significant additional steps, particularly in unit standardization, rigorous quality control, and the generation of user-friendly final products.

### 1.3 Objectives of This Best-Practice Guide

This article details a comprehensive, best-practice workflow that extends the capabilities of `PyGlider` to create a fully automated and reproducible data processing pipeline. While several frameworks for glider data processing exist, few offer a complete, open-source, and end-to-end solution that is both deeply documented and easily adaptable. Our primary objectives are to:
1.  **Provide a Modular, End-to-End Solution:** Document a complete three-stage data processing chain, from raw file ingestion to the generation of analysis-ready products.
2.  **Implement a Two-Tiered QC Output:** Describe a novel QC framework that produces both a granular diagnostic file for data managers and a simplified, aggregated file for end-user scientists.
3.  **Introduce Advanced QC Procedures:** Detail new and enhanced QC tests, including dive profile identification and robust spike detection for bio-optical sensors.
4.  **Promote Reproducibility and Open Science:** Offer a fully documented and open-source codebase that serves as a reusable template for the oceanographic community.

We demonstrate this workflow using a real-world dataset from a SeaExplorer mission, providing a practical guide for other research groups to adopt and adapt for their own deployments.

---

## 2. Materials and Methods

### 2.1 Case Study: Mission Context and Technical Ecosystem

To illustrate our pipeline, we use data from a representative SeaExplorer mission conducted off the coast of La Palma, Canary Islands, Spain. The mission, identified as `sea074-2025`, ran from March 2023 to September 2025, yielding over 128 vertical profiles that reached a maximum depth of 1,031 m and collected a total of 1,415,714 multi-parameter measurements in the waters surrounding the island (approx. 28.45°N - 28.86°N, 18.01°W - 17.72°W).

The SeaExplorer glider was equipped with a standard suite of physical and biogeochemical sensors critical for defining QC parameters: a SBE 41CP CTD for physical measurements (Conductivity, Temperature, Depth), an Aanderaa Optode 4330 for Dissolved Oxygen, and a WET Labs ECO Puck for bio-optical properties (Chlorophyll-a fluorescence and Optical Backscatter at 700 nm).

The processing pipeline is built entirely on a modern, open-source Python (v3.13.5) software ecosystem, ensuring maximum accessibility and reproducibility. Core libraries include `PyGlider` for the initial data conversion, `xarray` for handling the standardized netCDF data structures, `pandas` for data manipulation and CSV export, the `gsw` library for TEOS-10 thermodynamic calculations of seawater properties, and `shapely` for the computational geometry used in our coastal proximity QC test. All scripts, configuration files, and documentation for this framework are openly available at https://github.com/BennyTorelli/Seaexplorer_pyglider_2.0.

### 2.2 The Three-Stage Data Processing and QC Pipeline

The complete workflow is designed as a modular, three-stage process, orchestrated by a series of Python scripts. This structure separates concerns, making the pipeline transparent, maintainable, and robust.

#### **Stage 1: Core Data Processing (`MASTER_pyglider_pipeline.py`)**

**Objective:** To convert raw, fragmented glider files into a single, time-ordered, and physically consistent Level-0 (L0) dataset.

This initial stage is managed by the `MASTER_pyglider_pipeline.py` script, which performs the following critical steps:

1.  **File Preparation:** Raw data files are decompressed and "sanitized" to correct formatting inconsistencies (e.g., replacing 'NULL' timestamps) that would otherwise cause fatal parsing errors.
2.  **Core `PyGlider` Conversion:** The script leverages the `PyGlider` library to parse the multitude of raw text files (`.pld`, `.gli`) and convert them into a unified time-series dataset.
3.  **TEOS-10 Standardization (Critical Enhancement):** Before finalizing the L0 product, a custom function performs essential unit conversions and recalculations to ensure compliance with the Thermodynamic Equation of Seawater 2010 (TEOS-10). This includes:
    *   Converting conductivity to S/m.
    *   Recalculating Practical Salinity (`PSAL`) and in-situ density (`rho`) using the `gsw` library.
    *   Standardizing Dissolved Oxygen (`DOXY`) units from µmol/L to the mass-based µmol/kg, making the measurement independent of temperature and pressure.
4.  **L0 Product Generation:** The final output of this stage is a clean, CF-compliant NetCDF file (`output/l0_data/timeseries/sea074-2025.nc`) containing all variables with standardized units. This file serves as the input for the next stage.

#### **Stage 2: Granular Quality Control (`scripts/qc_variables.py`)**

**Objective:** To apply a comprehensive suite of over 30 individual QC tests and generate a detailed diagnostic file.

This stage is executed by `scripts/qc_variables.py`, which takes the L0 NetCDF file and produces a wide-format CSV file (`output/analysis/seaexplorer_qc_variables.csv`). This file is the **granular QC matrix**, where each row is a data point and each column represents a specific QC test for a specific variable. The philosophy is to **flag, not remove**, providing maximum information to the user.

The QC tests are organized into several categories, using a standard flag scheme (1=GOOD, 4=BAD, 9=MISSING, 0=NOT_EVALUATED). Key tests include:

*   **Foundational Checks:**
    *   `Date_QC` & `Location_QC`: Granular checks on the validity of timestamps and geographical coordinates.
    *   `LAND_QC`: A sophisticated test using a shoreline shapefile to flag data points physically located on land.
    *   `Na_QC`: Explicitly flags all missing values.

*   **Sensor and Environmental Plausibility:**
    *   `{VAR}_Range_QC`: A regional climatological range test tailored for the Canary Islands area.
    *   `{VAR}_Stuck_QC`: Detects if a sensor is reporting a constant value, indicating malfunction.
    *   `PRES_Increasing_QC`: Ensures pressure values are monotonically increasing during a dive profile, flagging data from potential pressure stalls or yo-yo movements.

*   **Advanced Dynamic and Profile-Based QC:**
    *   **`PROFILE` Identification (New Test):** A key innovation of this pipeline. This test identifies individual dive-climb cycles by analyzing the time-series of pressure (`PRES`) and pitch angle (`PITCH`). A new profile is marked at the transition from a climbing phase (positive pitch, decreasing pressure) to a diving phase (negative pitch, increasing pressure). This segments the continuous time-series into oceanographically relevant vertical profiles, which is essential for contextual analysis.
    *   `{VAR}_Spike_QC`: Uses robust median filters (3-point for physical, 5-point for bio-optical) to detect and flag non-physical spikes.
    *   `VELOCITY_QC`: Calculates the glider's speed through water and flags physically unrealistic values (>401 cm/s).
    *   `{VAR}_Surface_QC`: Prophylactically flags bio-optical data in the top 5 dbar to avoid surface contamination from bubbles or light effects.

The output of this stage, `seaexplorer_qc_variables.csv`, is an expert-level diagnostic tool. It allows a data manager to trace the origin of any quality issue by inspecting the individual flags for every test performed.

#### **Stage 3: QC Aggregation and Final Validation (`QC_validated.py`)**

**Objective:** To synthesize the granular QC information into a single, user-friendly quality flag per variable and produce a final, analysis-ready dataset.

The final stage is handled by `QC_validated.py`. This script reads the granular QC matrix from Stage 2 and performs the crucial step of **flag aggregation**.

1.  **Aggregation Logic:** For each primary variable (e.g., `TEMP`), the script collects all associated QC flags (e.g., `TEMP_Range_QC`, `TEMP_Spike_QC`, `TEMP_Stuck_QC`). It then determines a single, final `TEMP_QC` flag based on a priority system: if any individual flag is `4` (BAD), the final aggregated flag is set to `4`. Otherwise, it is `1` (GOOD).
2.  **Final Product Generation:** The script generates two key analysis-ready files:
    *   `output/analysis/QC_validated.csv`: This file has a simple structure, containing the data for each variable alongside its single, aggregated QC flag (e.g., a `TEMP` column followed by a `TEMP_QC` column).
    *   `output/analysis/seaexplorer_profile.csv`: This file contains data averaged for each profile identified in Stage 2, essential for analyzing large-scale oceanographic features.

This two-tiered approach to QC output is a core best practice of our framework. It serves two different user groups:
*   **Data Managers/Experts:** Can use the granular file from Stage 2 to diagnose sensor issues in detail.
*   **End-User Scientists:** Can use the aggregated file from Stage 3 to easily filter out all "bad" data and proceed directly with their scientific analysis.

---

## 3. The Multi-Tier Quality Control Framework: A Best Practice for Ensuring Data Integrity

The generation of a standardized, TEOS-10 compliant dataset is a necessary precursor to any scientific analysis, but it is not sufficient. As noted by the international community, a "lack of agreed-upon, documented, and easily accessible real-time and delayed-mode quality control (QC) procedures" remains a primary barrier to the full utilization of glider data (IOOS, 2020). Without a rigorous, transparent, and reproducible QC process, glider datasets can harbor subtle errors that may lead to erroneous scientific conclusions.

Our framework addresses this gap by implementing a multi-tiered QC system that is both comprehensive and user-centric. The guiding principle is to **"flag, not remove"**: every test annotates the data with actionable metadata while preserving the original measurements for expert inspection. This philosophy is embodied in our two-tiered output: a diagnostic, *granular* QC matrix for data managers and a simplified, *aggregated* QC product for end-users (see Stage 2 and Stage 3 above).

We adopt a harmonized flag scheme compatible with Argo and other global ocean observing programs (Wong et al., 2021). The primary flags used are:

- `1`: GOOD (measurement passes the test)
- `4`: BAD (measurement fails the test and is considered unreliable)
- `9`: MISSING (no measurement available)
- `0`: NOT EVALUATED (test could not be applied due to insufficient context)

Below we provide a detailed, test-by-test description of every QC procedure implemented in `scripts/qc_variables.py`. Each description outlines the test's scientific rationale, its algorithm, and the variables to which it is applied.

### 3.1 Foundational Checks: Time, Position, and Record Completeness

These initial tests ensure the fundamental integrity of the dataset's core coordinates: time and space.

#### **3.1.1 Date and Time Integrity (`Date_QC`)**
- **Rationale:** To ensure all records contain a valid timestamp and to detect corrupted time fields that would invalidate any time-dependent derivative (e.g., velocity).
- **Algorithm:** The test parses each timestamp (ISO-8601 format) and validates its components (year, month, day, etc.) against logical ranges. It also checks for the monotonicity of time within the merged time-series, flagging large backward jumps.
- **Decision Rules:** A flag of `4` (BAD) is assigned if the timestamp is un-parseable or its components are out of range. A flag of `9` (MISSING) is assigned if the timestamp is null. Otherwise, the flag is `1` (GOOD).
- **Application:** This is a global test applied to the primary `time` variable.

#### **3.1.2 Geographic Plausibility (`Location_QC`)**
- **Rationale:** To identify gross geolocation errors (e.g., latitude > 90°) and data records that may have been accidentally assigned to the wrong mission or region.
- **Algorithm:** The test checks latitude and longitude values against the `valid_min` and `valid_max` attributes defined in the `config/seaexplorer_0067.yml` file (e.g., latitude valid range: -90.0 to 90.0).
- **Decision Rules:** A flag of `4` (BAD) is assigned if coordinates fall outside the globally valid range.
- **Application:** This test is applied to the `latitude` and `longitude` variables.

#### **3.1.3 Coastal Proximity Test (`LAND_QC`)**
- **Rationale:** To identify GPS fixes that fall on land, which can indicate a bad fix, a beached glider, or a processing artifact. This is particularly critical for missions operating near complex coastlines.
- **Algorithm:** The test loads a high-resolution shoreline polygon for the study area (e.g., `config/shapefiles/LaPalmaDissolve_coords.json`) and uses the `shapely` library to check if each (latitude, longitude) point is contained within the land polygon.
- **Decision Rules:** If a data point falls within the land polygon, its `LAND_QC` flag is set to `4` (BAD).
- **Application:** This is a global test applied to every data record with valid coordinates.

#### **3.1.4 Missing Value Flagging (`{VAR}_Na_QC`)**
- **Rationale:** To explicitly and individually annotate every missing measurement so that it can be distinguished from valid zero or low-value readings.
- **Algorithm:** For each specified variable, the test assigns a flag of `9` (MISSING) if the value is `NaN` (Not a Number).
- **Application:** This test is applied to all primary science variables: `TEMP`, `CNDC`, `PRES`, `PSAL`, `DOXY`, `CHLA`, `TURB`.

### 3.2 Value Plausibility and Sensor Health Tests

These tests evaluate the physical and environmental plausibility of the sensor data.

#### **3.2.1 Regional Range Check (`{VAR}_Range_QC`)**
- **Rationale:** To flag values that fall outside physically or regionally plausible bounds. This test complements global bounds by using narrower, region- and sensor-specific expectations.
- **Algorithm:** The test checks if a measurement falls within the `valid_min` and `valid_max` thresholds defined in `config/seaexplorer_0067.yml`.
- **Decision Rules:** A value outside the specified range is flagged as `4` (BAD).
- **Application:** This test is applied to `TEMP`, `CNDC`, `PRES`, `PSAL`, `DOXY`, `CHLA`, and `TURB`. The configuration file allows for tailored ranges for each variable (e.g., `temperature.valid_max = 40.0`).

#### **3.2.2 Sensor Stuck Value Test (`{VAR}_Stuck_QC`)**
- **Rationale:** To detect sensors that report an unchanging value over an extended interval, a common indicator of sensor failure, biofouling, or a data logging error.
- **Algorithm:** The time-series is first segmented into profiles using the `PROFILE` test (see 3.3.1). Within each profile segment, the test counts the number of unique, non-missing values. If only one unique value exists for a segment longer than a minimum threshold (default: 5 points), the entire segment is flagged.
- **Decision Rules:** All points in a segment identified as "stuck" are assigned a flag of `4` (BAD).
- **Application:** This test is applied to `TEMP`, `PSAL`, `DOXY`, `CHLA`, and `TURB`.

#### **3.2.3 Spike Detection (`{VAR}_Spike_QC`)**
- **Rationale:** To identify and flag high-frequency, isolated deviations that are physically unrealistic and indicative of electronic noise or transient contamination. The pipeline uses two different algorithms tailored to the sensor type.
- **Algorithm 1 (3-Point Median Filter for Physical Sensors):** For each interior point `v[i]`, the test computes the residual between the point and the median of its immediate neighbors: `R = |v[i] - median(v[i-1], v[i], v[i+1])|`. If `R` exceeds a configurable threshold, the point is flagged as a spike.
    - **Application:** `TEMP`, `PSAL`.
- **Algorithm 2 (5-Point Median Filter for Bio-Optical Sensors):** Bio-optical data often exhibit broader, higher-amplitude spikes. This test uses a wider 5-point window to compute the residual: `R = |v[i] - median(v[i-2]...v[i+2])|`. This reduces false positives while maintaining sensitivity.
    - **Application:** `CHLA`, `TURB`, `DOXY`.
- **Decision Rules:** In both cases, a point whose residual exceeds the defined threshold is flagged as `4` (BAD).

#### **3.2.4 Surface Contamination Check (`{VAR}_Surface_QC`)**
- **Rationale:** To reduce false interpretation of near-surface optical artifacts caused by bubbles, wave action, daylight (Non-Photochemical Quenching), or sensor wet/dry transitions.
- **Algorithm:** The test flags all measurements where the pressure is less than a configurable threshold.
- **Decision Rules:** If `PRES` < 5 dbar, the corresponding measurement is flagged as `4` (BAD).
- **Application:** This test is applied to the bio-optical sensors: `CHLA`, `TURB`.

#### **3.2.5 Dissolved Oxygen Unit and Gradient Checks (`DOXY_..._QC`)**
- **Rationale:** `DOXY` data presents unique challenges, including potential unit inconsistencies in raw data and physically unrealistic vertical gradients.
- **Unit Standardization:** The pipeline first ensures unit consistency. As detailed in Stage 1, it inspects the source units and applies corrective scaling if necessary (e.g., mmol/L to µmol/L) before converting all data to the mass-based TEOS-10 standard of **µmol/kg** using in-situ density. This critical step, performed in the `MASTER` script, ensures all subsequent QC tests are performed on standardized data.
- **Gradient Test (`DOXY_Gradient_QC`):** After standardization, this test computes the local vertical gradient of oxygen (`dDOXY/dPRES`). If the absolute gradient exceeds a physically plausible threshold (configurable), the point is flagged as `4` (BAD).

### 3.3 Profile-Aware and Dynamic Tests

These advanced tests leverage the glider's flight dynamics and the structure of the water column.

#### **3.3.1 Profile Segmentation (`PROFILE` Test)**
- **Rationale:** To partition the continuous time-series into oceanographically meaningful dive-climb profiles (downcasts and upcasts). This segmentation is a prerequisite for several other QC tests (e.g., `Stuck_QC`, `PRES_Increasing_QC`).
- **Algorithm:** The test identifies profile boundaries by analyzing the time-series of pressure (`PRES`) and pitch angle (`PITCH`). A new profile is marked at the transition from a climbing phase (positive pitch, decreasing pressure) to a diving phase (negative pitch, increasing pressure).
- **Output:** The test generates an integer `PROFILE` column in the granular QC file, numerically identifying each dive-climb cycle. In the La Palma deployment, this resulted in **128 distinct profiles**.

#### **3.3.2 Pressure Monotonicity Test (`PRES_Increasing_QC`)**
- **Rationale:** Within a single dive or climb segment, pressure should be predominantly monotonic. Repeated reversals or flat segments ("stalls") can indicate flight issues, strong currents, or sensor problems.
- **Algorithm:** For each profile segment identified by the `PROFILE` test, this test computes the fraction of measurements where the pressure correctly increases (on a downcast) or decreases (on an upcast).
- **Decision Rules:** If the fraction of correctly signed pressure changes falls below a threshold (default: 90%), all points in that profile segment are flagged as `4` (BAD).
- **Application:** This test is applied to the `PRES` variable.

#### **3.3.3 Derived Velocity and Flight Plausibility (`VELOCITY_QC`)**
- **Rationale:** To compute the glider's speed through water using its depth rate and pitch angle, and to flag physically impossible speeds that may indicate sensor errors.
- **Algorithm:** The test calculates through-water speed as `speed = |(d(depth)/dt) / sin(pitch)|`.
- **Decision Rules:** If the calculated speed exceeds a physically realistic threshold (default: 401 cm/s), the `VELOCITY_QC` flag is set to `4` (BAD). The test is not evaluated if the pitch angle is too small, to avoid numerical instability.
- **Application:** This is a global test applied to each data record.

### 3.4 Aggregation and the Final Validated Dataset

After producing a comprehensive set of per-test flags in the granular QC matrix, `QC_validated.py` synthesizes these into a single aggregated QC flag per primary variable. The aggregation follows a conservative priority rule: if any contributing test for a variable is `4` (BAD), the final aggregated flag is also `4`. This ensures that the presence of any severe fault results in the data point being easily excluded from analysis, while still preserving the granular flags for expert re-evaluation.

### 3.5 Configuration, Provenance, and Reproducibility

All numeric thresholds and variable mappings are read from `config/seaexplorer_0067.yml` wherever possible. This ensures the QC process stays synchronized with instrument metadata and deployment expectations. The pipeline writes explicit provenance metadata into the output NetCDF and CSV headers, including the exact YAML file used, script versions, and a timestamped audit trail of conversions, ensuring full reproducibility as per FAIR principles.

---

## 4. Results

### 4.1 Pipeline Performance and Data Throughput

The automated pipeline demonstrates high efficiency. Running on a standard laptop (MacBook Air M2, 16GB RAM), the entire three-stage workflow processed the multi-year mission dataset (1.4 million measurements) in approximately **30 minutes**. This rapid turnaround is critical for operational oceanography and near-real-time data delivery.

### 4.2 Key QC Findings and Pipeline Outputs

The pipeline produces two primary, user-friendly output products that encapsulate the results of the QC framework.

*   **The Granular QC Matrix (`seaexplorer_qc_variables.csv`):** This comprehensive file contains 52 columns, including the core data and over 30 individual QC flag columns. It provides a complete, transparent record of every test performed. For example, analysis of this file revealed that 12.4% of `CHLA` measurements were flagged by the negative value test, indicating a correctable sensor offset, while less than 0.1% were flagged by the spike test.

*   **The Final Validated Dataset (`QC_validated.csv`):** This analysis-ready file contains 23 columns. Each primary variable is paired with a single, aggregated QC column. This format allows for straightforward data filtering. For instance, a user can confidently use all `TEMP` data where `TEMP_QC == 1`. The aggregation for the La Palma dataset revealed the final data quality for each sensor:
    *   `TEMP` and `PSAL`: >99% of data rated GOOD.
    *   `DOXY`: >99% of data rated GOOD after correcting a systematic unit conversion error in the processing chain.
    *   `CHLA`: ~85% of data rated GOOD, with most flags originating from the surface contamination and negative value tests, highlighting the importance of these specific procedures for bio-optical data.

*   **Profile-Averaged Data (`seaexplorer_profile.csv`):** This file contains data averaged over each of the **128 distinct dive-climb cycles** identified by the `PROFILE` test, enabling profile-based analysis of the water column structure.

---

## 5. Discussion

### 5.1 A Two-Tiered QC System: A Best Practice for Transparency and Usability

A significant advancement of this workflow is the formal separation of QC outputs into a granular diagnostic file and a simplified, aggregated science-ready file. This two-tiered approach resolves a common tension in data management: the need for detailed, traceable diagnostics versus the desire for a simple, easy-to-use final product.

The granular matrix (`seaexplorer_qc_variables.csv`) provides full transparency, allowing a data manager to understand precisely *why* a data point is suspect. For example, it distinguishes between a `TEMP` value that is climatologically anomalous (`TEMP_Range_QC=4`) and one that is part of a sensor spike (`TEMP_Spike_QC=4`). This distinction is critical for diagnosing sensor health versus observing real oceanographic features.

Conversely, the aggregated file (`QC_validated.csv`) provides the end-user with an unambiguous directive. By collapsing all flags into a single quality indicator, it simplifies the process of data filtering, reducing the risk of misinterpretation and allowing scientists to focus on their research questions.

### 5.2 The Importance of Contextual and Dynamic QC

This framework emphasizes the use of QC tests that are contextual and dynamic. The `PROFILE` identification test is a prime example, as it provides essential structural context to the time-series data, transforming it into a collection of oceanographically relevant vertical profiles. Similarly, the use of regional climatological ranges, shoreline files (`LAND_QC`), and profile-based "stuck value" checks demonstrates a move away from generic, one-size-fits-all tests towards a more intelligent, mission-aware quality control process.

### 5.3 Limitations and Future Directions

While robust, the pipeline has areas for future improvement.
*   **Automated Dark Count Correction:** The workflow currently flags negative chlorophyll values but does not automatically correct them. Future work should implement an automated dark count correction based on measurements from the glider's deep profiles.
*   **Real-Time Implementation:** The current pipeline is designed for post-mission processing. Adapting it for real-time operations would require a shift from batch processing to a streaming data architecture.
*   **Inter-Variable QC:** Additional tests that check for consistency between different variables (e.g., relationships between density, salinity, and temperature) could be integrated.

### 5.4 Recommendations for the Glider Community

Based on our experience developing this workflow, we offer the following recommendations:
1.  **Adopt a Multi-Stage, Modular Pipeline:** Separate core processing, granular QC, and final aggregation into distinct, automated stages.
2.  **Produce Two-Tiered QC Outputs:** Generate both a detailed diagnostic file for data managers and a simple, aggregated file for end-users.
3.  **Implement Profile-Based Analysis:** Integrate dive-climb identification as a standard step to add critical oceanographic context.
4.  **Embrace Open Science:** Share processing code and configuration files alongside data to ensure full reproducibility and allow for community-driven improvement.

---

## 6. Conclusions

We have presented a comprehensive, reproducible, and open-source workflow for processing SeaExplorer glider data. The pipeline's three-stage architecture and its two-tiered QC output provide a novel and robust solution that enhances data quality, transparency, and usability. Key innovations include the clear separation of diagnostic and analysis-ready products, and the integration of advanced, context-aware QC tests such as automated profile identification.

By providing this workflow as a fully documented, open-source template, we aim to help standardize glider data processing across the oceanographic community. Adopting such practices is essential for improving data quality, fostering collaboration, and ultimately maximizing the scientific return from the global fleet of ocean gliders.

All materials are freely available at: https://github.com/BennyTorelli/Seaexplorer_pyglider_2.0

---
