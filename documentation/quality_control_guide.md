# Guide to the Quality Control Framework for SeaExplorer Data

**Version:** 2.0  
**Last Updated:** November 6, 2025

---

## 1. Introduction to Quality Control (QC)

This document is a practical guide to the **13-tier Quality Control (QC) framework** applied to the SeaExplorer glider data. The goal of QC is to identify and flag questionable or erroneous data points, ensuring that scientific analyses are based on the most reliable measurements possible.

Our philosophy is to **flag, not remove**. This means we add information about data quality, but we never delete the original data. This allows you, the end-user, to make informed decisions about which data to use for your specific analysis.

The output of this process is a comprehensive CSV file (`seaexplorer_qc_variables.csv`) where every measurement has a corresponding set of QC flags.

---

## 2. Understanding the QC Flags

We use a simple flagging system based on the international Argo program. Every data point for each QC test receives one of the following flags:

| Flag | Name | Description | Interpretation & Action |
|------|------|-------------|----------------|
| **1** | **GOOD** | The data point has passed the test. | Data is considered reliable for this specific test. |
| **4** | **BAD** | The data point has failed the test. | Data is considered erroneous or highly suspect. **Action:** You should probably exclude this data point from your analysis. |
| **9** | **MISSING** | The data point was not recorded (`NaN`). | Represents a data gap. This is not an error, just an absence of information. |
| **0** | **NOT EVALUATED** | The test was not performed on this data point. | This usually happens at the edges of a profile or if other required data is missing. It doesn't mean the data is bad. |

---

## 3. The New Data Columns: PITCH and VELOCITY

In addition to the QC flags, the new pipeline version adds two important columns to the main data files, providing insight into the glider's flight dynamics.

### `PITCH`
*   **What it is:** The pitch angle of the glider in degrees. A positive pitch indicates the glider is pointing upwards (ascending), while a negative pitch indicates it is pointing downwards (descending).
*   **Why it's useful:** It tells you about the glider's orientation and is essential for calculating the glider's speed through the water.

### `VELOCITY`
*   **What it is:** The calculated speed of the glider through the water, in **centimeters per second (cm/s)**.
*   **How it's calculated:** It is derived from the change in `DEPTH` over `TIME`, corrected by the `PITCH` angle. The absolute value is taken, so it always represents a positive speed.
    $$ VELOCITY_{cm/s} = \left| \frac{\Delta \text{DEPTH} / \Delta \text{TIME}}{\sin(\text{PITCH}_{\text{radians}})} \right| \times 100 $$
*   **Why it's useful:** This variable is a powerful diagnostic tool. It can help identify issues with the glider's propulsion, unusual flight behavior, or strong underlying currents.

---

## 4. Detailed QC Test Descriptions

Here is a breakdown of every test applied to the data.

### 4.1 Tier 1: Granular Temporal Validity Test (`Date_QC`)

*   **Purpose:** To ensure every timestamp is valid and makes sense.
*   **How it works:** This is an enhanced test that checks each part of the timestamp (`YYYY-MM-DD HH:MM:SS`) individually. It verifies:
    *   The **Year** is reasonable (e.g., after 1990).
    *   The **Month** is between 1 and 12.
    *   The **Day** is valid for that specific month (e.g., no 31st of April, handles leap years).
    *   The **Hour**, **Minute**, and **Second** are all within their correct ranges.
*   **Flagging:**
    *   `1` (GOOD): The timestamp is perfectly valid.
    *   `4` (BAD): Any part of the timestamp is incorrect.

### 4.2 Tier 2: Geographic Position Tests

#### 4.2.1 Gross Location Test (`Location_QC`)

*   **Purpose:** To catch major GPS errors.
*   **How it works:** Checks that the glider's latitude and longitude are within a valid range for the mission's operational area (e.g., around the Canary Islands).
*   **Flagging:**
    *   `1` (GOOD): The position is within the expected box.
    *   `4` (BAD): The position is far outside the expected area.

#### 4.2.2 Coastal Proximity Test (`LAND_QC`)

*   **Purpose:** To check if the glider is reporting a position on land.
*   **How it works:** Compares every GPS coordinate with a high-resolution map of La Palma. If the point falls inside the land polygon, it's flagged as bad. This is very important for missions close to the coast.
*   **Flagging:**
    *   `1` (GOOD): The glider is in the ocean.
    *   `4` (BAD): The glider's reported position is on land.
    *   `9` (MISSING): GPS coordinates are not available for this point.

### 4.3 Tier 3: Missing Value Test (`{VAR}_Na_QC`)

*   **Purpose:** To clearly label all data gaps.
*   **How it works:** This simple test checks if a measurement for a variable is `NaN` (Not a Number). It applies to `TEMP`, `CNDC`, `DOXY`, `CHLA`, and `TURB`.
*   **Flagging:**
    *   `1` (GOOD): A valid number is present.
    *   `9` (MISSING): The value is `NaN`.

### 4.4 Tier 4: Sensor Physical Limits Test (`{VAR}_Sensor_QC`)

*   **Purpose:** To check if a sensor is reporting values that are physically impossible for it to measure.
*   **How it works:** Compares each measurement against the manufacturer's specified range. For example, the temperature sensor cannot read values below -5°C or above 42°C. A value outside this range indicates a serious sensor or data corruption issue.
*   **Flagging:**
    *   `1` (GOOD): The value is within the sensor's technical limits.
    *   `4` (BAD): The value is outside the sensor's limits.

### 4.5 Tier 5: Regional Range Test (`{VAR}_Range_QC`)

*   **Purpose:** To check if measurements are plausible for the specific ocean region (La Palma).
*   **How it works:** This test is more specific than the sensor limits test. It uses a narrower range based on historical data for the Canary Islands. For example, while the sensor *can* read 35°C, a temperature that high would be extremely unlikely in this region and would be flagged. This helps detect unusual oceanographic events or sensor drift.
*   **Flagging:**
    *   `1` (GOOD): The value is within the expected range for the La Palma region.
    *   `4` (BAD): The value is outside the regional range.

### 4.6 Tier 6: Spike Detection (Physical Variables)

*   **Purpose:** To find and flag isolated, sharp spikes in the physical sensor data (`TEMP`, `PSAL`, `DOXY`).
*   **How it works:** It looks at a small window of three data points. If the middle point is a dramatic outlier compared to its two neighbors, it's flagged as a spike. This is good for catching random electronic noise.
*   **Flagging:**
    *   `1` (GOOD): The point is consistent with its neighbors.
    *   `4` (BAD): The point is an anomalous spike.

### 4.7 Tier 7: Enhanced Spike Detection (Bio-Optical)

*   **Purpose:** A specialized test to find spikes in the noisier bio-optical sensors (`CHLA`, `TURB`).
*   **How it works:** This enhanced test uses a wider 5-point window and calculates the median. It's particularly good at finding sharp, negative spikes. A key feature is the **strict neighbor requirement**: if any of the neighbors in the window are `NaN`, the test is not performed on that point. This avoids incorrectly flagging points near data gaps.
*   **Flagging:**
    *   `1` (GOOD): The point is consistent with the local median.
    *   `4` (BAD): The point is a significant negative anomaly.
    *   `0` (NOT EVALUATED): There are not enough valid neighbors to perform the test.

### 4.8 Tier 8: Surface Contamination Test (Bio-Optical)

*   **Purpose:** To flag bio-optical data near the surface that might be contaminated.
*   **How it works:** Measurements from `CHLA` and `TURB` sensors in the top 5 meters (pressure <= 5 dbar) are often unreliable due to air bubbles, intense sunlight, or surface film. This test simply flags all bio-optical data in this shallow layer as bad.
*   **Flagging:**
    *   `1` (GOOD): The measurement was taken deeper than 5 meters.
    *   `4` (BAD): The measurement was taken at or shallower than 5 meters.

### 4.9 Tier 9: Vertical Gradient Test (DOXY)

*   **Purpose:** To detect unrealistically sharp changes in the vertical profile of dissolved oxygen.
*   **How it works:** This test checks if the difference in oxygen between one point and the next is too large to be physically plausible. It helps identify issues with sensor response time.
*   **Flagging:**
    *   `1` (GOOD): The change in oxygen is reasonable.
    *   `4` (BAD): The change in oxygen is too abrupt.

### 4.10 Tier 10: Pressure Integrity Tests

#### 4.10.1 Maximum Pressure Test (`PRES_Max_QC`)

*   **Purpose:** To ensure the glider did not exceed its maximum depth rating.
*   **How it works:** Flags any pressure reading greater than 1100 dbar (a safe margin for a 1000m-rated glider).
*   **Flagging:**
    *   `1` (GOOD): Pressure is within the safe limit.
    *   `4` (BAD): Pressure exceeds the safe limit.

#### 4.10.2 Pressure Increasing Test (`PRES_Increasing_QC`)

*   **Purpose:** To check if the pressure sensor is "stuck."
*   **How it works:** A moving glider should always see its pressure change. This test flags any point where the pressure is identical to the previous measurement.
*   **Flagging:**
    *   `1` (GOOD): The pressure value has changed.
    *   `4` (BAD): The pressure value is identical to the previous one.

### 4.11 Tier 11: Density Inversion Test

*   **Purpose:** To detect physically unstable water column conditions.
*   **How it works:** In the ocean, denser water should be below lighter water. If the data shows a significant amount of lighter water below denser water, it's called a density inversion. This usually points to a problem with the salinity measurement. This test flags the `TEMP`, `PSAL`, and `PRES` data points that cause such an inversion.
*   **Flagging:**
    *   `1` (GOOD): The water column is stable.
    *   `4` (BAD): A significant density inversion was detected.

### 4.12 Tier 12: Stuck Sensor Test

*   **Purpose:** To identify if a sensor is stuck on a single value for an entire profile.
*   **How it works:** This is a more advanced version of the pressure test. It looks at a whole ascent or descent. If a sensor (like `TEMP` or `CHLA`) reports the exact same number for the entire profile, it's clearly broken. The test flags all points in that profile for that variable as bad.
*   **Flagging:**
    *   `1` (GOOD): The sensor shows varying values throughout the profile.
    *   `4` (BAD): The sensor is stuck on one value for the whole profile.

### 4.13 Tier 13: Velocity Plausibility Test (`VELOCITY_QC`)

*   **Purpose:** To check if the calculated glider speed is physically realistic.
*   **How it works:** This test flags any `VELOCITY` value greater than **401 cm/s**. This threshold is based on the sum of the glider's maximum possible speed (~51 cm/s) and the speed of the strongest known ocean currents (~350 cm/s). A speed higher than this is physically impossible and indicates an error in the `DEPTH`, `TIME`, or `PITCH` values used to calculate it.
*   **Flagging:**
    *   `1` (GOOD): The velocity is within a plausible range.
    *   `4` (BAD): The velocity is physically impossible.

---

## 5. Conclusion

This guide provides a clear overview of every quality check performed on the SeaExplorer data. By understanding what each flag means, you can confidently filter the dataset to use only the highest-quality data for your research.
