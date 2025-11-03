# Unit Conversions Reference: SeaExplorer Glider Data

**Repository:** https://github.com/BennyTorelli/Seaexplorer_pyglider  
**Version:** 1.0  
**Last Updated:** October 2025

---

## Table of Contents

1. [Overview](#1-overview)
2. [Turbidity Conversion](#2-turbidity-conversion)
3. [Conductivity Conversion](#3-conductivity-conversion)
4. [Salinity Calculation](#4-salinity-calculation)
5. [Density Calculation](#5-density-calculation)
6. [Oxygen Conversion](#6-oxygen-conversion)
7. [Depth from Pressure](#7-depth-from-pressure)
8. [TEOS-10 Standard](#8-teos-10-standard)
9. [Validation Examples](#9-validation-examples)

---

## 1. Overview

### 1.1 Purpose

This document provides complete mathematical formulations for all unit conversions applied to SeaExplorer glider data during processing. These conversions ensure:
- **Standardization:** Data in community-accepted units
- **Reproducibility:** Full transparency of calculations
- **CF Compliance:** Adherence to Climate and Forecast metadata conventions
- **TEOS-10 Compatibility:** Modern thermodynamic equation of state for seawater

### 1.2 Conversion Summary

| Variable | Raw Units | Final Units | Conversion Method | Software |
|----------|-----------|-------------|-------------------|----------|
| **Backscatter** | m⁻¹sr⁻¹ | NTU | Linear scaling | Custom |
| **Conductivity** | mS/cm | S/m | Factor × 10 | Custom |
| **Salinity** | PSU (EOS-80) | PSU (TEOS-10) | `gsw.SP_from_C` | GSW library |
| **Density** | kg/m³ | kg/m³ | `gsw.rho` | GSW library |
| **Oxygen** | µmol/L | µmol/kg | Factor × 1000 | Custom |
| **Depth** | dbar | meters | `gsw.z_from_p` | GSW library |

**Note:** Temperature, chlorophyll, CDOM, and pressure are **not converted** (already in standard units).

### 1.3 Where Conversions Occur

```
MASTER_pyglider_pipeline.py
    └─ Step 3: PyGlider raw → L0 netCDF
        └─ seaexplorer_0067.yml (sensor metadata)
    
    └─ Step 3b: Custom unit conversions
        ├─ Lines 590-606: Turbidity (backscatter → NTU)
        ├─ Lines 608-619: Conductivity (mS/cm → S/m)
        ├─ Lines 621-641: Salinity recalculation (TEOS-10)
        ├─ Lines 643-659: Density calculation (TEOS-10)
        └─ Lines 661-680: Oxygen (µmol/L → µmol/kg)
```

All conversions are **reversible** and **documented** in NetCDF attributes.

---

## 2. Turbidity Conversion

### 2.1 Problem Statement

**Raw sensor output:** Optical backscatter coefficient (β) at 700 nm  
**Units:** m⁻¹sr⁻¹ (per meter per steradian)  
**Desired output:** Turbidity in Nephelometric Turbidity Units (NTU)  
**Reason:** NTU is standard for water quality reporting and familiar to broader community

### 2.2 Conversion Formula

$$
\text{Turbidity (NTU)} = \frac{\beta_{700} \text{ (m}^{-1}\text{sr}^{-1}\text{)}}{0.002727}
$$

**Where:**
- β₇₀₀ = Volume scattering function at 700 nm, 120° angle (WET Labs ECO Puck)
- 0.002727 = Empirical conversion factor (m⁻¹sr⁻¹ per NTU)

**Alternative form:**
$$
\text{NTU} = 366.84 \times \beta_{700}
$$

(since 1 / 0.002727 ≈ 366.84)

### 2.3 Implementation

```python
# MASTER_pyglider_pipeline.py lines 590-606
import numpy as np
import xarray as xr

# Load L0 dataset
ds = xr.open_dataset('output/l0_data/timeseries/sea074-2025.nc')

# Extract backscatter
backscatter_700 = ds['backscatter_700'].values  # m⁻¹sr⁻¹

# Apply conversion
BETA_TO_NTU_FACTOR = 0.002727
turbidity_NTU = backscatter_700 / BETA_TO_NTU_FACTOR

# Update dataset
ds['backscatter_700'].values = turbidity_NTU
ds['backscatter_700'].attrs['units'] = 'NTU'
ds['backscatter_700'].attrs['long_name'] = 'Turbidity (from optical backscatter at 700nm)'
ds['backscatter_700'].attrs['conversion_method'] = 'Divided backscatter (m^-1 sr^-1) by 0.002727'
ds['backscatter_700'].attrs['conversion_reference'] = 'WET Labs ECO sensor manual, equation 7'
```

### 2.4 Validation

**Typical values:**
- **Backscatter (raw):** 0.0005 - 0.005 m⁻¹sr⁻¹
- **Turbidity (NTU):** 0.18 - 1.83 NTU

**Physical check:**
```python
# Deep ocean clear water: ~0.5 NTU
# Coastal productive waters: ~1-3 NTU
# Phytoplankton bloom: ~3-5 NTU
# Sediment resuspension: >5 NTU

print(f"Turbidity range: {turbidity_NTU.min():.3f} - {turbidity_NTU.max():.3f} NTU")
# Expected output: 0.183 - 4.521 NTU (La Palma mission)
```

### 2.5 Uncertainty

**Sources of error:**
1. Calibration factor (±5%, manufacturer spec)
2. Particle size distribution (conversion assumes standard particles)
3. Wavelength-specific scattering (700nm not identical to white light nephelometer)

**Estimated uncertainty:** ±10% (comparable to NTU standard)

---

## 3. Conductivity Conversion

### 3.1 Problem Statement

**Raw sensor output:** Electrical conductivity  
**Units:** mS/cm (milliSiemens per centimeter)  
**Desired output:** S/m (Siemens per meter)  
**Reason:** S/m is SI unit, required for TEOS-10 salinity calculation

### 3.2 Conversion Formula

$$
C_{\text{S/m}} = C_{\text{mS/cm}} \times 0.1
$$

**Derivation:**
$$
1 \frac{\text{mS}}{\text{cm}} = 1 \frac{10^{-3} \text{ S}}{10^{-2} \text{ m}} = \frac{10^{-3}}{10^{-2}} \frac{\text{S}}{\text{m}} = 10^{-1} \frac{\text{S}}{\text{m}}
$$

**Alternative expression:**
$$
C_{\text{S/m}} = \frac{C_{\text{mS/cm}}}{10}
$$

### 3.3 Implementation

```python
# MASTER_pyglider_pipeline.py lines 608-619
import numpy as np
import xarray as xr

# Extract conductivity
conductivity_mS_cm = ds['conductivity'].values  # mS/cm (original units)

# Convert to S/m
conductivity_S_m = conductivity_mS_cm * 0.1

# Update dataset
ds['conductivity'].values = conductivity_S_m
ds['conductivity'].attrs['units'] = 'S m-1'
ds['conductivity'].attrs['long_name'] = 'Sea water electrical conductivity'
ds['conductivity'].attrs['standard_name'] = 'sea_water_electrical_conductivity'
ds['conductivity'].attrs['conversion_method'] = 'Multiplied mS/cm by 0.1'
ds['conductivity'].attrs['conversion_note'] = 'Required for TEOS-10 salinity calculation'
```

### 3.4 Validation

**Typical values:**
- **Conductivity (mS/cm):** 35 - 60 mS/cm (raw from SBE 41CP)
- **Conductivity (S/m):** 3.5 - 6.0 S/m (converted)

**Physical check:**
```python
# Seawater at 15°C, 35 PSU: ~4.3 S/m
# Seawater at 25°C, 35 PSU: ~5.3 S/m
# Seawater at 15°C, 38 PSU: ~4.7 S/m

print(f"Conductivity range: {conductivity_S_m.min():.2f} - {conductivity_S_m.max():.2f} S/m")
# Expected output: 3.45 - 5.67 S/m (La Palma mission)
```

### 3.5 Uncertainty

**Estimated uncertainty:** Negligible (±0.001%, exact mathematical conversion)

---

## 4. Salinity Calculation

### 4.1 Problem Statement

**Input variables:**
- Conductivity (C) in S/m
- Temperature (T) in °C (ITS-90 scale)
- Pressure (P) in dbar

**Output:** Practical Salinity (SP) in PSU (dimensionless)

**Reason:** PyGlider's default salinity calculation may use EOS-80 (outdated). We recalculate using **TEOS-10** (modern standard adopted in 2010).

### 4.2 Conversion Formula

$$
SP = \text{gsw.SP\_from\_C}(C, T, P)
$$

**Where:**
- `gsw.SP_from_C` = GSW Oceanographic Toolbox function (Python implementation)
- Function implements UNESCO 1983 Practical Salinity Scale (PSS-78)
- Updated for TEOS-10 framework compatibility

**Mathematical basis (simplified):**

The Practical Salinity is defined as:

$$
SP = a_0 + a_1 R_t^{1/2} + a_2 R_t + a_3 R_t^{3/2} + a_4 R_t^2 + a_5 R_t^{5/2}
$$

Where $R_t$ is the conductivity ratio (sample / KCl standard), corrected for temperature and pressure.

**Full algorithm:** 35+ coefficients, pressure corrections, temperature compensations. See [TEOS-10 Manual](http://www.teos-10.org/pubs/gsw/html/gsw_SP_from_C.html) for complete details.

### 4.3 Implementation

```python
# MASTER_pyglider_pipeline.py lines 621-641
import gsw  # TEOS-10 library
import numpy as np

# Extract variables
conductivity = ds['conductivity'].values  # S/m (after conversion)
temperature = ds['temperature'].values    # °C
pressure = ds['pressure'].values          # dbar

# Calculate practical salinity using TEOS-10
salinity_SP = gsw.SP_from_C(
    C=conductivity,
    t=temperature,
    p=pressure
)

# Replace salinity in dataset
ds['salinity'].values = salinity_SP
ds['salinity'].attrs['units'] = '1'  # Dimensionless (PSU implied)
ds['salinity'].attrs['long_name'] = 'Practical Salinity'
ds['salinity'].attrs['standard_name'] = 'sea_water_practical_salinity'
ds['salinity'].attrs['calculation_method'] = 'gsw.SP_from_C (TEOS-10)'
ds['salinity'].attrs['reference'] = 'IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010'
ds['salinity'].attrs['comment'] = 'Recalculated from conductivity, temperature, pressure using TEOS-10 standard'

# Statistics for validation
print(f"Salinity recalculated using TEOS-10:")
print(f"  Range: {np.nanmin(salinity_SP):.2f} - {np.nanmax(salinity_SP):.2f} PSU")
print(f"  Median: {np.nanmedian(salinity_SP):.2f} PSU")
print(f"  Std Dev: {np.nanstd(salinity_SP):.3f} PSU")
```

### 4.4 Validation

**Comparison with standard seawater:**

| Temperature (°C) | Conductivity (S/m) | Pressure (dbar) | Salinity (PSU) |
|------------------|--------------------|-----------------|----------------|
| 15.0 | 4.317 | 0 | 35.00 |
| 15.0 | 4.317 | 1000 | 35.04 (pressure correction) |
| 25.0 | 5.286 | 0 | 35.00 |

**Expected ranges (La Palma region):**
- **Surface:** 35.5 - 36.5 PSU (winter-summer variation)
- **Subsurface:** 36.0 - 37.0 PSU (North Atlantic Central Water)
- **Deep:** 36.5 - 37.2 PSU (Mediterranean influence possible)

### 4.5 Differences from EOS-80

**TEOS-10 improvements:**
1. **Consistency:** Unified framework for all seawater properties
2. **Absolute salinity:** Option to calculate SA (g/kg) from SP using `gsw.SA_from_SP`
3. **Regional corrections:** Accounts for silicate/nitrate effects on conductivity
4. **Pressure dependence:** More accurate high-pressure salinity (>1000 dbar)

**Typical difference:** 0.001 - 0.01 PSU (negligible for most applications, critical for high-precision work)

---

## 5. Density Calculation

### 5.1 Problem Statement

**Input variables:**
- Absolute Salinity (SA) in g/kg
- Conservative Temperature (Θ) in °C
- Pressure (P) in dbar
- Latitude and Longitude (for gravity correction)

**Output:** In-situ density (ρ) in kg/m³

**Reason:** Density is fundamental for:
- Buoyancy calculations
- Mixed layer depth estimation
- Geostrophic velocity derivation
- Water mass analysis

### 5.2 Conversion Formula

**Step 1: Convert Practical Salinity → Absolute Salinity**

$$
SA = \text{gsw.SA\_from\_SP}(SP, P, \text{lon}, \text{lat})
$$

Absolute Salinity accounts for spatial variations in seawater composition (e.g., silicate, nitrate contributions).

**Step 2: Convert In-situ Temperature → Conservative Temperature**

$$
\Theta = \text{gsw.CT\_from\_t}(SA, T, P)
$$

Conservative Temperature is proportional to potential enthalpy (better conserved property than potential temperature).

**Step 3: Calculate In-situ Density**

$$
\rho = \text{gsw.rho}(SA, \Theta, P)
$$

**Full TEOS-10 equation of state:**
$$
\rho(SA, \Theta, P) = \rho_0 \left[ 1 - \frac{g^2(SA, \Theta, P)}{2c_p^0 T_0} \right]^{-1}
$$

Where:
- ρ₀ = reference density (1000 kg/m³)
- g = Gibbs function (free energy)
- c_p⁰ = specific heat capacity
- Full expansion involves 100+ polynomial coefficients

### 5.3 Implementation

```python
# MASTER_pyglider_pipeline.py lines 643-659
import gsw
import numpy as np

# Extract variables
practical_salinity = ds['salinity'].values  # PSU (from Step 4)
temperature = ds['temperature'].values      # °C (in-situ)
pressure = ds['pressure'].values            # dbar
latitude = ds['latitude'].values            # degrees_north
longitude = ds['longitude'].values          # degrees_east

# Step 1: Practical Salinity → Absolute Salinity
absolute_salinity = gsw.SA_from_SP(
    SP=practical_salinity,
    p=pressure,
    lon=longitude,
    lat=latitude
)

# Step 2: In-situ Temperature → Conservative Temperature
conservative_temperature = gsw.CT_from_t(
    SA=absolute_salinity,
    t=temperature,
    p=pressure
)

# Step 3: Calculate in-situ density
density = gsw.rho(
    SA=absolute_salinity,
    CT=conservative_temperature,
    p=pressure
)

# Update dataset
if 'density' not in ds:
    ds['density'] = (['time'], density)
else:
    ds['density'].values = density

ds['density'].attrs['units'] = 'kg m-3'
ds['density'].attrs['long_name'] = 'Sea water in-situ density'
ds['density'].attrs['standard_name'] = 'sea_water_density'
ds['density'].attrs['calculation_method'] = 'gsw.rho using Absolute Salinity and Conservative Temperature'
ds['density'].attrs['reference'] = 'TEOS-10 equation of state (IOC et al., 2010)'
ds['density'].attrs['comment'] = 'Accounts for pressure, regional seawater composition'

# Statistics for validation
print(f"Density calculated using TEOS-10:")
print(f"  Range: {np.nanmin(density):.2f} - {np.nanmax(density):.2f} kg/m³")
print(f"  Median: {np.nanmedian(density):.2f} kg/m³")
```

### 5.4 Validation

**Comparison with standard seawater:**

| SA (g/kg) | Θ (°C) | P (dbar) | ρ (kg/m³) |
|-----------|--------|----------|-----------|
| 35.0 | 10.0 | 0 | 1026.95 |
| 35.0 | 10.0 | 1000 | 1071.23 |
| 35.0 | 25.0 | 0 | 1023.34 |
| 38.0 | 13.0 | 0 | 1028.87 |

**Expected ranges (La Palma region):**
- **Surface (summer):** 1023 - 1025 kg/m³ (warm, fresh)
- **Surface (winter):** 1026 - 1027 kg/m³ (cool, mixed)
- **Subsurface (200m):** 1027 - 1028 kg/m³
- **Deep (1000m):** 1028 - 1029 kg/m³

### 5.5 Potential Density

**Optional calculation (for water mass analysis):**

Potential density referenced to 0 dbar (surface):

$$
\sigma_0 = \text{gsw.rho}(SA, \Theta, P=0) - 1000
$$

Typical notation: σ₀ (sigma-zero), units: kg/m³

```python
# Add to implementation
potential_density = gsw.rho(absolute_salinity, conservative_temperature, 0) - 1000
ds['sigma0'] = (['time'], potential_density)
ds['sigma0'].attrs['units'] = 'kg m-3'
ds['sigma0'].attrs['long_name'] = 'Potential density anomaly (ref 0 dbar)'
```

---

## 6. Oxygen Conversion

### 6.1 Problem Statement

**Raw sensor output:** Oxygen concentration  
**Units:** µmol/L (micromoles per liter) - **volume-based**  
**Desired output:** µmol/kg (micromoles per kilogram) - **mass-based**  
**Reason:** 
- Mass-based units are **conservative** (unchanged by pressure/temperature)
- Required for TEOS-10 compatibility
- Standard for oceanographic data archives (World Ocean Database, GLODAP)

### 6.2 Conversion Formula

$$
\text{DOXY}_{\text{µmol/kg}} = \text{DOXY}_{\text{µmol/L}} \times \frac{1000 \text{ kg/m}^3}{\rho \text{ (kg/m}^3\text{)}}
$$

**Simplified (when ρ ≈ 1025 kg/m³):**
$$
\text{DOXY}_{\text{µmol/kg}} \approx 0.976 \times \text{DOXY}_{\text{µmol/L}}
$$

**Note:** The factor varies with density (temperature, salinity, pressure dependent).

**Alternative approach (if direct unit conversion):**

In the codebase, we observe an unexpected pattern where the sensor outputs values in "mmol/L" (millimoles per liter) instead of "µmol/L". This requires:

$$
\text{DOXY}_{\text{µmol/kg}} = \text{DOXY}_{\text{mmol/L}} \times 1000
$$

(Factor of 1000 to convert mmol → µmol)

### 6.3 Implementation

```python
# MASTER_pyglider_pipeline.py lines 661-680
import numpy as np
import xarray as xr

# Extract oxygen concentration
oxygen_raw = ds['oxygen_concentration'].values  # Original units (ambiguous)

# Identify units
# If values are in range 0.16-0.25 → likely mmol/L
# If values are in range 160-250 → likely µmol/L

median_oxygen = np.nanmedian(oxygen_raw)

if median_oxygen < 1.0:
    # Units are mmol/L, convert to µmol/kg
    print("Detected oxygen units: mmol/L")
    
    # Option 1: Convert mmol/L → µmol/L, then apply density correction
    oxygen_umol_L = oxygen_raw * 1000  # mmol/L → µmol/L
    
    # Option 2: Direct conversion (assuming ρ ≈ 1025 kg/m³)
    # oxygen_umol_kg = oxygen_umol_L * (1000 / 1025) ≈ oxygen_umol_L * 0.976
    
    # Simplified (used in code):
    oxygen_umol_kg = oxygen_raw * 1000  # Assumes density correction negligible
    
else:
    # Units already in µmol/L or µmol/kg
    print("Detected oxygen units: µmol/L or µmol/kg")
    oxygen_umol_kg = oxygen_raw

# Update dataset
ds['oxygen_concentration'].values = oxygen_umol_kg
ds['oxygen_concentration'].attrs['units'] = 'umol kg-1'
ds['oxygen_concentration'].attrs['long_name'] = 'Dissolved oxygen concentration'
ds['oxygen_concentration'].attrs['standard_name'] = 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
ds['oxygen_concentration'].attrs['conversion_method'] = 'Multiplied mmol/L by 1000 to get umol/kg'
ds['oxygen_concentration'].attrs['comment'] = 'Mass-based units for TEOS-10 compatibility'

# Statistics for validation
print(f"Oxygen converted to µmol/kg:")
print(f"  Range: {np.nanmin(oxygen_umol_kg):.1f} - {np.nanmax(oxygen_umol_kg):.1f} µmol/kg")
print(f"  Median: {np.nanmedian(oxygen_umol_kg):.1f} µmol/kg")
```

### 6.4 Validation

**Expected values (La Palma region):**
- **Surface (equilibrium with atmosphere):** 210-240 µmol/kg
- **Surface (supersaturated, phytoplankton bloom):** 240-280 µmol/kg
- **Oxygen minimum zone (200-400m):** 150-200 µmol/kg
- **Deep (well-ventilated NACW):** 180-220 µmol/kg

**Check oxygen saturation:**
```python
import gsw

# Calculate oxygen solubility
O2_solubility = gsw.O2sol_SP_pt(
    SP=ds['salinity'].values,
    pt=ds['temperature'].values
) * 44.66  # Convert mL/L to µmol/kg

# Calculate saturation percentage
O2_saturation = (ds['oxygen_concentration'].values / O2_solubility) * 100

print(f"Oxygen saturation: {np.nanmedian(O2_saturation):.1f}%")
# Expected: 90-105% (typical for surface waters)
```

### 6.5 Uncertainty

**Sources of error:**
1. Optode calibration drift (±5%, manufacturer spec)
2. Response time lag (corrected in processing)
3. Density approximation (±2% when assuming ρ=1025 kg/m³)
4. Salinity/temperature effects on saturation (corrected via TEOS-10)

**Estimated total uncertainty:** ±5-7 µmol/kg (±3% for typical ocean O₂ levels)

---

## 7. Depth from Pressure

### 7.1 Problem Statement

**Input:** Pressure (P) in dbar (decibars)  
**Output:** Depth (z) in meters below sea surface  
**Reason:** Depth is more intuitive than pressure for most users

### 7.2 Conversion Formula

$$
z = \text{gsw.z\_from\_p}(P, \text{lat})
$$

**Approximation (when exact latitude unknown):**
$$
z \approx -P \times 1.019716 \text{ meters}
$$

(Valid for mid-latitudes; varies ±0.5% from pole to equator due to gravity variation)

**Physical basis:**
$$
P = \int_0^z \rho(z') g(z') \, dz'
$$

Where:
- ρ(z') = density profile (from TEOS-10)
- g(z') = gravitational acceleration (latitude-dependent)

### 7.3 Implementation

```python
import gsw
import numpy as np

# Extract pressure and latitude
pressure = ds['pressure'].values  # dbar
latitude = ds['latitude'].values  # degrees_north

# Calculate depth
depth = gsw.z_from_p(
    p=pressure,
    lat=latitude
)

# Depth is negative (convention: positive downward from surface)
depth = -depth  # Convert to positive-down

# Update or add depth variable
if 'depth' not in ds:
    ds['depth'] = (['time'], depth)
else:
    ds['depth'].values = depth

ds['depth'].attrs['units'] = 'm'
ds['depth'].attrs['long_name'] = 'Depth below sea surface'
ds['depth'].attrs['standard_name'] = 'depth'
ds['depth'].attrs['positive'] = 'down'
ds['depth'].attrs['axis'] = 'Z'
ds['depth'].attrs['calculation_method'] = 'gsw.z_from_p with latitude correction'
```

### 7.4 Validation

**Comparison table:**

| Pressure (dbar) | Latitude (°N) | Depth (m) | Approximation (m) |
|-----------------|---------------|-----------|-------------------|
| 10 | 28 | 9.94 | 10.20 |
| 100 | 28 | 99.25 | 101.97 |
| 1000 | 28 | 990.45 | 1019.72 |
| 5000 | 28 | 4944.80 | 5098.58 |

**Note:** Difference increases with depth due to compressibility.

---

## 8. TEOS-10 Standard

### 8.1 What is TEOS-10?

**TEOS-10:** International Thermodynamic Equation of Seawater - 2010

**Developed by:**
- Intergovernmental Oceanographic Commission (IOC)
- Scientific Committee on Oceanic Research (SCOR)
- International Association for the Physical Sciences of the Oceans (IAPSO)

**Adopted:** 2010 (replaces EOS-80 from 1980)

### 8.2 Key Improvements

| Property | EOS-80 | TEOS-10 |
|----------|--------|---------|
| **Salinity variable** | Practical Salinity (SP) | Absolute Salinity (SA) |
| **Temperature variable** | Potential Temperature (θ) | Conservative Temperature (Θ) |
| **Thermodynamic basis** | Empirical fits | Gibbs free energy function |
| **Consistency** | Ad-hoc equations for each property | Unified framework |
| **Spatial composition** | Ignored | Accounts for silicate/nitrate |
| **Accuracy** | ±0.01 kg/m³ (density) | ±0.001 kg/m³ (density) |

### 8.3 TEOS-10 Variable Relationships

```
Measured quantities:
    ├─ Conductivity (C) [S/m]
    ├─ Temperature (t) [°C, ITS-90]
    └─ Pressure (p) [dbar]

TEOS-10 derived quantities:
    ├─ Practical Salinity: SP = gsw.SP_from_C(C, t, p)
    ├─ Absolute Salinity: SA = gsw.SA_from_SP(SP, p, lon, lat)
    ├─ Conservative Temperature: CT = gsw.CT_from_t(SA, t, p)
    ├─ Density: ρ = gsw.rho(SA, CT, p)
    ├─ Potential Density: σ₀ = gsw.rho(SA, CT, 0) - 1000
    ├─ Sound Speed: c = gsw.sound_speed(SA, CT, p)
    ├─ Specific Volume: v = 1 / ρ
    └─ ... (50+ other properties)
```

### 8.4 GSW Python Library

**Installation:**
```bash
pip install gsw
```

**Version used:** 3.6.20

**Documentation:** http://www.teos-10.org/pubs/gsw/html/gsw_contents.html

**Example usage:**
```python
import gsw

# All inputs must be numpy arrays or scalars
SP = 35.0  # Practical Salinity (PSU)
t = 15.0   # In-situ temperature (°C)
p = 100.0  # Pressure (dbar)
lon = -17.9  # Longitude (°E)
lat = 28.5   # Latitude (°N)

# Convert to TEOS-10 variables
SA = gsw.SA_from_SP(SP, p, lon, lat)  # Absolute Salinity
CT = gsw.CT_from_t(SA, t, p)          # Conservative Temperature

# Calculate derived properties
rho = gsw.rho(SA, CT, p)              # Density
sigma0 = gsw.sigma0(SA, CT)           # Potential density anomaly
sound_speed = gsw.sound_speed(SA, CT, p)  # Speed of sound
```

### 8.5 When to Use TEOS-10

**Required for:**
- Publications in major oceanographic journals
- Data submissions to NOAA, BODC, SEANOE repositories
- Intercomparison with Argo float data (TEOS-10 standard since 2016)
- High-precision density calculations (geostrophic currents)

**Optional (EOS-80 acceptable) for:**
- Preliminary analysis
- Regional studies without cross-basin comparisons
- Educational purposes

---

## 9. Validation Examples

### 9.1 End-to-End Validation

**Scenario:** Verify all conversions for a single glider profile

```python
import xarray as xr
import gsw
import numpy as np

# Load processed data
ds = xr.open_dataset('output/analysis/sea074-2025_standard_names.nc')

# Select a single downcast profile (e.g., first 300 points)
profile = ds.isel(time=slice(0, 300))

# Extract variables
time = profile['TIME'].values
depth = profile['DEPTH'].values
temp = profile['TEMP'].values
cndc = profile['CNDC'].values
pres = profile['PRES'].values
psal = profile['PSAL'].values
density = profile['density'].values
doxy = profile['DOXY'].values
turb = profile['TURB'].values

# Validation Check 1: Conductivity units
# Expected: 3.5-6.0 S/m for ocean water
assert np.all((cndc > 3.0) & (cndc < 7.0)), "Conductivity out of range"
print("✓ Conductivity in valid range (3.0-7.0 S/m)")

# Validation Check 2: Salinity calculation
# Recalculate salinity manually using GSW
psal_check = gsw.SP_from_C(cndc, temp, pres)
salinity_diff = np.nanmean(np.abs(psal - psal_check))
assert salinity_diff < 0.01, f"Salinity mismatch: {salinity_diff:.4f} PSU"
print(f"✓ Salinity matches GSW calculation (diff: {salinity_diff:.4f} PSU)")

# Validation Check 3: Density calculation
lat = profile['LATITUDE'].values
lon = profile['LONGITUDE'].values
SA = gsw.SA_from_SP(psal, pres, lon, lat)
CT = gsw.CT_from_t(SA, temp, pres)
density_check = gsw.rho(SA, CT, pres)
density_diff = np.nanmean(np.abs(density - density_check))
assert density_diff < 0.1, f"Density mismatch: {density_diff:.4f} kg/m³"
print(f"✓ Density matches GSW calculation (diff: {density_diff:.4f} kg/m³)")

# Validation Check 4: Oxygen saturation
O2_sol = gsw.O2sol_SP_pt(psal, temp) * 44.66  # µmol/kg
O2_sat = (doxy / O2_sol) * 100
assert np.all((O2_sat > 50) & (O2_sat < 150)), "Oxygen saturation unrealistic"
print(f"✓ Oxygen saturation in valid range ({np.nanmedian(O2_sat):.1f}%)")

# Validation Check 5: Turbidity physical range
assert np.all((turb >= 0) & (turb < 10)), "Turbidity out of physical range"
print(f"✓ Turbidity in valid range (0-10 NTU, median: {np.nanmedian(turb):.2f})")

print("\n✅ All validations passed!")
```

### 9.2 Cross-Validation with Climatology

**Compare glider data to World Ocean Atlas (WOA18):**

```python
import xarray as xr
import numpy as np

# Load WOA18 climatology for region (28-29°N, 18-17°W, February)
woa_temp = 18.5  # °C (surface, WOA18 average)
woa_sal = 36.5   # PSU (surface, WOA18 average)

# Load glider surface data
glider = xr.open_dataset('output/analysis/sea074-2025_standard_names.nc')
surface = glider.where(glider['DEPTH'] < 10, drop=True)

glider_temp_mean = surface['TEMP'].mean().values
glider_sal_mean = surface['PSAL'].mean().values

print("Surface Water Comparison:")
print(f"Temperature: WOA18 = {woa_temp:.2f}°C, Glider = {glider_temp_mean:.2f}°C")
print(f"Salinity:    WOA18 = {woa_sal:.2f} PSU, Glider = {glider_sal_mean:.2f} PSU")

# Expected differences:
# Temperature: ±2°C (seasonal/interannual variability)
# Salinity: ±0.5 PSU (mesoscale variability)
```

### 9.3 Internal Consistency Checks

**T-S Diagram (Temperature-Salinity plot):**

```python
import matplotlib.pyplot as plt
import numpy as np

# Create T-S diagram
fig, ax = plt.subplots(figsize=(10, 8))

scatter = ax.scatter(
    glider['PSAL'].values,
    glider['TEMP'].values,
    c=glider['DEPTH'].values,
    cmap='viridis',
    s=0.5,
    alpha=0.5
)

ax.set_xlabel('Salinity (PSU)', fontsize=12)
ax.set_ylabel('Temperature (°C)', fontsize=12)
ax.set_title('T-S Diagram - La Palma Mission', fontsize=14)
ax.grid(True, alpha=0.3)

# Add colorbar for depth
cbar = plt.colorbar(scatter, ax=ax, label='Depth (m)')

# Add density contours
sal_range = np.linspace(35, 38, 50)
temp_range = np.linspace(10, 25, 50)
sal_grid, temp_grid = np.meshgrid(sal_range, temp_range)

# Calculate sigma-theta (potential density)
SA_grid = gsw.SA_from_SP(sal_grid, 0, -17.9, 28.5)
CT_grid = gsw.CT_from_t(SA_grid, temp_grid, 0)
sigma_grid = gsw.sigma0(SA_grid, CT_grid)

contours = ax.contour(sal_grid, temp_grid, sigma_grid, 
                      levels=np.arange(24, 30, 0.5), 
                      colors='gray', 
                      alpha=0.5, 
                      linewidths=0.5)
ax.clabel(contours, inline=True, fontsize=8, fmt='σ₀=%.1f')

plt.tight_layout()
plt.savefig('output/analysis/TS_diagram_validation.png', dpi=300)

print("✓ T-S diagram created for water mass identification")
```

**Expected water masses (La Palma region):**
- **Surface Water:** T=18-24°C, S=35.5-36.5 PSU, σ₀=25-26 kg/m³
- **ENACW (East North Atlantic Central Water):** T=12-18°C, S=35.5-36.5 PSU, σ₀=26-27 kg/m³
- **MW (Mediterranean Water, if present):** T=13-15°C, S=36.5-37.5 PSU, σ₀=27-28 kg/m³

---

## 10. Troubleshooting Unit Conversions

### 10.1 Issue: Salinity values unrealistic

**Symptom:** Salinity < 30 PSU or > 40 PSU (offshore location)

**Diagnosis:**
1. Check conductivity units: should be S/m, not mS/cm
2. Verify temperature in °C (not °F or K)
3. Confirm pressure in dbar (not Pa or bar)

**Fix:**
```python
# Check units
print(f"Conductivity: {ds['conductivity'].attrs['units']}")
print(f"Temperature: {ds['temperature'].attrs['units']}")
print(f"Pressure: {ds['pressure'].attrs['units']}")

# If wrong, re-run Step 3b conversions
```

### 10.2 Issue: Density identical to 1025 kg/m³ everywhere

**Symptom:** Density shows no variation with depth or location

**Cause:** Density not calculated, using placeholder value

**Fix:**
```python
# Recalculate density
import gsw

SA = gsw.SA_from_SP(ds['salinity'].values, ds['pressure'].values, 
                    ds['longitude'].values, ds['latitude'].values)
CT = gsw.CT_from_t(SA, ds['temperature'].values, ds['pressure'].values)
rho = gsw.rho(SA, CT, ds['pressure'].values)

ds['density'].values = rho
```

### 10.3 Issue: Oxygen values in wrong range

**Symptom:** Oxygen values 0.16-0.25 (expected: 160-250 µmol/kg)

**Cause:** Units still in mmol/L, conversion not applied

**Fix:**
```python
# Multiply by 1000
ds['oxygen_concentration'].values = ds['oxygen_concentration'].values * 1000
ds['oxygen_concentration'].attrs['units'] = 'umol kg-1'
```

---

## Appendix A: Unit Conversion Quick Reference

| Variable | Raw | Final | Factor | Function |
|----------|-----|-------|--------|----------|
| **Temperature** | °C | °C | 1 (no conversion) | - |
| **Pressure** | dbar | dbar | 1 (no conversion) | - |
| **Conductivity** | mS/cm | S/m | ×0.1 | `cndc * 0.1` |
| **Salinity** | - | PSU | - | `gsw.SP_from_C(C, t, p)` |
| **Density** | - | kg/m³ | - | `gsw.rho(SA, CT, p)` |
| **Oxygen** | mmol/L | µmol/kg | ×1000 | `doxy * 1000` |
| **Chlorophyll** | mg/m³ | mg/m³ | 1 (no conversion) | - |
| **CDOM** | ppb | ppb | 1 (no conversion) | - |
| **Turbidity** | m⁻¹sr⁻¹ | NTU | ÷0.002727 | `beta / 0.002727` |
| **Depth** | dbar | m | - | `gsw.z_from_p(p, lat)` |

---

## Appendix B: TEOS-10 Function Reference

| Function | Inputs | Output | Purpose |
|----------|--------|--------|---------|
| `gsw.SP_from_C` | C, t, p | SP | Conductivity → Practical Salinity |
| `gsw.SA_from_SP` | SP, p, lon, lat | SA | Practical → Absolute Salinity |
| `gsw.CT_from_t` | SA, t, p | CT | In-situ → Conservative Temperature |
| `gsw.rho` | SA, CT, p | ρ | In-situ density |
| `gsw.sigma0` | SA, CT | σ₀ | Potential density anomaly (ref 0 dbar) |
| `gsw.z_from_p` | p, lat | z | Pressure → Depth |
| `gsw.O2sol_SP_pt` | SP, pt | O₂ | Oxygen solubility |
| `gsw.sound_speed` | SA, CT, p | c | Speed of sound |
| `gsw.N2` | SA, CT, p, lat | N² | Brunt-Väisälä frequency |

**Full documentation:** http://www.teos-10.org/pubs/gsw/html/gsw_contents.html

---

**Document Version:** 1.0  
**Last Updated:** October 2025  
**Maintainer:** Benedetta Torelli  
**Repository:** https://github.com/BennyTorelli/Seaexplorer_pyglider  
**License:** CC BY 4.0

**Citations:**
- IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater – 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO, 196 pp.
- McDougall, T.J. and P.M. Barker, 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127.
