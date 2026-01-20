# Geodetic Inversion Package for MATLAB

This package performs **kinematic finite-fault slip inversion** using multi-source geodetic observations:

- **InSAR phase & phase gradients**
- **POT/MAI** (azimuth offsets)
- **GNSS** displacements
- **Optical** displacement measurements

The framework supports flexible data combinations and complex fault geometries while maintaining computational efficiency.

---

## Directory Structure
```text
Geodetic_Inversion/
├── example/
│   └── Ridgecrest/
│       ├── input/          # GNSS, InSAR, phase-gradient observations
│       └── model/          # Inversion configuration files and slip model
└── src/
    ├── Subsample.m            # Data subsampling
    ├── main_inv.m             # Main inversion entry
    └── main_inv_ridgecrest.m
```
## 📊 Data Preparation
*Perform quality control and preprocessing before inversion. All coordinates use WGS84 datum.*

### 📡 InSAR Phase/Phase Gradient Data
**File format (space/tab-delimited):**
| Column | Description                          | Units |
|--------|--------------------------------------|-------|
| 1      | Longitude                            | deg   |
| 2      | Latitude                             | deg   |
| 3      | Topography                           | m     |
| 4      | Look vector (East component)         | -     |
| 5      | Look vector (North component)        | -     |
| 6      | Look vector (Vertical component)     | -     |
| 7      | LOS displacement or phase gradient   | mm/~  |
| 8      | Sampling ratio (relative data weight)| -     |

> **Note:** if not specified Sampling ratio (relative data weight), set 1 is ok

### 🧭 POT/MAI Data
**File format (space/tab-delimited):**
| Column | Description               | Units |
|--------|---------------------------|-------|
| 1      | Longitude                 | deg   |
| 2      | Latitude                  | deg   |
| 3      | Topography                | m     |
| 4      | Sampling ratio or std. error | mm    |

### 📡 GNSS Data
**File format (space/tab-delimited):**
| Column | Description          | Units |
|--------|----------------------|-------|
| 1      | Longitude            | deg   |
| 2      | Latitude             | deg   |
| 3      | East displacement    | mm    |
| 4      | North displacement   | mm    |
| 5      | Vertical displacement| mm    |
| 6      | East error           | mm    |
| 7      | North error          | mm    |
| 8      | Vertical error       | mm    |
| 9      | Topography           | m     |

---

## ⚙️ Configuration (`config.inv`)
### Dataset Definition
```text
num_des_sources      = 2   # Descending InSAR datasets
num_asc_sources      = 1   # Ascending InSAR datasets
num_x_grd_sources    = 0   # West-East phase gradients
num_y_grd_sources    = 1   # South-North phase gradients
num_azi_sources      = 1   # POT/MAI datasets
num_gps_sources      = 3   # GNSS stations
phi1                 = 349.25 # Satellite track angle (degrees)
```

### GNSS Component Selection
```text
# Enable (1) or disable (0) components per dataset:
gps1h = 1   # Horizontal (E,N) for GNSS #1
gps1v = 0   # Vertical component disabled
gps2h = 1
gps2v = 1
```

### Data File Assignment
Assign paths in the `{data_files}` block of `config.inv`:
```text
des1 = t71_des.lltnde
asc1 = t64_asc_edit.lltnde
asc2 = t65_asc_alos2_new.lltnde
asc3 = t66_asc_alos2.lltnde

dec_grdx = Rc_ph_grdx_dec_masked.lltnde
asc_grdx = Rc_ph_grdx_asc_masked.lltnde
dec_grdy = Rc_ph_grdy_dec_masked.lltnde
asc_grdy = Rc_ph_grdy_asc_masked.lltnde

azi1 = t65_mai.llde

gps1 = GPS71_new.lldet
gps2 = GPS71_new.lldet

gov1 = optical_offset.lldeat
```
---

## 🌋 Fault Geometry (`model` block)
{model_params}
   num_of_sources = 8
   poisson_ratio = 0.25
{origin}
   xo = -117.5990  
   yo = 35.7700
{trace1}
       x = -8.379680e+03
       y = 1.216606e+04
       z = 0.000000e+00
     len = 1.686743e+04
     wid = 2.000000e+04
     dip = 90
  strike = 1.357613e+02

### Basic Parameters
```text
num_of_source = 3    # Fault segments
xo = 142.365         # Epicenter longitude (deg)
yo = 38.322          # Epicenter latitude (deg)
```

### Geometry Specification (Choose **ONE** method)
#### Option 1: Manual UTM Definition
Directly define patches in local UTM coordinates (advanced users).

#### Option 2: `trace2patch` Utility (Recommended)
```matlab
% Define fault trace endpoints (geographic coordinates)
lon_endpoints = [lon1, lon2, lon2, lon3, ..., lon_n];
lat_endpoints = [lat1, lat2, lat2, lat3, ..., lat_n];
lon_epi = lon0;  % Epicenter longitude
lat_epi = lat0;  % Epicenter latitude

% Generate patches (default: 20km width, 90° dip)
[patch_params] = trace2patch(lon_endpoints, lat_endpoints, lon_epi, lat_epi);
```
> **Defaults:** 20 km fault width, 90° dip angle. Modify in output files as needed.

---

## 🔧 Inversion Parameters (`inversion` block)
### Discretization & Slip Constraints
```text
top_patch_width        = 500    # Min patch width (m)
top_patch_length       = 500    # Min patch length (m)
patch_increment_factor = 1.5    # Resolution increase with depth

strike_slip = 1    # Enable strike-slip
dip_slip    = 1    # Enable dip-slip
normal_slip = 0    # Disable normal-slip (typical for most earthquakes)
```

### Regularization
```text
positivity_max          = 1000    # Max slip magnitude (m)
bottom_zero_constraint  = 1     # Enforce zero slip at bottom (1=yes)
smooth_factor           = 0.1   # Smoothing strength (1e-2 to 1e1)
smooth_between_segments = 1     # Cross-segment smoothing
smooth_dip_over_strike  = 0.5   # Dip vs. strike smoothing ratio
```

### Data Weighting & Corrections
```text
weight_phase  = 1.0   # InSAR phase weight
weight_ph_grd = 8000   # Phase gradient weight
weight_azi    = 0.1   # POT/MAI weight
weight_gps    = 0.1   # GNSS weight

remove_ramp          = 1  # Remove orbital/ionospheric ramps
consider_topography  = 0  # Ignore topography effects (0=no)
switch_phase         = 1  # Phase sign convention (+1/-1)
```

---

## 🔗 Advanced Constraints
### Inter-Segment Smoothness (`smooth` block)
```text
num_seg_smooth = 2       # Adjacent segment pairs
smo1 = [1, 2]            # Smooth between seg1 & seg2
smo2 = [2, 3]            # Smooth between seg2 & seg3

num_inter_smooth = 1     # Complex fault connections
smoi1 = [1, 3, 2]        # Connect seg1-seg3 with type-2 smoothing (artificial connection)
```

### Fault Edge Constraints (`edge_constraints` block)
```text
bot = 1       # Zero slip at bottom edge
top = 0       # No constraint at top edge
side = 1      # Enable side constraints
num_side = 2  # Segments with side constraints
side1 = [1, 1, 0]  # Segment1: left-edge=constrained, right-edge=free
side2 = [3, 0, 1]  # Segment3: left-edge=free, right-edge=constrained
```

---

## ▶️ Running the Inversion
1. Configure `config.inv` with your parameters
2. Set file paths in `main_inv.m`:
   ```matlab
    earthquake = 'Ridgecrest';  %
    input_path = ['/home/path/to/Geodetic_Inversion/example/',earthquake,'/input/'];
    model_path = ['/home/path/to/Geodetic_Inversion/example/',earthquake,'/model/'];
    addpath '/home/path/to/Geodetic_Inversion/src'
   ```
3. Execute the inversion:
   ```matlab
   main_inv.m; % Start computation
   ```
---

## 📚 References
1. **Core Methodology**  
   Xu, X., Tong, X., Sandwell, D. T., Milliner, C. W., Dolan, J. F., Hollingsworth, J., ... & Ayoub, F. (2016). [Refining the shallow slip deficit](https://doi.org/10.1093/gji/ggv532). *Geophysical Journal International*, 204(3), 1867-1886.

2. **Phase-Gradient Constraints**  
   Zhang, Y., & Xu, X. (2024). [Constraining shallow slip deficit with phase gradient data](https://doi.org/10.1093/gji/ggaf427). *Geophysical Journal International*, 236(1), 123-145.  
