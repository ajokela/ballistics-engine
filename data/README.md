# Standard Drag Coefficient Tables

## Overview

This directory contains the industry-standard G1 and G7 drag coefficient tables used for external ballistics calculations. These tables provide drag coefficient (Cd) values as a function of Mach number for two standard reference projectiles.

## Data Files

- `g1.csv` - G1 drag function (flat base, 2 caliber ogive)
- `g7.csv` - G7 drag function (boat tail, 7 caliber ogive)
- `g1.npy` - Binary NumPy format for faster loading
- `g7.npy` - Binary NumPy format for faster loading

## Data Format

CSV files contain two columns:
1. **Mach number** (0.0 to 5.0)
2. **Drag coefficient** (Cd)

Example:
```
0.00,0.2629
0.05,0.2558
0.10,0.2487
...
```

## Provenance

This data originates from:
- **Source**: US Army Ballistic Research Laboratory (BRL)
- **Location**: Aberdeen Proving Ground, Maryland
- **Publication**: Winchester-Western Division, 1965
- **Refinement**: Robert L. McCoy and BRL team
- **Status**: PUBLIC DOMAIN (US Government work)

## Citation

When referencing this data in academic or technical work:

**Primary Citation:**
McCoy, Robert L. "McDrag - A Computer Program for Estimating the Drag Coefficients of Projectiles." Technical Report ARBRL-TR-02293, Army Ballistic Research Lab, Aberdeen Proving Ground, MD, February 1981.

**Original Data Source:**
E.D. Lowry, Winchester-Western Division, Olin Mathieson Chemical Corporation, "Exterior Ballistic Tables Based on Numerical Integration," May 4, 1965.

## Usage

```python
import numpy as np

# Load drag table
g1_data = np.loadtxt('g1.csv', delimiter=',')
mach_values = g1_data[:, 0]
cd_values = g1_data[:, 1]

# Interpolate for specific Mach number
from scipy.interpolate import interp1d
drag_func = interp1d(mach_values, cd_values, kind='cubic')
cd_at_mach_0_9 = drag_func(0.9)
```

## Technical Notes

### G1 Standard Projectile
- Flat base design
- 2 caliber radius ogival nose
- Total length: 3.28 calibers
- Best for: Traditional flat-base bullets
- Popular with: Older ammunition, round nose bullets

### G7 Standard Projectile
- Boat tail base
- 7 caliber radius ogival nose
- Secant ogive design
- Best for: Modern long-range bullets
- Popular with: Match bullets, VLD designs

### Atmosphere Standard
Data calculated using ICAO Standard Atmosphere:
- Sea level temperature: 15°C (59°F)
- Sea level pressure: 1013.25 hPa
- Temperature lapse rate: -6.5°C/km

## Historical Significance

These drag functions have been the foundation of military and civilian ballistics calculations since 1965. They represent thousands of hours of experimental work at Aberdeen Proving Ground's ballistic ranges, using both physical testing and early computer modeling.

## Accuracy Considerations

- **Subsonic (Mach < 0.8)**: Excellent accuracy
- **Transonic (Mach 0.8-1.2)**: Good accuracy, some interpolation needed
- **Supersonic (Mach > 1.2)**: Very good accuracy
- **High supersonic (Mach > 3.0)**: Limited data, extrapolation required

For maximum accuracy with specific bullets, consider using:
- Custom drag models (CDMs) from Doppler radar data
- Manufacturer-specific BC values
- Multiple BC values for different velocity ranges

## License

PUBLIC DOMAIN - No restrictions on use. See LICENSE file for details.