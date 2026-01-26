# Online Service Terms and Privacy

This document describes the terms of use and privacy practices for the `--online` feature of the ballistics-engine CLI.

## Overview

The ballistics-engine includes an optional `--online` flag that sends trajectory calculations to a cloud-based API for ML-enhanced predictions. **This feature is entirely optional** - all core functionality works locally without any network communication.

## Open Source vs. Proprietary

| Component | License | Network Required |
|-----------|---------|------------------|
| ballistics-engine library | MIT (Open Source) | No |
| CLI tool | MIT (Open Source) | No |
| Local physics solver | MIT (Open Source) | No |
| FFI bindings | MIT (Open Source) | No |
| **Online API service** | **Proprietary** | **Yes** |
| **ML models** | **Proprietary** | **Yes** |

The online service (API, ML models, cloud infrastructure) is **proprietary software** and is not covered by the MIT license. The source code for the online service will not be made available.

## Data Collection (--online flag only)

When you use the `--online` flag, the following data is transmitted to `api.ballistics.7.62x51mm.sh`:

### Trajectory Parameters
- Ballistic coefficient (value and type: G1, G7, etc.)
- Bullet specifications (mass, diameter)
- Muzzle velocity and target distance
- Zero distance and shooting angle
- Wind conditions (speed, direction)
- Atmospheric conditions (temperature, pressure, humidity, altitude)
- Geographic coordinates (latitude, longitude) if provided
- Weather model settings (wind shear model, weather zones)

### Technical Data
- Your IP address
- Client version (e.g., `ballistics-cli/0.13.29`)
- Request timestamp
- Response status and latency

## Data Usage

Your data is used **solely for service delivery**:
- Processing your trajectory calculation request
- Returning results to your client
- Maintaining operational logs for debugging

**We do NOT:**
- Use your data to train or improve ML models
- Share your data with third parties
- Use your data for marketing or advertising
- Build user profiles from your requests

## Data Retention

- Request logs are retained for **30 days**
- Logs are automatically deleted after 30 days by Google Cloud Platform
- We do not maintain any permanent database of trajectory calculations

## Data Location

The online service is hosted on Google Cloud Platform (Cloud Run) in the United States. All data is transmitted over HTTPS.

## Opting Out

To use ballistics-engine without any network communication:

### Option 1: Don't use --online
Simply omit the `--online` flag. All local calculations work without network access.

```bash
# Local calculation (no network)
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000

# Online calculation (sends data to API)
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000 --online
```

### Option 2: Build without online feature
```bash
cargo install ballistics-engine --no-default-features
```

This removes the `--online` flag entirely from the CLI.

## No Warranty

The online service is provided "AS IS" without warranty of any kind. Ballistics calculations are for informational purposes only. You are solely responsible for verifying calculations and for the safe and legal use of firearms.

## Service Availability

We reserve the right to modify, suspend, or discontinue the online service at any time. The open source CLI will continue to function for local calculations regardless of online service status.

## Contact

For questions about the online service:
- Email: email@tinycomputers.io
- Website: https://tinycomputers.io
- Full terms: https://ballistics.rs/terms.html
- Full privacy policy: https://ballistics.rs/privacy.html
