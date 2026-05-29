# SedimentChronology_Workflow

A guided R workflow for dating lake sediment cores using short- and long-term radioisotope chronologies.

## Repository Structure

```
SedimentChronology_Workflow/
├── RERCA/   # Recent rates of carbon accumulation (Pb-210 dating)
│   ├── data/
│   │   ├── template_pb210_data.csv   # Fill this in with your core data
│   │   └── example_pb210_data.csv    # Example dataset from a real core
│   ├── 01_pb210_CRS.R                # CRS model using the pb210 package
│   ├── 02_rplum.R                    # Bayesian Pb-210 with rplum
│   ├── 03_serac.R                    # CRS/CIC/CFCS models using serac
│   └── 04_bayesian_pb210.R           # Manual Bayesian approach (Dunnington 2019)
└── LORCA/   # Long-term rates of carbon accumulation (C-14 dating) — coming soon
```

## RERCA — Pb-210 Short-term Chronology

Four complementary approaches are provided. Run them in order, or jump to the
method you need. Each script is self-contained and walks you through data
loading, background estimation, model fitting, and output.

| Script | Package | Method |
|--------|---------|--------|
| `01_pb210_CRS.R` | [pb210](https://github.com/paleolimbot/pb210) | CRS (Constant Rate of Supply) |
| `02_rplum.R` | [rplum](https://CRAN.R-project.org/package=rplum) | Bayesian Pb-210 |
| `03_serac.R` | [serac](https://github.com/rosalieb/serac) | CRS / CIC / CFCS |
| `04_bayesian_pb210.R` | base R + Stan | Manual Bayesian (Dunnington 2019) |

## Quick Start

1. Copy `RERCA/data/template_pb210_data.csv` and fill in your measurements.
2. Open the script for the method you want to use.
3. Set `data_file` at the top of the script to point to your filled-in CSV.
4. Run the script section-by-section following the step comments.

## Data Requirements

At minimum you need:
- Upper and lower depth (cm) for each sediment slice
- Dry bulk density (g/cm³)
- Total Pb-210 activity (DPM/g dry weight) with counting error (±1 SD)

Optional but recommended:
- Supported Pb-210 / Ra-226 activity — improves background estimation
- Cs-137 activity — independent time marker for model validation
- Loss on ignition (%) — needed for carbon accumulation rates
- Water content / % moisture

## LORCA — C-14 Long-term Chronology

Placeholder — will be added in a future update.
