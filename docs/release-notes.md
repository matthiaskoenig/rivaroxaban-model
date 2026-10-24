# Release notes

## Release notes for rivaroxaban-model 0.8.0

Made from pkdb_models 3a0b49144b with pkpdsim 0.11.0.

The changes since 0.7.0.

### Core libraries

The versions of the core libraries cannot be compared, there is no lock file in the previous release. The release uses:

- pkpdsim 0.11.0
- sbmlsim 0.13.0
- sbmlutils 0.17.0

### Libraries

The versions cannot be compared, there is no lock file in the previous release.

### Data

- Added: `apixaban/Frost2014`, `apixaban/Kreutz2017`, `edoxaban/Lenard2024`, `edoxaban/Lenard2025`, `edoxaban/Metze2021`, `edoxaban/Rohr2024`, `rivaroxaban/Alalawneh2024`, `rivaroxaban/Ding2020`, `rivaroxaban/GenisNajera2022`, `rivaroxaban/Graff2007`, `rivaroxaban/Greenblatt2018`, `rivaroxaban/Jiang2010`, `rivaroxaban/Kubitza2005`, `rivaroxaban/Kubitza2005a`, `rivaroxaban/Kubitza2006`, `rivaroxaban/Kubitza2013`, `rivaroxaban/Liu2024`, `rivaroxaban/Moore2014`, `rivaroxaban/Santos2024`, `rivaroxaban/Tang2025`, `rivaroxaban/Weinz2009`, `rivaroxaban/Zhao2009`

### Analysis

#### Models

- `rivaroxaban_body.xml`: added 2 species, 2 reactions, 5 parameters; removed 7 parameters; changed values `EC50_riv_Xa` 0.000295895 → 0.000558136, `EC50_riv_aPTT` 0.000408918 → 0.000845487, `Emax_Xa` 0.686093 → 1, `Emax_aPTT` 0.950019 → 1.59525
- `rivaroxaban_body_flat.xml`: added 2 species, 4 reactions, 35 parameters; removed 1 species, 19 parameters; changed values `EC50_riv_Xa` 0.000295895 → 0.000558136, `EC50_riv_aPTT` 0.000408918 → 0.000845487, `Emax_Xa` 0.686093 → 1, `Emax_aPTT` 0.950019 → 1.59525, `GU__F_riv_abs` 0.82 → 0.93, `GU__RIVABS_k` 0.0520279 → 0.0818778, `GU__RXEXC_k` 0.1 → 0.000462963, `LI__RXEX_BI_k` 0.0001 → 0.179531
- `rivaroxaban_coagulation.xml`: added 1 parameter; removed 2 parameters; changed values `EC50_riv_Xa` 0.000295895 → 0.000558136, `EC50_riv_aPTT` 0.000408918 → 0.000845487, `Emax_Xa` 0.686093 → 1, `Emax_aPTT` 0.950019 → 1.59525
- `rivaroxaban_intestine.xml`: added 2 reactions, 20 parameters; removed 1 species, 2 parameters; changed values `F_riv_abs` 0.82 → 0.93, `RIVABS_k` 0.0520279 → 0.0818778, `RXEXC_k` 0.1 → 0.000462963
- `rivaroxaban_kidney.xml`: added 6 parameters; removed 6 parameters
- `rivaroxaban_liver.xml`: added 5 parameters; removed 3 parameters; changed values `RXEX_BI_k` 0.0001 → 0.179531

#### Experiments

- Added: `Alalawneh2024`, `CoagulationExperiment`, `Ding2020`, `DoseDependencyExperiment`, `Frost2014`, `GenisNajera2022`, `Graff2007`, `Greenblatt2018`, `Jiang2010`, `Kreutz2017`, `Kubitza2005`, `Kubitza2005a`, `Kubitza2006`, `Kubitza2013`, `Lenard2024`, `Lenard2025`, `Liu2024`, `Metze2021`, `Moore2014`, `RivaroxabanParameterScan`, `Rohr2024`, `Santos2024`, `Tang2025`, `Weinz2009`, `Zhao2009`

#### Fits

- `PD`: `PD_20261018_082509__94e87`, no fit record in the previous release
- `PK`: `PK_20261018_075718__82026`, no fit record in the previous release

## Release notes for rivaroxaban-model v0.7.1

We are pleased to release a new version of the rivaroxaban-model including the 
following changes:

- model preprint added
- FAIR assessment added
- uv update
- updated github actions

Your rivaroxaban-model team

## Release notes for rivaroxaban-model v0.7.0

We are pleased to release a new version of the rivaroxaban-model including the 
following changes:

- updated pharmacokinetics
- updated pharmacodynamics
- simulations
- parameter optimization

Your rivaroxaban-model team

## Release notes for rivaroxaban-model v0.1.1

We are pleased to release a new version of the rivaroxaban-model including the 
following changes:

- initial model release
- zenodo integration

Your rivaroxaban-model team

## Release notes for rivaroxaban-model v0.1.0

We are pleased to release a new version of the rivaroxaban-model including the 
following changes:

- initial model release

Your rivaroxaban-model team

