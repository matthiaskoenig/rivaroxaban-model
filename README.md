[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15157735.svg)](https://doi.org/10.5281/zenodo.15157735)
[![Release](https://img.shields.io/github/v/release/matthiaskoenig/rivaroxaban-model)](https://github.com/matthiaskoenig/rivaroxaban-model/releases)
[![Documentation](https://github.com/matthiaskoenig/rivaroxaban-model/actions/workflows/docs.yml/badge.svg)](https://matthiaskoenig.github.io/rivaroxaban-model/)
[![CI](https://github.com/matthiaskoenig/rivaroxaban-model/actions/workflows/ci.yml/badge.svg)](https://github.com/matthiaskoenig/rivaroxaban-model/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/docker/v/matthiaskoenig/rivaroxaban?sort=semver&label=docker)](https://hub.docker.com/r/matthiaskoenig/rivaroxaban)
[![Python](https://img.shields.io/badge/python-3.14-blue)](https://www.python.org)
[![Code license: MIT](https://img.shields.io/badge/code%20license-MIT-blue)](https://opensource.org/license/MIT)
[![Model license: CC BY-SA 4.0](https://img.shields.io/badge/model%20license-CC%20BY--SA%204.0-lightgrey)](https://creativecommons.org/licenses/by-sa/4.0/)

# Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of rivaroxaban

This repository provides the rivaroxaban physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model: the dose and food dependent dissolution and the absorption, the tissue distribution, the hepatic metabolism (CYP3A4, CYP2J2, hydrolysis), the renal excretion by filtration and active secretion, the excretion of the metabolites and the anticoagulant effects (PT, aPTT, inhibition of factor Xa).

Documentation: <https://matthiaskoenig.github.io/rivaroxaban-model/>

## Quick start

```bash
git clone https://github.com/matthiaskoenig/rivaroxaban-model.git && cd rivaroxaban-model
uv sync
uv run run_simulation pkdb_models.models.rivaroxaban --action all -r results
```

## How to cite

> Casabianca, E., Myshkina, M. & König, M. (2026). *Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of rivaroxaban*. Zenodo. https://doi.org/10.5281/zenodo.15157735

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

&copy; 2025-2026 Elisabetta Casabianca, Mariia Myshkina and Matthias König
