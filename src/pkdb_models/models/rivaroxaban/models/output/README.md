[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15157735.svg)](https://doi.org/10.5281/zenodo.15157735)
[![GitHub Action](https://github.com/matthiaskoenig/rivaroxaban-model/actions/workflows/python.yml/badge.svg)](https://github.com/matthiaskoenig/rivaroxaban-model/actions/workflows/python.yml)
[![GitHub Action](https://github.com/matthiaskoenig/rivaroxaban-model/actions/workflows/docker.yml/badge.svg)](https://github.com/matthiaskoenig/rivaroxaban-model/actions/workflows/docker.yml)

# rivaroxaban model
This repository provides the rivaroxaban physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model: the oral dosing as tablets, whose dissolution saturates with the dose and is raised by food, and as solutions, the tissue distribution, the hepatic metabolism (CYP3A4, CYP2J2, hydrolysis) with CYP3A4 activity for interactions, the renal excretion by filtration and active secretion relative to the glomerular filtration rate (renal impairment), the excretion of the metabolites in urine and feces, cirrhosis and the anticoagulant effects of rivaroxaban (prothrombin time, activated partial thromboplastin time and the inhibition of factor Xa).

The model is distributed as [SBML](http://sbml.org) format available from [`rivaroxaban_body_flat.xml`](./models/rivaroxaban_body_flat.xml) with 
corresponding [SBML4humans model report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Frivaroxaban-model%2Fmain%2Fmodels%2Frivaroxaban_body_flat.xml) and [model equations](./models/rivaroxaban_body_flat.md).

The COMBINE archive is available from [`rivaroxaban_model.omex`](./rivaroxaban_model.omex).
The FAIR assessment is available from [`rivaroxaban_model_fair.xlsx`](./rivaroxaban_model_fair.xlsx).

![model overview](./figures/rivaroxaban_model.png)

### Comp submodels
* **liver** submodel [`rivaroxaban_liver.xml`](./models/rivaroxaban_liver.xml) with [SBML4humans report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Frivaroxaban-model%2Fmain%2Fmodels%2Frivaroxaban_liver.xml) and [equations](./models/rivaroxaban_liver.md).
* **kidney** submodel [`rivaroxaban_kidney.xml`](./models/rivaroxaban_kidney.xml) with [SBML4humans report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Frivaroxaban-model%2Fmain%2Fmodels%2Frivaroxaban_kidney.xml) and [equations](./models/rivaroxaban_kidney.md).
* **intestine** submodel [`rivaroxaban_intestine.xml`](./models/rivaroxaban_intestine.xml) with [SBML4humans report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Frivaroxaban-model%2Fmain%2Fmodels%2Frivaroxaban_intestine.xml) and [equations](./models/rivaroxaban_intestine.md).
* **coagulation** submodel [`rivaroxaban_coagulation.xml`](./models/rivaroxaban_coagulation.xml) with [SBML4humans report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Frivaroxaban-model%2Fmain%2Fmodels%2Frivaroxaban_coagulation.xml) and [equations](./models/rivaroxaban_coagulation.md).
* **whole-body** submodel [`rivaroxaban_body.xml`](./models/rivaroxaban_body.xml) with [SBML4humans report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Frivaroxaban-model%2Fmain%2Fmodels%2Frivaroxaban_body.xml) and [equations](./models/rivaroxaban_body.md).

## How to cite
To cite the model repository

> Casabianca, E. Myshkina, M. & König, M. (2026).
> *Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of rivaroxaban.*
> Zenodo. [https://doi.org/10.5281/zenodo.15157735](https://doi.org/10.5281/zenodo.15157735)

To cite the preprint

> Elisabetta Casabianca, Mariia Myshkina, Matthias König
> A Systems Pharmacology Approach to Rivaroxaban: Physiologically Based Modeling of Pharmacokinetics and Coagulation Dynamics
> Preprints 2025, 202507.1945. (preprint). doi:[10.20944/preprints202507.1945.v1](https://doi.org/10.20944/preprints202507.1945.v1)

## License

* Source Code: [MIT](https://opensource.org/licenses/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

## Run simulations
### python
Clone the repository 
```bash
git clone https://github.com/matthiaskoenig/rivaroxaban-model.git
cd rivaroxaban-model
```

#### uv
Run the complete analysis with uv (https://docs.astral.sh/uv/getting-started/installation/):
```bash
uv sync
uv run run_simulation pkdb_models.models.rivaroxaban -a all -r results
```

#### pip
If you use pip install the package via
```bash
pip install -e .
```
Run the complete analysis in the environment via:
```bash
run_simulation pkdb_models.models.rivaroxaban -a all -r results
```

### docker
Simulations can also be run within a docker container:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/rivaroxaban:latest /bin/bash
```

Run the complete analysis:
```bash
uv run run_simulation pkdb_models.models.rivaroxaban -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```


## Funding
Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B).
Mariia Myshkina was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within Priority Programme SPP 2311 by grant number 465194077 (Subproject SimLivA).

© 2025-2026 Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)
