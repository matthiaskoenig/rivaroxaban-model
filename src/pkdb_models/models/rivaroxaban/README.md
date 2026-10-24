# Rivaroxaban

Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of the factor Xa inhibitor rivaroxaban: its oral dosing as tablets, whose dissolution in the intestine saturates with the dose and is raised by food (which also slows the gastric emptying), and as oral solutions, and its intravenous dosing, its distribution into the tissues, its hepatic metabolism (CYP3A4, CYP2J2 and hydrolysis, the metabolites lumped in rivaroxaban equivalents) with the CYP3A4 activity of the liver and the P-gp/BCRP activity of the kidney as factors of the interactions, the renal excretion of rivaroxaban by filtration and active secretion relative to the glomerular filtration rate, the excretion of the metabolites in urine and, via the bile, in feces, and the anticoagulant effects of the plasma concentration (prothrombin time linear, activated partial thromboplastin time and inhibition of factor Xa as Emax models). The model describes the effects of the dose, the form, food, renal impairment, cirrhosis and inhibitors and inducers of CYP3A4 and P-gp.

The published version of the model is [rivaroxaban-model](https://github.com/matthiaskoenig/rivaroxaban-model) ([doi:10.5281/zenodo.15157735](https://doi.org/10.5281/zenodo.15157735)).

## Structure

- `__init__.py`: the paths of the model and its model object `MODEL` with its release metadata and its adopted fits, which the commands below run (pkpdsim).
- `models/`: the tissue models (`model_intestine.py`, `model_liver.py`, `model_kidney.py`), the anticoagulant effects (`model_coagulation.py`) and the whole-body model (`model_body.py`). `models/factory.py` creates the SBML models, the flattened whole-body model which is simulated, the ODEs and the COMBINE archive in `models/output/`; the fitted values come from the adopted fits (`models/fits.py`, `models/output/fits/`).
- `experiments/`: the simulation experiments of the clinical studies (`studies/`), the dose dependency and the coagulation (`misc/`) and the parameter scans of the dose, the renal function, the cirrhosis and the fed state (`scans/`). The fit mappings are described with the vocabulary of `pkpdsim.metadata`.
- `fitting/`: the fit problems of the pharmacokinetics (`PK`) and of the anticoagulant effects on the adopted pharmacokinetics (`PD`), their parameters and the selection of the training, validation (interactions, doses, populations, limitations), aggregate and outlier data.
- `sensitivity/`: the sensitivity analysis of the pharmacokinetics of rivaroxaban and its metabolites and of the maxima of the effects.
- `rivaroxaban_pk.py`: the pharmacokinetic and pharmacodynamic parameters of the scans, the non-compartmental analysis of `pkpdutils`.
- `tolerances.py`: the absolute tolerances of the integrator per state.
- `simulations.py`: the simulation experiments and their groups, `MODEL` runs them.
- `data/`: the data of the studies, a subset of the PK-DB studies (`data/collect_data.py`).
- `TODO.md`: the open issues of the model.

## Commands

```bash
# models, simulations
run_simulation pkdb_models.models.rivaroxaban --action factory
run_simulation pkdb_models.models.rivaroxaban --action list_experiments
run_simulation pkdb_models.models.rivaroxaban --action simulate --experiments all
run_simulation pkdb_models.models.rivaroxaban --action simulate --experiments all --uncertainty

# fits, adoption of a fit, reports of stored parameters, identifiability
run_fit pkdb_models.models.rivaroxaban --subset=PK --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.rivaroxaban --subset=PD --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.rivaroxaban --adopt PK_<fit id>
run_report pkdb_models.models.rivaroxaban --subset=PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.rivaroxaban --subset=PK results/fit/PK_<fit id>/parameters.json

# sensitivity analysis (local and sampling by default)
python -m pkdb_models.models.rivaroxaban.sensitivity.sensitivity_analysis local
```

The results are written into `results/` of the model.

## Adopted fits and uncertainty

The values of the fitted parameters come from the adopted fits of the model,
not from numbers in the code: `PK_20261018_075718__82026` (the
pharmacokinetics) and `PD_20261018_082509__94e87` (PT, aPTT and the inhibition
of factor Xa). `run_fit pkdb_models.models.rivaroxaban --adopt <fit id>`
copies the record of a finished fit (its parameters, problem, profile
likelihood and Fisher information) into `models/output/fits/<fit id>/`, which is
tracked; the ids of the adopted fits are listed in `fits` of the model object
(`__init__.py`). The factory takes the fitted values and their uncertainty from
the records (`models/fits.py`) and writes the uncertainty into the SBML
(`distrib:uncertainty`). `run_simulation ... --uncertainty` runs the
experiments over draws of these distributions and draws the simulations as the
median with a band of quantiles into `results/simulation_uncertainty`. A refit
is: fit, adopt, replace the id in `fits`, rebuild the model. See the
documentation of the repository on fitting (`docs/fitting.md`).
