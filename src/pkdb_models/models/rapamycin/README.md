# Rapamycin

Physiologically based pharmacokinetic (PBPK) model of the mTOR inhibitor rapamycin (sirolimus): its oral absorption (dissolution, fraction absorbed, metabolism by CYP3A4/5 in the gut wall and efflux by P-glycoprotein), its distribution in whole blood (bound in the erythrocytes, the plasma is the blood over the blood to plasma ratio) and in the tissues (with a slow binding which gives the terminal half-life of about 60 hr), its hepatic metabolism by CYP3A4/5 and the excretion of its metabolites with the bile, the feces and the urine. The model describes the effects of the CYP3A4 and CYP3A5 genotypes (CYP3A5 *3/*3 is the reference), of food, of hepatic and renal impairment and of coadministered inhibitors and inducers of CYP3A and P-glycoprotein such as cyclosporine, ritonavir and rifampin.

The published version of the model is [rapamycin-model](https://github.com/matthiaskoenig/rapamycin-model) ([doi:10.5281/zenodo.17091694](https://doi.org/10.5281/zenodo.17091694)).

## Structure

- `__init__.py`: the paths of the model and its model object `MODEL` with its release metadata and its adopted fits, which the commands below run (pkpdsim).
- `models/`: the tissue models (`model_intestine.py`, `model_liver.py`, `model_kidney.py`) and the whole-body model (`model_body.py`). `models/factory.py` creates the SBML models, the flattened whole-body model which is simulated, the ODEs and the COMBINE archive in `models/output/`; the fitted values come from the adopted fits (`models/fits.py`, `models/output/fits/`).
- `experiments/`: the simulation experiments of the clinical studies (`studies/`), the dose dependency (`misc/`) and the parameter scans of dose, food, hepatic and renal function, CYP3A4 and CYP3A5 activity and P-glycoprotein (`scans/`). The fit mappings are described with the vocabulary of `pkpdsim.metadata`.
- `fitting/`: the fit problem of the pharmacokinetics (`PK`), its parameters and the selection of the training, validation and outlier data.
- `sensitivity/`: the sensitivity analysis of the pharmacokinetics.
- `rapamycin_pk.py`: the pharmacokinetic parameters of the scans, the non-compartmental analysis of `pkpdutils`.
- `tolerances.py`: the absolute tolerances of the integrator per state.
- `simulations.py`: the simulation experiments and their groups, `MODEL` runs them.
- `data/`: the data of the studies, a subset of the PK-DB studies (`data/collect_data.py`).
- `TODO.md`: the open issues of the model.

## Commands

```bash
# models, simulations
run_simulation pkdb_models.models.rapamycin --action factory
run_simulation pkdb_models.models.rapamycin --action list_experiments
run_simulation pkdb_models.models.rapamycin --action simulate --experiments all
run_simulation pkdb_models.models.rapamycin --action simulate --experiments all --uncertainty

# fits, adoption of a fit, reports of stored parameters, identifiability
run_fit pkdb_models.models.rapamycin --subset=PK --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.rapamycin --adopt PK_<fit id>
run_report pkdb_models.models.rapamycin --subset=PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.rapamycin --subset=PK results/fit/PK_<fit id>/parameters.json

# sensitivity analysis (local by default)
python -m pkdb_models.models.rapamycin.sensitivity.sensitivity_analysis
```

The results are written into `results/` of the model.

## Adopted fits and uncertainty

The values of the fitted parameters come from the adopted fits of the model,
not from numbers in the code. `run_fit pkdb_models.models.rapamycin --adopt <fit id>`
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
