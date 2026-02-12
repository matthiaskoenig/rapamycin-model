[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17091694.svg)](https://doi.org/10.5281/zenodo.17091694)
[![GitHub Action](https://github.com/matthiaskoenig/rapamycin-model/actions/workflows/python.yml/badge.svg)](https://github.com/matthiaskoenig/rapamycin-model/actions/workflows/python.yml)
[![GitHub Action](https://github.com/matthiaskoenig/rapamycin-model/actions/workflows/docker.yml/badge.svg)](https://github.com/matthiaskoenig/rapamycin-model/actions/workflows/docker.yml)

# rapamycin model
This repository provides the rapamycin physiologically based pharmacokinetics/ pharmacodynamics (PBPK/PD) model.

The model is distributed as [SBML](http://sbml.org) format available from [`rapamycin_body_flat.xml`](./models/rapamycin_body_flat.xml) with 
corresponding [SBML4humans model report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/rapamycin-model/main/models/rapamycin_body_flat.xml) and [model equations](./models/rapamycin_body_flat.md).

The COMBINE archive is available from [`rapamycin_model.omex`](./rapamycin_model.omex).

![model overview](./figures/rapamycin_model.png)

### Comp submodels
* **kidney** submodel [`rapamycin_liver.xml`](./models/rapamycin_liver.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/rapamycin-model/main/models/rapamycin_liver.xml) and [equations](./models/rapamycin_liver.md).
* **kidney** submodel [`rapamycin_kidney.xml`](./models/rapamycin_kidney.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/rapamycin-model/main/models/rapamycin_kidney.xml) and [equations](./models/rapamycin_kidney.md).
* **intestine** submodel [`rapamycin_intestine.xml`](./models/rapamycin_intestine.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/rapamycin-model/main/models/rapamycin_intestine.xml) and [equations](./models/rapamycin_intestine.md).
* **whole-body** submodel [`rapamycin_body.xml`](./models/rapamycin_body.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/rapamycin-model/main/models/rapamycin_body.xml) and [equations](./models/rapamycin_body.md).

## How to cite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17091694.svg)](https://doi.org/10.5281/zenodo.17091694)

> Jesionek, A. & König, M. (2026).
> *Physiologically based pharmacokinetic (PBPK) model of rapamycin.*   
> Zenodo. [https://doi.org/10.5281/zenodo.17091694](https://doi.org/10.5281/zenodo.17091694)

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.


## Run simulations
### python
Clone the repository 
```bash
git clone https://github.com/matthiaskoenig/rapamycin-model.git
cd rapamycin-model
```

#### uv
Run the complete analysis with uv (https://docs.astral.sh/uv/getting-started/installation/):
```bash
uv run run_rapamycin -a all -r results
```

#### pip
If you use pip install the package via
```bash
pip install -e .
```
Run the complete analysis in the environment via:
```bash
run run_rapamycin -a all -r results
```

### docker
Simulations can also be run within a docker container:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/rapamycin:latest /bin/bash
```

Run the complete analysis:
```bash
uv run run_rapamycin -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```

## Funding
Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B).

© 2025-2026 Monika Jesionek and Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)
