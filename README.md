[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17091694.svg)](https://doi.org/10.5281/zenodo.17091694)
[![Release](https://img.shields.io/github/v/release/matthiaskoenig/rapamycin-model)](https://github.com/matthiaskoenig/rapamycin-model/releases)
[![Documentation](https://github.com/matthiaskoenig/rapamycin-model/actions/workflows/docs.yml/badge.svg)](https://matthiaskoenig.github.io/rapamycin-model/)
[![CI](https://github.com/matthiaskoenig/rapamycin-model/actions/workflows/ci.yml/badge.svg)](https://github.com/matthiaskoenig/rapamycin-model/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/docker/v/matthiaskoenig/rapamycin?sort=semver&label=docker)](https://hub.docker.com/r/matthiaskoenig/rapamycin)
[![Python](https://img.shields.io/badge/python-3.14-blue)](https://www.python.org)
[![Code license: MIT](https://img.shields.io/badge/code%20license-MIT-blue)](https://opensource.org/license/MIT)
[![Model license: CC BY-SA 4.0](https://img.shields.io/badge/model%20license-CC%20BY--SA%204.0-lightgrey)](https://creativecommons.org/licenses/by-sa/4.0/)

# Physiologically based pharmacokinetic (PBPK) model of rapamycin

This repository provides the rapamycin (sirolimus) physiologically based pharmacokinetic (PBPK) model: its oral absorption, its metabolism by CYP3A4 and CYP3A5 in the intestine and the liver, its distribution in whole blood and the tissues and the excretion of its metabolites.

Documentation: <https://matthiaskoenig.github.io/rapamycin-model/>

## Quick start

```bash
git clone https://github.com/matthiaskoenig/rapamycin-model.git && cd rapamycin-model
uv sync
uv run run_simulation pkdb_models.models.rapamycin --action all -r results
```

## How to cite

> Jesionek, M. & König, M. (2026). *Physiologically based pharmacokinetic (PBPK) model of rapamycin*. Zenodo. https://doi.org/10.5281/zenodo.17091694

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

&copy; 2025-2026 Monika Jesionek and Matthias König
