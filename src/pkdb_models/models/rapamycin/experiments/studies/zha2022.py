from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)
from pkdb_models.models.rapamycin.experiments.metadata import (
    RapamycinMappingMetaData,
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, Coadministration
)

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rapamycin.helpers import run_experiments


class Zha2022(RapamycinSimulationExperiment):
    """Simulation experiment of Zha2022.

    The	ritonavir-containing 3D	regimen	resulted in a significant increase in sirolimus exposure, consistent with
    the	known strong inhibitory	effect of ritonavir on CYP3A requiring dose and/or frequency
    modification when co-administered with each other. - most important here RITONAVIR"""

    interventions = ["RAP2", "RAP05_3D"]
    colors = {
        "RAP2": "black",
        "RAP05_3D": "#e3256b",
    }
    doses = {
        "RAP2": 2,
        "RAP05_3D": 0.5,
    }
    f_cyp3a4 = {
        "RAP2": 1,
        "RAP05_3D": 0.08, #changed, effect was correct, but more inhibiting - lower decrease probably cause of the bounding
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                dset.unit_conversion("mean", 1 / self.Mr.rap)
                dsets[label] = dset

        #console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=500 * 60,  # [min] 48h
                    steps=3000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(77.2, "kg"),
                        "PODOSE_rap": Q_(self.doses[intervention], "mg"),

                        # CYP3A4 activity (inhibition by ritonavir)
                        "GU__f_cyp3a4": Q_(self.f_cyp3a4[intervention], "dimensionless"),
                        "LI__f_cyp3a4": Q_(self.f_cyp3a4[intervention], "dimensionless"),
                    },
                )]
            )

        #console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}

        for intervention in self.interventions:
            mappings[f"fm_rap_{intervention}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"rapamycin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_rap_{intervention}", xid="time", yid="[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.NR,
                    application_form=ApplicationForm.NR,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE if intervention == "RAP2" else Coadministration.RITONAVIR
                )
            )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time))
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)
        for intervention in self.interventions:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=intervention,
                color=self.colors[intervention],
            )
                # data
            plots[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=intervention,
                color=self.colors[intervention],
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    run_experiments(Zha2022, output_dir=Zha2022.__name__)
