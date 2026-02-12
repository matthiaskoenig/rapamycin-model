from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData

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

class Tortorici2014(RapamycinSimulationExperiment):
    """Simulation experiment of Tortorici2014.

    Rifampin is a potent CYP3A4 inducer and speed up the metabolism of sirolimus.
    """

    interventions = ["RAP20", "RAP20_RIF600"]
    labels = {
        "RAP20": "rapamycin_RAP20",
        "RAP20_RIF600": "rapamycin_RAP20_RIF600",
    }
    colors = {
        "RAP20": "black",
        "RAP20_RIF600":"#e3256b",
    }
    f_cyp3a4 = {
        "RAP20": 1.0,  # no change in activity
        "RAP20_RIF600": 2.0  # speed up the metabolism (example value)
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                if not label.startswith("rapamycin_"):
                    continue
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
                    end=150 * 60,  # [min] 48h
                    steps=3000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(77.1, "kg"),
                        "PODOSE_rap": Q_(20, "mg"),

                        # CYP3A4 activity (inhibition by diltiazem)
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
                    self, task=f"task_rap_{intervention}", xid="time", yid=f"[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE if intervention == "RAP20" else Coadministration.RIFAMPIN
                )
            )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True
        )
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for intervention in self.interventions:
                # simulation
            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid=f"[Cveblood_rap]",
                label=intervention,
                color=self.colors[intervention],
            )
                # data
            plots[0].add_data(
                dataset=self.labels[intervention],
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=intervention,
                color=self.colors[intervention],
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    run_experiments(Tortorici2014, output_dir=Tortorici2014.__name__)
