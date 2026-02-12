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

class Zimmerman2003(RapamycinSimulationExperiment):
    """Simulation experiment of Zimmerman2003.

    Cyclosporine blocks the P-GP transporter - that slows down the elimination of rapamycin.
    """

    interventions = ["RAP10", "RAP10_CSA300", "RAP10_CSA300B"]
    labels = {
        "RAP10": "rapamycin_RAP10",
        "RAP10_CSA300": "rapamycin_RAP10_CSA300",
        "RAP10_CSA300B": "rapamycin_RAP10_CSA300B",
    }
    colors = {
        "RAP10": "black",
        "RAP10_CSA300": "#e3256b",
        "RAP10_CSA300B": "#d65282",
    }
    f_cyp3a4 = {
        "RAP10": 1.0,  # no change
        "RAP10_CSA300": 0.6,  # less metabolism
        "RAP10_CSA300B": 0.8  # in between
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
                    steps=2000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(74.6, "kg"),
                        "PODOSE_rap": Q_(10, "mg"),

                        # CYP3A4 activity (inhibition by cyclosporine)
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
                    coadministration=Coadministration.NONE if intervention == "RAP10" else Coadministration.CYCLOSPORINE
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
    run_experiments(Zimmerman2003, output_dir=Zimmerman2003.__name__)
