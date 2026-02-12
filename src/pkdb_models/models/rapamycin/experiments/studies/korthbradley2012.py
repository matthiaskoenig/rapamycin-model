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


class KorthBradley2012(RapamycinSimulationExperiment):
    """Simulation experiment of KorthBradley2012."""

    interventions = ["RAP_T5", "RAP_NS5", "RAP_S5"]
    colors = {
        "RAP_T5": "black",
        "RAP_NS5": "#FFB98C",
        "RAP_S5": "#D28F66",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("rapamycin_"):
                   dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[label] = dset

        # console.print(dsets)
        #console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=150 * 60,  # [min]
                    steps=2000,
                    changes={

                        **self.default_changes(),

                        # physiological changes
                        "BW": Q_(79, "kg"),

                        # dose
                        "PODOSE_rap": Q_(5, "mg"),
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
                    count="count"
                ),
                observable=FitData(
                    self, task=f"task_rap_{intervention}", xid="time", yid=f"[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET if intervention == "RAP_T5" else ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE
                )
            )

        return mappings

    def figures(self) -> Dict[str, Figure]:

        fig1 = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig1.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)

        plots[0].set_yaxis(self.label_rap, unit=self.unit_rap, scale="linear")

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
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=intervention,
                color=self.colors[intervention],
            )
        return {
            fig1.sid: fig1
        }


if __name__ == "__main__":
    run_experiments(KorthBradley2012, output_dir=KorthBradley2012.__name__)
