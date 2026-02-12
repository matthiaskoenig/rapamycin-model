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


class Brattstrom2000(RapamycinSimulationExperiment):
    """Simulation experiment of Brattstrom2000."""

    interventions = ["RAP03", "RAP1", "RAP3", "RAP5", "RAP8"]
    colors = {
        "RAP03": "black",
        "RAP1": "#FFE299",
        "RAP3": "#FFB347",
        "RAP5": "#FF944D",
        "RAP8": "#D2691E",
    }
    doses = {  # [mg]
        "RAP03": 0.3,
        "RAP1": 1,
        "RAP3": 3,
        "RAP5": 5,
        "RAP8": 8,
    }
    bodyweights = {  # [mg]
        "RAP03": 74.2,
        "RAP1": 82.3,
        "RAP3": 74.5,
        "RAP5": 75.8,
        "RAP8": 74.8,
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
            dose = self.doses[intervention]
            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=350 * 60,  # [min]
                    steps=2000,
                    changes={
                        **self.default_changes(),
                        # physiological changes
                        "BW": Q_(self.bodyweights[intervention], "kg"),
                        # dose
                        "PODOSE_rap": Q_(dose, "mg"),
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
                    self, task=f"task_rap_{intervention}", xid="time",yid="[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE
                )
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",

        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)

        fig_12 = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (12h)",
        )
        plots_12 = fig_12.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots_12[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap),
        plots_12[0].set_xaxis(self.label_time, unit=self.unit_time,min=0, max=12)


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

            plots_12[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=intervention,
                color=self.colors[intervention],
            )
            plots_12[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=intervention,
                color=self.colors[intervention],
            )

        return {
            fig.sid: fig,
            fig_12.sid: fig_12,
        }


if __name__ == "__main__":
    run_experiments(Brattstrom2000, output_dir=Brattstrom2000.__name__)
