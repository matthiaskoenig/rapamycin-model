from typing import Dict
import matplotlib.pyplot as plt
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


class Leelahavanichkul2005(RapamycinSimulationExperiment):
    """Simulation experiment of Leelahavanichkul2005."""

    interventions = ["RAP6"]
    colors = {
        "RAP6": "black",
    }
    subjects = [f"S{k+1}" for k in range(12)]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab3A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("rapamycin_"):
                    if fig_id == "Fig1":
                        dset.unit_conversion("value", 1 / self.Mr.rap)
                    elif fig_id == "Tab3A":
                        dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[label] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=30 * 60,  # [min]
                    steps=1000,
                    changes={

                        **self.default_changes(),

                        # physiological changes
                        "BW": Q_(64.73, "kg"),  # weighted average men and women

                        # dose
                        "PODOSE_rap": Q_(6, "mg"),
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
                    yid_sd=None,
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
                    coadministration=Coadministration.NONE
                )
            )

            for ks, subject in enumerate(self.subjects):
                mappings[f"fm_rap_{intervention}_{subject}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"rapamycin_{intervention}_{subject}",
                        xid="time",
                        yid="value",
                        yid_sd=None,
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
                        coadministration=Coadministration.NONE
                    )
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        figures = {}

        cmap = plt.get_cmap("copper")
        subject_colors = {subj: cmap(i / (len(self.subjects)-1)) for i, subj in enumerate(self.subjects)}

        # individual data
        fig = Figure(
            experiment=self,
            sid=f"Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        for ks, subject in enumerate(self.subjects):
            plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

            for intervention in self.interventions:
                # simulation
                plots[0].add_data(
                    task=f"task_rap_{intervention}",
                    xid="time",
                    yid=f"[Cveblood_rap]",
                    label=intervention if ks == 0 else None,
                    color=subject_colors[subject],
                )
                # data
                plots[0].add_data(
                    dataset=f"rapamycin_{intervention}_{subject}",
                    xid="time",
                    yid="value",
                    count="count",
                    yid_sd=None,
                    label=subject,
                    color=subject_colors[subject],
                )

            figures[fig.sid] = fig

        # mean data
        fig3a = Figure(
            experiment=self,
            sid="Tab3A",
            name=f"{self.__class__.__name__}",
        )
        plots3a = fig3a.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )

        plots3a[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)

        for intervention in self.interventions:
            # simulation
            plots3a[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=intervention,
                color=self.colors[intervention],
            )
            plots3a[0].add_data(
                dataset="rapamycin_RAP6",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=intervention,
                color=self.colors[intervention],
            )

        figures[fig3a.sid] = fig3a

        return figures



if __name__ == "__main__":
    run_experiments(Leelahavanichkul2005, output_dir=Leelahavanichkul2005.__name__)
