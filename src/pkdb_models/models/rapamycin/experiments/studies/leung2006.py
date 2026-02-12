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


class Leung2006(RapamycinSimulationExperiment):
    """Simulation experiment of Leung2006."""

    interventions = ["RAP40"]

    sids = {
        "blood": "[Cveblood_rap]",
        "plasma": "[Cve_rap]",
    }
    tissues = list(sids.keys())

    colors = {
        "blood": "tab:red",
        "plasma": "tab:blue",
    }
    markers = {
        "Fig1": "s",
        "Fig2": "o",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for tissue in self.tissues:
                df_sub = df[(df["measurement_type"] == "concentration") & (df["tissue"] == tissue)]
                dset = DataSet.from_df(df_sub.copy(), self.ureg)
                if tissue in ("blood", "plasma"):
                    dset.unit_conversion("mean", 1/self.Mr.rap)
                dsets[f"{fig_id}_rapamycin_RAP40_{tissue}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims: Dict[str, TimecourseSim] = {}

        for intervention in self.interventions:
            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=150 * 60,  # [min]
                    steps=2000,
                    changes={

                        **self.default_changes(),

                        # physiological changes
                        "BW": Q_(73, "kg"),  # weighted average men and women

                        # dose
                        "PODOSE_rap": Q_(40, "mg"),
                    },
                )]
            )

        #console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}

        for intervention in self.interventions:
            for tissue in self.tissues:
                for fig_id in ["Fig2"]:  # ["Fig1", "Fig2"]:
                    mappings[f"fm_rap_{intervention}_{tissue}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{fig_id}_rapamycin_RAP40_{tissue}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self, task=f"task_rap_{intervention}", xid="time", yid=self.sids[tissue],
                        ),
                        metadata=RapamycinMappingMetaData(
                            tissue=Tissue.BLOOD if tissue == "blood" else Tissue.PLASMA,
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
            sid="Fig1_2",
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time, max=48), legend=True)
        plots[0].set_yaxis(self.label_rap, unit=self.unit_rap)

        for intervention in self.interventions:
            for tissue in self.tissues:

                # simulation
                plots[0].add_data(
                    task=f"task_rap_{intervention}",
                    xid="time",
                    yid=self.sids[tissue],
                    label=tissue,
                    color=self.colors[tissue],
                )

                # data
                for fig_id in ["Fig2"]:  # ["Fig1", "Fig2"]:
                    plots[0].add_data(
                        dataset=f"{fig_id}_rapamycin_RAP40_{tissue}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=f"{tissue} ({fig_id})",
                        color=self.colors[tissue],
                        marker=self.markers[fig_id],
                    )
        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Leung2006, output_dir=Leung2006.__name__)
