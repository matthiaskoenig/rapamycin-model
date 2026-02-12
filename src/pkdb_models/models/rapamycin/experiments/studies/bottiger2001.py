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


class Bottiger2001(RapamycinSimulationExperiment):
    """Simulation experiment of Bottiger2001.

    Diltiazem inhibits the human CYP3A4 enzyme that metabolizes
    cyclosporine in the intestines, as well as in the liver.
    """

    interventions = ["RAP10", "RAP10_DIL120"]
    colors = {
        "RAP10": "black",
        "RAP10_DIL120": "#e3256b",
    }

    f_cyp3a4 = {
        "RAP10": 1.0,  # no change in activity
        "RAP10_DIL120": 0.5,  # inhibition CYP3A4 by diltiazem
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                if label.startswith("rapamycin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[label] = dset

        # console.print(dsets.keys())
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
                        "BW": Q_((82*12 + 64*6)/18, "kg"),  # weighted average
                        "PODOSE_rap": Q_(10, "mg"),

                        # CYP3A4 activity (inhibition by diltiazem)
                        "GU__f_cyp3a4": Q_(self.f_cyp3a4[intervention], "dimensionless"),
                        "LI__f_cyp3a4": Q_(self.f_cyp3a4[intervention], "dimensionless"),
                    },
                )]
            )
        # console.print(tcsims.keys())
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
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE if intervention == "RAP10" else Coadministration.DILTIAZEM
                )
            )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit = self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for intervention in self.interventions:

            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=intervention,
                color=self.colors[intervention],
            )
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
        fig.sid: fig,
    }


if __name__ == "__main__":
    run_experiments(Bottiger2001, output_dir=Bottiger2001.__name__)
