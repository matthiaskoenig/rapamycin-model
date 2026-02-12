from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData

from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)
from pkdb_models.models.rapamycin.experiments.metadata import (
    RapamycinMappingMetaData,
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, Coadministration,
    Genotype
)

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rapamycin.helpers import run_experiments

class Renders2007(RapamycinSimulationExperiment):
    """Simulation experiment of Renders2007."""


    groups = ["CYP3A5_1_1", "CYP3A5_1_3", "CYP3A5_3_3"]
    colors = {
        "CYP3A5_1_1": "black",
        "CYP3A5_1_3": "#CC3333",
        "CYP3A5_3_3": "#FF6666",

    }
    f_cyp3a5 = {
        "CYP3A5_1_1": RapamycinSimulationExperiment.cyp3a5_activity["*1/*1"],
        "CYP3A5_1_3": RapamycinSimulationExperiment.cyp3a5_activity["*1/*3"],
        "CYP3A5_3_3": RapamycinSimulationExperiment.cyp3a5_activity["*3/*3"],
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig4"]:
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
        for group in self.groups:
            tcsims[f"rap_{group}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=25 * 60,  # [min]
                    steps=2000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(70, "kg"),  # FIXME: no bodyweight in article
                        "PODOSE_rap": Q_(2.9, "mg"),

                        # CYP3A5 activity
                        "GU__f_cyp3a5": Q_(self.f_cyp3a5[group], "dimensionless"),
                        "LI__f_cyp3a5": Q_(self.f_cyp3a5[group], "dimensionless"),
                    },
                )]
            )

        #console.print(tcsims.keys())
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}

        for group in self.groups:
            if group == "CYP3A5_1_1":
                continue

            genotype_tokens = group.split("_")
            genotype = Genotype(f"{genotype_tokens[0]} *{genotype_tokens[1]}/*{genotype_tokens[2]}")

            mappings[f"fm_rap_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"rapamycin_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_rap_{group}", xid="time", yid="[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.NR,
                    dosing=Dosing.SINGLE,
                    health=Health.RENAL_TRANSPLANT,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.NONE,
                    genotype=genotype,
                )
            )

        return mappings


    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig4",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True
        )
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for group in self.groups:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{group}",
                xid="time",
                yid="[Cveblood_rap]",
                label=group,
                color=self.colors[group],
            )
            # data
            if group == "CYP3A5_1_1":
                continue
            plots[0].add_data(
                dataset=f"rapamycin_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=group,
                color=self.colors[group],
                )

        return {fig.sid: fig}


if __name__ == "__main__":
    run_experiments(Renders2007, output_dir=Renders2007.__name__)
