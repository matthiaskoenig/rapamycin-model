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


class BasaDenes2019(RapamycinSimulationExperiment):
    """BasaDenes2019 experiment.

    Simulation experiment of BasaDenes2019. Co-administration of famotidine, no effect on CYPs
    or p-gp, but it enhanced sirolimus stability and absorption,
    reducing risk of degradation by gastric acid in the stomach.
    """

    interventions = ["RAP05", "RAP2", "RAP10", "RAP40", "RAP40_FD"]
    colors = {
        "RAP05": "#FFB347",
        "RAP2": "#FF944D",
        "RAP10": "#D2691E",
        "RAP40": "black",
        "RAP40_FD": "purple",
    }
    doses = {
        "RAP05": 0.5,
        "RAP2": 2,
        "RAP10": 10,
        "RAP40": 40,
        "RAP40_FD": 40,
    }
    f_oatp = {
        "RAP05": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP2": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP10": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP40": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP40_FD": RapamycinSimulationExperiment.fasting_map["fed"],
    }
    intervention_subjects = {
        "RAP05":[f"S{i}" for i in range(1,9)],
        "RAP2":[f"A{i}" for i in range(1,9)],
        "RAP10": [f"B{i}" for i in range(1,9)],
        "RAP40":[f"C{i}" for i in range(1,9)],
        "RAP40_FD":[f"D{i}" for i in range(1,9)],
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1A", "Fig1C"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("rapamycin_"):
                    if fig_id == "Fig1A":
                        dset.unit_conversion("mean", 1 / self.Mr.rap)
                    elif fig_id == "Fig1C":
                        dset.unit_conversion("value", 1 / self.Mr.rap)

                dsets[label] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            dose = self.doses[intervention]
            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=51 * 60,  # [min]
                    steps=2000,
                    changes={

                        **self.default_changes(),

                        # physiological changes
                        "BW": Q_(85.03, "kg"),  # weighted average men and women

                        # dose
                        "PODOSE_rap": Q_(dose, "mg"),

                        # fasting/fed
                        "GU__f_oatp": Q_(self.f_oatp[intervention], "dimensionless"),
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
                    self, task=f"task_rap_{intervention}", xid="time", yid="[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FED if intervention == "RAP40_FD" else Fasting.FASTED,
                    coadministration=Coadministration.FAMOTIDINE
                )
            )

            for subject in self.intervention_subjects[intervention]:
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
                        self, task=f"task_rap_{intervention}", xid="time", yid="[Cveblood_rap]",
                    ),
                    metadata=RapamycinMappingMetaData(
                        tissue=Tissue.BLOOD,
                        route=Route.PO,
                        application_form=ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED if intervention == "RAP40_FD" else Fasting.FASTED,
                        coadministration=Coadministration.FAMOTIDINE
                    )
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        figures = {}

        all_subjects = [s for subjects in self.intervention_subjects.values() for s in subjects]
        subject_colors = {subj: plt.get_cmap("tab10")(i % 10) for i, subj in enumerate(all_subjects)}

        for intervention in self.interventions:
            fig = Figure(
                experiment=self,
                sid=f"Fig1C_{intervention}",
                name=f"{self.__class__.__name__}_{intervention}",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)
                # color = subject_colors[subject]

            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=intervention,
                color=self.colors[intervention],
                )
            for ks, subject in enumerate(self.intervention_subjects[intervention]):
                    # data
                plots[0].add_data(
                    dataset=f"rapamycin_{intervention}_{subject}",
                    xid="time",
                    yid="value",
                    count="count",
                    yid_sd=None,
                    label=subject,
                    color=self.colors[intervention],
                )

            figures[fig.sid] = fig

        # mean data (dose dependency)
        fig1a = Figure(
            experiment=self,
            sid="Fig1A_dose",
            name=f"{self.__class__.__name__} (dose)",
        )
        plots1a = fig1a.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots1a[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)

        for intervention in self.interventions:
            if intervention.endswith("_FD"):
                continue

            # simulation
            plots1a[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=intervention,
                color=self.colors[intervention],
            )
            plots1a[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=intervention,
                color=self.colors[intervention],
            )
        figures[fig1a.sid] = fig1a

        # mean data (fasted vs. food)
        fig1a = Figure(
            experiment=self,
            sid="Fig1A_food",
            name=f"{self.__class__.__name__} (food)",
        )
        plots1a = fig1a.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time))
        plots1a[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)

        for intervention in self.interventions:
            if not intervention.startswith("RAP40"):
                continue
            color = "purple" if intervention == "RAP40_FD" else "black"
            label = "fed" if intervention == "RAP40_FD" else "fasted"

            # simulation
            plots1a[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=label,
                color=color,
            )
            plots1a[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color=color,
            )
        figures[fig1a.sid] = fig1a

        return figures



if __name__ == "__main__":
    run_experiments(BasaDenes2019, output_dir=BasaDenes2019.__name__)
