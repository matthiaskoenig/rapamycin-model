"""Simulation experiment of the study BasaDenes2019."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.rapamycin import MODEL
from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)


class BasaDenes2019(RapamycinSimulationExperiment):
    """BasaDenes2019 experiment.

    Healthy subjects (4 cohorts of 8) took a single dose of 0.5, 2, 10 or 40 mg
    of a novel nano-amorphous formulation of sirolimus (a powder reconstituted
    in water, not Rapamune) after an overnight fast; the 40 mg cohort took the
    dose again after a high fat meal (crossover). All took 40 mg famotidine 3 h
    before the dose, to raise the gastric pH: famotidine does not inhibit CYP3A
    or P-glycoprotein, it reduces the degradation of sirolimus by gastric acid.
    The formulation has a higher bioavailability than Rapamune (Cmax per dose
    about twice the one of the solution) and its absorption is less than dose
    proportional from 10 to 40 mg, neither of which the model describes: a
    limitation of the model. The means and SD (Fig1A) are mapped, the
    individual curves of the same subjects (Fig1C, 0-24 h) are plotted.
    """

    interventions: ClassVar[list[str]] = ["RAP05", "RAP2", "RAP10", "RAP40", "RAP40_FD"]
    colors: ClassVar[dict[str, str]] = {
        "RAP05": "#FFB347",
        "RAP2": "#FF944D",
        "RAP10": "#D2691E",
        "RAP40": "black",
        "RAP40_FD": "purple",
    }
    labels_intervention: ClassVar[dict[str, str]] = {
        "RAP05": "0.5 mg",
        "RAP2": "2 mg",
        "RAP10": "10 mg",
        "RAP40": "40 mg",
        "RAP40_FD": "40 mg fed",
    }
    doses: ClassVar[dict[str, float]] = {
        "RAP05": 0.5,
        "RAP2": 2,
        "RAP10": 10,
        "RAP40": 40,
        "RAP40_FD": 40,
    }
    #: body weight [kg] and height [cm] of the cohorts (Table 1), the fed arm
    #: is the 40 mg cohort
    bodyweights: ClassVar[dict[str, float]] = {
        "RAP05": 79.75,
        "RAP2": 87.25,
        "RAP10": 86.95,
        "RAP40": 86.18,
        "RAP40_FD": 86.18,
    }
    heights: ClassVar[dict[str, float]] = {
        "RAP05": 169.5,
        "RAP2": 177.8,
        "RAP10": 174.0,
        "RAP40": 174.5,
        "RAP40_FD": 174.5,
    }
    #: fasted, the 40 mg cohort a second time after a high fat meal
    fed: ClassVar[dict[str, float]] = {
        "RAP05": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP2": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP10": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP40": RapamycinSimulationExperiment.fasting_map["fasted"],
        "RAP40_FD": RapamycinSimulationExperiment.fasting_map["fed"],
    }
    #: the subjects of the individual curves in the data (S: 0.5 mg, A: 2 mg,
    #: B: 10 mg, C: 40 mg, D: 40 mg fed)
    intervention_subjects: ClassVar[dict[str, list[str]]] = {
        "RAP05": [f"S{i}" for i in range(1, 9)],
        "RAP2": [f"A{i}" for i in range(1, 9)],
        "RAP10": [f"B{i}" for i in range(1, 9)],
        "RAP40": [f"C{i}" for i in range(1, 9)],
        "RAP40_FD": [f"D{i}" for i in range(1, 9)],
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1A", "Fig1C"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("rapamycin_"):
                    if fig_id == "Fig1A":
                        dset.unit_conversion("mean", 1 / self.Mr.rap)
                    elif fig_id == "Fig1C":
                        dset.unit_conversion("value", 1 / self.Mr.rap)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=51,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[intervention], "kg"),
                    "HEIGHT": Q(self.heights[intervention], "cm"),
                    # dose
                    "PODOSE_rap": Q(self.doses[intervention], "mg"),
                    # fasted or fed
                    "GU__fed": Q(self.fed[intervention], "dimensionless"),
                },
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:

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
                    self,
                    task=f"task_rap_{intervention}",
                    xid="time",
                    yid="[Cveblood_rap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    # nano-amorphous sirolimus reconstituted in water
                    application_form=ApplicationForm.NANOPARTICLES,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FED
                    if intervention == "RAP40_FD"
                    else Fasting.FASTED,
                    coadministration=Coadministration.FAMOTIDINE,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        figures = {}

        # individual data (Fig1C)
        fig1c = Figure(
            experiment=self,
            sid="Fig1C",
            num_cols=len(self.interventions),
            name=f"{self.__class__.__name__} (individuals)",
        )
        plots1c = fig1c.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=0, max=25),
            legend=True,
        )
        for kp, intervention in enumerate(self.interventions):
            plots1c[kp].set_title(self.labels_intervention[intervention])
            plots1c[kp].set_yaxis(
                self.label_rap_blood, unit=self.unit_rap, scale="log", min=0.1
            )
            plots1c[kp].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.labels_intervention[intervention],
                color=self.colors[intervention],
            )
            for subject in self.intervention_subjects[intervention]:
                plots1c[kp].add_data(
                    dataset=f"rapamycin_{intervention}_{subject}",
                    xid="time",
                    yid="value",
                    count="count",
                    yid_sd=None,
                    label=subject,
                    color=self.colors[intervention],
                )
        figures |= self.figures_by_sid(fig1c)

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
        plots1a[0].set_yaxis(
            self.label_rap_blood, unit=self.unit_rap, scale="log", min=0.1
        )

        for intervention in self.interventions:
            if intervention.endswith("_FD"):
                continue

            # simulation
            plots1a[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.labels_intervention[intervention],
                color=self.colors[intervention],
            )
            plots1a[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.labels_intervention[intervention],
                color=self.colors[intervention],
            )
        figures |= self.figures_by_sid(fig1a)

        # mean data (fasted vs. fed)
        fig1a = Figure(
            experiment=self,
            sid="Fig1A_food",
            name=f"{self.__class__.__name__} (food)",
        )
        plots1a = fig1a.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots1a[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)

        for intervention in ["RAP40", "RAP40_FD"]:
            # simulation
            plots1a[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.labels_intervention[intervention],
                color=self.colors[intervention],
            )
            plots1a[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.labels_intervention[intervention],
                color=self.colors[intervention],
            )
        figures |= self.figures_by_sid(fig1a)

        return figures


if __name__ == "__main__":
    MODEL.run_experiments(BasaDenes2019, output_dir=BasaDenes2019.__name__)
