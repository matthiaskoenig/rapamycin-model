"""Simulation experiment of the study Zhang2017."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Genotype,
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


class Zhang2017(RapamycinSimulationExperiment):
    """Simulation experiment of Zhang2017.

    31 healthy Chinese Han men took a single 5 mg tablet of sirolimus after an
    overnight fast. The study reports the timecourses of the same 31 subjects
    twice, split by the CYP3A4 genotype (*1/*1 n=18, *1/*1G n=10, *1G/*1G n=3,
    Fig1A) and by the CYP3A5 genotype (*1/*1 n=2, *1/*3 n=14, *3/*3 n=15,
    Fig1B). A group of one split has the genotype of its split, the genotype of
    the other gene is the distribution of the cohort (the genotypes of the other
    gene within a group are not reported), i.e., its mean activity.
    """

    interventions: ClassVar[list[str]] = ["RAP5"]

    #: the subjects per genotype of the cohort (Table 2)
    counts: ClassVar[dict[str, int]] = {
        "CYP3A4_1_1": 18,
        "CYP3A4_1_1G": 10,
        "CYP3A4_1G_1G": 3,
        "CYP3A5_1_1": 2,
        "CYP3A5_1_3": 14,
        "CYP3A5_3_3": 15,
    }
    cyp3a4_groups: ClassVar[list[str]] = ["CYP3A4_1_1", "CYP3A4_1_1G", "CYP3A4_1G_1G"]
    cyp3a5_groups: ClassVar[list[str]] = ["CYP3A5_1_1", "CYP3A5_1_3", "CYP3A5_3_3"]
    groups: ClassVar[list[str]] = cyp3a4_groups + cyp3a5_groups
    genotypes: ClassVar[dict[str, str]] = {
        "CYP3A4_1_1": "*1/*1",
        "CYP3A4_1_1G": "*1/*1G",
        "CYP3A4_1G_1G": "*1G/*1G",
        "CYP3A5_1_1": "*1/*1",
        "CYP3A5_1_3": "*1/*3",
        "CYP3A5_3_3": "*3/*3",
    }
    colors: ClassVar[dict[str, str]] = {
        "CYP3A4_1_1": "black",
        "CYP3A4_1_1G": "tab:blue",
        "CYP3A4_1G_1G": "darkblue",
        "CYP3A5_1_1": "#CC3333",
        "CYP3A5_1_3": "#FF6666",
        "CYP3A5_3_3": "black",
    }
    bodyweights: ClassVar[dict[str, float]] = {  # [kg]
        "CYP3A4_1_1": 63.6,
        "CYP3A4_1_1G": 64.3,
        "CYP3A4_1G_1G": 72.3,
        "CYP3A5_1_1": 74.0,
        "CYP3A5_1_3": 65.1,
        "CYP3A5_3_3": 63.1,
    }
    heights: ClassVar[dict[str, float]] = {  # [cm]
        "CYP3A4_1_1": 171.2,
        "CYP3A4_1_1G": 172.7,
        "CYP3A4_1G_1G": 180.0,
        "CYP3A5_1_1": 181.0,
        "CYP3A5_1_3": 172.5,
        "CYP3A5_3_3": 171.5,
    }

    @classmethod
    def mean_activity(cls, groups: list[str], activity: dict[str, float]) -> float:
        """Get the mean activity of an enzyme over the genotypes of the cohort.

        Args:
            groups: the groups of the genotypes of the enzyme.
            activity: the activity per genotype of the enzyme.

        Returns:
            The activity weighted by the subjects per genotype.
        """
        n = sum(cls.counts[group] for group in groups)
        return (
            sum(cls.counts[group] * activity[cls.genotypes[group]] for group in groups)
            / n
        )

    def f_cyp3a(self, group: str) -> tuple[float, float]:
        """Get the activity of CYP3A4 and CYP3A5 of a group.

        Args:
            group: the group.

        Returns:
            The activity of CYP3A4 and of CYP3A5, the genotype of the group for
            its enzyme and the mean of the cohort for the other enzyme.
        """
        if group in self.cyp3a4_groups:
            return (
                self.cyp3a4_activity[self.genotypes[group]],
                self.mean_activity(self.cyp3a5_groups, self.cyp3a5_activity),
            )
        return (
            self.mean_activity(self.cyp3a4_groups, self.cyp3a4_activity),
            self.cyp3a5_activity[self.genotypes[group]],
        )

    def group_label(self, group: str) -> str:
        """Get the label of a group, its genotype."""
        return f"{group.split('_')[0]} {self.genotypes[group]}"

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1A", "Fig1B"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("rapamycin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for group in self.groups:
            f_cyp3a4, f_cyp3a5 = self.f_cyp3a(group)
            simulations[f"rap_RAP5_{group}"] = Simulation(
                time_unit="hr",
                end=150,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[group], "kg"),
                    "HEIGHT": Q(self.heights[group], "cm"),
                    # dose
                    "PODOSE_rap": Q(5, "mg"),
                    # fasted
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    # CYP3A4 and CYP3A5 activity
                    "GU__f_cyp3a4": Q(f_cyp3a4, "dimensionless"),
                    "LI__f_cyp3a4": Q(f_cyp3a4, "dimensionless"),
                    "GU__f_cyp3a5": Q(f_cyp3a5, "dimensionless"),
                    "LI__f_cyp3a5": Q(f_cyp3a5, "dimensionless"),
                },
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:

        mappings = {}
        for group in self.groups:
            # the genotype of the split of the group
            genotype = Genotype(f"{group.split('_')[0]} {self.genotypes[group]}")

            mappings[f"fm_rap_RAP5_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"rapamycin_RAP5_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_rap_RAP5_{group}",
                    xid="time",
                    yid="[Cveblood_rap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE,
                    genotype=genotype,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        figs = {}
        for gene, groups in [
            ("CYP3A4", self.cyp3a4_groups),
            ("CYP3A5", self.cyp3a5_groups),
        ]:
            fig = Figure(
                experiment=self,
                sid=f"Fig1_{gene}",
                name=f"{self.__class__.__name__} (5 mg, {gene})",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

            for group in groups:
                # simulation
                plots[0].add_data(
                    task=f"task_rap_RAP5_{group}",
                    xid="time",
                    yid="[Cveblood_rap]",
                    label=self.group_label(group),
                    color=self.colors[group],
                )
                # data
                plots[0].add_data(
                    dataset=f"rapamycin_RAP5_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.group_label(group),
                    color=self.colors[group],
                )
            figs |= self.figures_by_sid(fig)

        return figs


if __name__ == "__main__":
    MODEL.run_experiments(Zhang2017, output_dir=Zhang2017.__name__)
