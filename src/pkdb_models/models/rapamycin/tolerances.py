"""The tolerances of the integrator for the rapamycin model."""

from sbmlsim.model.tolerances import AbsoluteTolerance


def absolute_tolerance(value: float, drug: float) -> AbsoluteTolerance:
    """Get the absolute tolerance of the integrator for the rapamycin model.

    The species of rapamycin and its metabolites have whole blood
    concentrations of about 1e-10 mM (0.1 ng/ml) in the terminal phase, which
    the data reach after several days; the doses (the targets of the rate
    rules) are milligrams.

    Args:
        value: tolerance of the targets of the rate rules, i.e. the doses.
        drug: tolerance of the species, of the amount and of the concentration
            [mmole, mM]; the tolerance of a concentration is multiplied by the
            volume of its compartment.

    Returns:
        The absolute tolerance per state.
    """
    return AbsoluteTolerance(amount=drug, concentration=drug, other=value)


#: the absolute tolerance of the simulation experiments
ABSOLUTE_TOLERANCE = absolute_tolerance(1e-10, drug=1e-14)
