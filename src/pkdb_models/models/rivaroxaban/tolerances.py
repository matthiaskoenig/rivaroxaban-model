"""The tolerances of the integrator for the rivaroxaban model."""

from sbmlsim.model.tolerances import AbsoluteTolerance


def absolute_tolerance(value: float, species: float) -> AbsoluteTolerance:
    """Get the absolute tolerance of the integrator for the rivaroxaban model.

    The species of rivaroxaban and its metabolites are small: the plasma
    concentrations are about 1e-6 to 1e-3 mM (1e-6 mM after the microdose of
    25 µg of the cocktail studies Lenard2024, Lenard2025 and Rohr2024), the
    amounts about 1e-6 to 1e-1 mmole. The other states are the doses in mg (the
    targets of the rate rules).

    Args:
        value: tolerance of the other states [mg].
        species: tolerance of the species, of their concentrations [mM] and
            amounts [mmole].

    Returns:
        The absolute tolerance per kind of state.
    """
    return AbsoluteTolerance(amount=species, concentration=species, other=value)


#: the absolute tolerance of the simulation experiments
ABSOLUTE_TOLERANCE = absolute_tolerance(1e-10, species=1e-14)
