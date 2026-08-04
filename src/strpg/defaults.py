"""Single source of truth for public STR-PG runtime defaults."""

from __future__ import annotations

from dataclasses import dataclass

from .alignment import AlignmentSoftmaxParameters
from .phmm import PHMMParameters

DEFAULT_LIKELIHOOD_BACKEND = "sw"
SUPPORTED_LIKELIHOOD_BACKENDS = ("sw", "phmm")
READ_CACHE_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class GenotypeDefaults:
    likelihood_backend: str = DEFAULT_LIKELIHOOD_BACKEND
    read_selection_backend_independent: bool = True
    use_population_frequency: bool = True
    use_balding_nichols: bool = True
    use_length_smoothing: bool = True
    use_mixture_model: bool = True


GENOTYPE_DEFAULTS = GenotypeDefaults()
SW_DEFAULTS = AlignmentSoftmaxParameters(
    match=2,
    mismatch=-5,
    gap=-5,
    band=500,
    temperature=12.0,
)
# These are the audited pre-calibration defaults. The development-grid P05
# parameters did not generalize and are intentionally not production defaults.
PHMM_DEFAULTS = PHMMParameters(
    mismatch_rate=0.01,
    flank_indel_open=0.003,
    repeat_indel_open=0.005,
    flank_indel_extend=0.10,
    repeat_indel_extend=0.30,
)
