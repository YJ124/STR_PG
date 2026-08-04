"""Public likelihood backend factory."""

from __future__ import annotations

from ..alignment import AlignmentSoftmaxParameters
from ..phmm import PHMMParameters
from .base import AlleleLikelihoodBackend
from .phmm import PHMMLikelihoodBackend
from .sw import SWLikelihoodBackend


def create_likelihood_backend(
    name: str,
    *,
    sw_params: AlignmentSoftmaxParameters | None = None,
    phmm_params: PHMMParameters | None = None,
) -> AlleleLikelihoodBackend:
    normalized = name.strip().lower()
    if normalized == "sw":
        return SWLikelihoodBackend(sw_params)
    if normalized == "phmm":
        return PHMMLikelihoodBackend(phmm_params)
    raise ValueError(
        f"Unsupported likelihood backend {name!r}; expected 'sw' or 'phmm'"
    )


__all__ = [
    "AlleleLikelihoodBackend",
    "SWLikelihoodBackend",
    "PHMMLikelihoodBackend",
    "create_likelihood_backend",
]
