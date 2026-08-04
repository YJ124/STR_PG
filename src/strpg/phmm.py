from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - exact Python fallback remains available
    njit = None

from .models import AlleleRecord
from .utils import reverse_complement

NEG_INF = -math.inf


def _logadd2(a: float, b: float) -> float:
    if a == NEG_INF:
        return b
    if b == NEG_INF:
        return a
    if a < b:
        a, b = b, a
    return a + math.log1p(math.exp(b - a))


def _logadd3(a: float, b: float, c: float) -> float:
    return _logadd2(_logadd2(a, b), c)


if njit is not None:
    _logadd2_numba = njit(cache=True)(_logadd2)

    @njit(cache=True)
    def _logadd3_numba(a: float, b: float, c: float) -> float:
        return _logadd2_numba(_logadd2_numba(a, b), c)
else:  # pragma: no cover
    _logadd2_numba = _logadd2
    _logadd3_numba = _logadd3


def _logsum(values: list[float]) -> float:
    total = NEG_INF
    for value in values:
        total = _logadd2(total, value)
    return total


def _forward_array(
    read: np.ndarray,
    template: np.ndarray,
    repeat_start: int,
    repeat_end: int,
    flank_logs: tuple[float, float, float, float, float],
    repeat_logs: tuple[float, float, float, float, float],
    match_log: float,
    mismatch_log: float,
    ins_emit_log: float,
) -> float:
    n, m = len(read), len(template)
    if n == 0 or m == 0:
        return NEG_INF
    start_log = -math.log(m + 1.0)
    m_prev = np.full(m + 1, start_log, dtype=np.float64)
    i_prev = np.full(m + 1, NEG_INF, dtype=np.float64)
    d_prev = np.full(m + 1, NEG_INF, dtype=np.float64)
    for i in range(n):
        m_cur = np.full(m + 1, NEG_INF, dtype=np.float64)
        i_cur = np.full(m + 1, NEG_INF, dtype=np.float64)
        d_cur = np.full(m + 1, NEG_INF, dtype=np.float64)
        for j in range(1, m + 1):
            logs = repeat_logs if repeat_start <= j - 1 < repeat_end else flank_logs
            log_mm, log_mi, log_md, log_gap_to_m, log_gap_ext = logs
            emit = match_log if read[i] == template[j - 1] else mismatch_log
            m_cur[j] = emit + _logadd3_numba(
                m_prev[j - 1] + log_mm,
                i_prev[j - 1] + log_gap_to_m,
                d_prev[j - 1] + log_gap_to_m,
            )
            i_cur[j] = ins_emit_log + _logadd2_numba(
                m_prev[j] + log_mi, i_prev[j] + log_gap_ext
            )
            d_cur[j] = _logadd2_numba(
                m_cur[j - 1] + log_md, d_cur[j - 1] + log_gap_ext
            )
        m_prev, i_prev, d_prev = m_cur, i_cur, d_cur
    total = NEG_INF
    for state in (m_prev, i_prev, d_prev):
        for j in range(1, m + 1):
            total = _logadd2_numba(total, state[j])
    return total


_forward_array_jit = njit(cache=True)(_forward_array) if njit is not None else None


@dataclass(frozen=True)
class PHMMParameters:
    mismatch_rate: float = 0.01
    flank_indel_open: float = 0.003
    repeat_indel_open: float = 0.005
    flank_indel_extend: float = 0.10
    repeat_indel_extend: float = 0.30
    insertion_background: float = 0.25

    def validate(self) -> None:
        for name, value in self.__dict__.items():
            if not (0.0 < value < 1.0):
                raise ValueError(f"{name} must be in (0,1), got {value}")
        if 2 * self.repeat_indel_open >= 1 or 2 * self.flank_indel_open >= 1:
            raise ValueError("Twice the indel-open probability must be <1")


class MotifAwarePHMM:
    """Semi-global motif-aware profile/pair HMM evaluated by forward DP.

    The registered allele sequence instantiates the motif-phase-aware repeat profile.
    Match, insertion and deletion transitions are position-specific: repeat-body
    positions receive higher indel probabilities than retained flanking anchors.
    """

    def __init__(self, params: PHMMParameters | None = None):
        self.params = params or PHMMParameters()
        self.params.validate()
        p = self.params
        self._flank_logs = (
            math.log(1.0 - 2.0 * p.flank_indel_open),
            math.log(p.flank_indel_open),
            math.log(p.flank_indel_open),
            math.log(1.0 - p.flank_indel_extend),
            math.log(p.flank_indel_extend),
        )
        self._repeat_logs = (
            math.log(1.0 - 2.0 * p.repeat_indel_open),
            math.log(p.repeat_indel_open),
            math.log(p.repeat_indel_open),
            math.log(1.0 - p.repeat_indel_extend),
            math.log(p.repeat_indel_extend),
        )
        self._match_log = math.log(1.0 - p.mismatch_rate)
        self._mismatch_log = math.log(p.mismatch_rate / 3.0)
        self._ins_emit_log = math.log(p.insertion_background)

    def _forward_log_likelihood_python(self, read: str, allele: AlleleRecord) -> float:
        read = read.upper()
        template = allele.sequence.upper()
        n, m = len(read), len(template)
        if n == 0 or m == 0:
            return NEG_INF

        start_log = -math.log(m + 1.0)
        m_prev = [start_log] * (m + 1)
        i_prev = [NEG_INF] * (m + 1)
        d_prev = [NEG_INF] * (m + 1)

        for read_base in read:
            m_cur = [NEG_INF] * (m + 1)
            i_cur = [NEG_INF] * (m + 1)
            d_cur = [NEG_INF] * (m + 1)
            for j in range(1, m + 1):
                logs = self._repeat_logs if allele.repeat_start <= j - 1 < allele.repeat_end else self._flank_logs
                log_mm, log_mi, log_md, log_gap_to_m, log_gap_ext = logs
                emit = self._match_log if read_base == template[j - 1] else self._mismatch_log
                m_cur[j] = emit + _logadd3(
                    m_prev[j - 1] + log_mm,
                    i_prev[j - 1] + log_gap_to_m,
                    d_prev[j - 1] + log_gap_to_m,
                )
                i_cur[j] = self._ins_emit_log + _logadd2(m_prev[j] + log_mi, i_prev[j] + log_gap_ext)
                d_cur[j] = _logadd2(m_cur[j - 1] + log_md, d_cur[j - 1] + log_gap_ext)
            m_prev, i_prev, d_prev = m_cur, i_cur, d_cur

        # The uniform start distribution above already assigns probability
        # mass across alternative alignment starts. Terminal positions are
        # mutually exclusive paths and therefore must be summed, not averaged.
        # Dividing this sum by m applies a second template-length prior and
        # creates a non-generative allele-length bias.
        return _logsum(m_prev[1:] + i_prev[1:] + d_prev[1:])

    def forward_log_likelihood(self, read: str, allele: AlleleRecord) -> float:
        if _forward_array_jit is None:
            return self._forward_log_likelihood_python(read, allele)
        read_array = np.frombuffer(read.upper().encode("ascii"), dtype=np.uint8)
        template_array = np.frombuffer(
            allele.sequence.upper().encode("ascii"), dtype=np.uint8
        )
        return float(
            _forward_array_jit(
                read_array,
                template_array,
                allele.repeat_start,
                allele.repeat_end,
                self._flank_logs,
                self._repeat_logs,
                self._match_log,
                self._mismatch_log,
                self._ins_emit_log,
            )
        )

    def orientation_marginal_log_likelihood(self, read: str, allele: AlleleRecord) -> float:
        fwd = self.forward_log_likelihood(read, allele)
        rev = self.forward_log_likelihood(reverse_complement(read), allele)
        return _logadd2(fwd, rev) - math.log(2.0)
