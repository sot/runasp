"""Smoothers reproducing the upstream C pipeline for `process_fids.py`.

- `sg_smooth_per_slot` mirrors `aca_filter_centr` / `asp_centr_smooth`:
  per-slot Savitzky-Golay (polyorder 4, window from t_smooth / median dt).
- `hann_smooth` mirrors `aca_corr_fid` / `asp_slide_smooth`: normalized
  Hann convolution with replicated-endpoint padding.
"""

import numpy as np
from scipy.signal import savgol_filter
from scipy.signal.windows import hann


def _odd(n):
    n = int(np.ceil(n))
    return n + 1 if n % 2 == 0 else n


def sg_smooth_per_slot(*, time, slot, ang_y, ang_z, t_smooth=40.0, polyorder=4):
    """Savitzky-Golay smooth ang_y / ang_z per slot.

    Window length per slot is `ceil_to_odd(t_smooth / median(dt))`, capped
    so it never exceeds the slot's sample count (must also stay > polyorder).
    Edge handling is `mode="nearest"`, which matches the C tool's
    nearest-good-sample replacement at endpoints.
    """
    time = np.asarray(time, dtype=float)
    slot = np.asarray(slot)
    y_sm = np.asarray(ang_y, dtype=float).copy()
    z_sm = np.asarray(ang_z, dtype=float).copy()

    for s in np.unique(slot):
        sel = slot == s
        t_s = time[sel]
        if t_s.size <= polyorder + 1:
            continue
        dt = np.median(np.diff(t_s))
        window = min(_odd(t_smooth / dt), _odd(t_s.size - 1))
        if window <= polyorder:
            continue
        y_sm[sel] = savgol_filter(y_sm[sel], window, polyorder, mode="nearest")
        z_sm[sel] = savgol_filter(z_sm[sel], window, polyorder, mode="nearest")

    return y_sm, z_sm


def hann_smooth(values, *, windowlen=152):
    """Hann-windowed smooth with replicated-endpoint padding."""
    values = np.asarray(values, dtype=float)
    kernel = hann(windowlen, sym=True)
    kernel /= kernel.sum()
    pad = windowlen // 2
    padded = np.pad(values, pad, mode="edge")
    return np.convolve(padded, kernel, mode="same")[pad:pad + values.size]
