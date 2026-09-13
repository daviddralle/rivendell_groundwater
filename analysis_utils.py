"""Shared, explicitly bounded operations used by the manuscript notebooks."""
import re
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid


def interpolate_short_gaps(frame, max_missing):
    """Fill only complete interior missing runs of at most max_missing samples."""
    result = frame.copy()
    for name in frame:
        missing = frame[name].isna()
        groups = missing.ne(missing.shift()).cumsum()
        run_length = missing.groupby(groups).transform('sum')
        candidate = frame[name].interpolate(method='time', limit_area='inside')
        allowed = missing & (run_length <= max_missing)
        result.loc[allowed, name] = candidate.loc[allowed]
    return result


def integrate_profile(depth, theta, lower, upper, max_gap=0.61):
    """Integrate a piecewise-linear profile on fixed bounds; never extrapolate."""
    depth, theta = np.asarray(depth), np.asarray(theta)
    if len(depth) < 2 or np.any(np.diff(depth) <= 0):
        raise ValueError('Profile depths must be unique and increasing')
    if not np.isfinite(depth).all() or not np.isfinite(theta).all():
        raise ValueError('Profile contains nonfinite values')
    if depth[0] > lower or depth[-1] < upper:
        raise ValueError('Profile does not span the fixed interval')
    relevant = (depth[:-1] < upper) & (depth[1:] > lower)
    if (np.diff(depth)[relevant] > max_gap + 1e-8).any():
        raise ValueError('Internal sampling gap exceeds 0.61 m')
    knots = np.r_[lower, depth[(depth > lower) & (depth < upper)], upper]
    return float(trapezoid(np.interp(knots, depth, theta), knots) * 1000)


def surveyed_elevation(description):
    """Read the original survey Z field embedded in the well point layer."""
    match = re.search(r'<td>Z</td>\s*<td>([-\d.]+)</td>', description.replace('\\n', '\n'))
    if match is None:
        raise ValueError('Well layer has no surveyed Z elevation')
    return float(match.group(1))


WELL_IDS = [1, 2, 3, 5, 6, 7, 10, 12, 13, 14, 15, 16]


def well_colors():
    import seaborn as sns
    return dict(zip(WELL_IDS, sns.color_palette('tab20', len(WELL_IDS))))
