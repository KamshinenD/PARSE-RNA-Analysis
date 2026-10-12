"""Shared estimator for the selection-direction (regression-to-the-mean) controls.

One unit, one vote. Both cohorts are re-modelling comparisons in which a handful of
ribosomes carry a quarter of all matched pairs, so a pair-weighted mean is a mean
over ribosomes, not over molecules -- and the rest of the paper reports these
cohorts per molecule. The two directions are also PAIRED on the unit: a unit only
contributes if it has pairs selected in both directions, otherwise the comparison
silently swaps the sample between directions.

Reading the two directions together:

    forward = T + R      (selected in model A, measured in B)
    reverse = R - T      (selected in B, measured in A)

where R is the regression effect the selection guarantees and T the genuine,
direction-specific improvement, which works AGAINST the reverse direction. Hence
T = (forward - reverse) / 2, and the asymmetry is what the control actually tests.
"""
import numpy as np
import pandas as pd

N_BOOT, SEED = 4000, 42


def per_unit(df, unit, sel, frm, to, cut):
    sub = df[df[sel] < cut]
    return sub.groupby(unit).apply(lambda g: (g[to] - g[frm]).mean(), include_groups=False)


def paired_control(df, unit, col_a, col_b, cut, n_boot=N_BOOT, seed=SEED, r_floor=None):
    """Forward/reverse/asymmetry per unit, paired on units seen in both directions."""
    fwd = per_unit(df, unit, col_a, col_a, col_b, cut)
    rev = per_unit(df, unit, col_b, col_b, col_a, cut)
    both = sorted(set(fwd.index) & set(rev.index))
    # per-unit agreement between the two models: near-zero means the two records
    # are not really the same pairs, which mimics a pure-regression null
    def _r(g):
        if len(g) < 3 or g[col_a].nunique() < 2 or g[col_b].nunique() < 2:
            return np.nan          # undefined, not zero: too few points or no variance
        return g[col_a].corr(g[col_b])

    corr = df.groupby(unit).apply(_r, include_groups=False)
    keep = [u for u in both if r_floor is None or (pd.notna(corr.get(u)) and corr[u] >= r_floor)]
    F, R = fwd[keep].to_numpy(), rev[keep].to_numpy()
    d = F - R
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    bf, br, bd = F[idx].mean(1), R[idx].mean(1), d[idx].mean(1)
    ci = lambda a: [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]
    return {
        "n_units": int(len(keep)),
        "n_units_forward_only": int(len(set(fwd.index) - set(rev.index))),
        "n_units_reverse_only": int(len(set(rev.index) - set(fwd.index))),
        "n_pairs_forward": int((df[df[col_a] < cut][unit].isin(keep)).sum()),
        "n_pairs_reverse": int((df[df[col_b] < cut][unit].isin(keep)).sum()),
        "forward": {"mean": float(F.mean()), "ci95": ci(bf)},
        "reverse": {"mean": float(R.mean()), "ci95": ci(br)},
        "asymmetry": {"mean": float(d.mean()), "ci95": ci(bd),
                      "p_le_0": float((bd <= 0).mean()),
                      "units_with_forward_greater": int((d > 0).sum())},
        "implied_true_improvement": {"mean": float(d.mean() / 2),
                                     "ci95": [ci(bd)[0] / 2, ci(bd)[1] / 2],
                                     "definition": "(forward - reverse) / 2"},
        "median_unit_correlation": float(corr[keep].median()),
        "r_floor": r_floor,
    }


def permutation_null(df, unit, col_a, col_b, cut, n=200, seed=0):
    """Shuffle the partner score within each unit: destroys the pairing, keeps the
    marginals. Pure regression is direction-symmetric, so this is what 'no real
    correspondence' looks like."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        d = df.copy()
        d[col_b] = d.groupby(unit)[col_b].transform(lambda s: rng.permutation(s.to_numpy()))
        f = d[d[col_a] < cut]
        r = d[d[col_b] < cut]
        out.append(((f[col_b] - f[col_a]).mean(), (r[col_a] - r[col_b]).mean()))
    a = np.array(out)
    return {"forward_mean": float(a[:, 0].mean()), "reverse_mean": float(a[:, 1].mean()),
            "n_shuffles": n,
            "note": "pair-weighted; symmetric by construction when pairing is destroyed"}
