"""Self-checks of the reference implementation (NumPy only):  python test_reference.py
 1. EXCL uploads at most one unit per block; Madow returns floor or ceil of the recording mass per recording.
 2. The marginal upload frequency of every unit equals its probability (within Monte-Carlo error) for the three draws.
 3. V1 of the block totals equals the estimate V-hat of the paper; the shore estimate is unbiased for a fixed predictor.
"""
import numpy as np
import samplers as S
import estimator as E
from demo_synthetic import make_population

rng = np.random.default_rng(5)
s_dev, y_dev, _ = make_population(4000, rng)
s, y, gaps = make_population(3000, rng)
cal = S.calibrate(s_dev, y_dev); jb = S.bin_of(cal, s); b = 0.1
forced = S.block_flags(cal['q'][jb], float(y_dev.mean()), 2, gaps)
p = S.external_probs('AM', s, cal, b, s_dev)
R = 4000; n = len(p)
excl = S.simulate_excl_openloop(p, forced, R, np.random.default_rng(1), 200)
madow = S.simulate_madow(p, gaps, forced, R, np.random.default_rng(2), 200)
indep, _ = S.draw_independent(p, R, np.random.default_rng(3))
starts = np.flatnonzero(excl['NEW'][0])
assert np.add.reduceat(excl['S'].astype(int), starts, axis=1).max() <= 1
seg = S.segments(gaps, n); cnt = np.add.reduceat(madow['S'].astype(int), np.flatnonzero(seg), axis=1)
mass = madow['segmass']
assert (cnt >= np.floor(mass - 1e-9)).all() and (cnt <= np.ceil(mass + 1e-9)).all()
for name, Sx in (('independent', indep), ('EXCL', excl['S']), ('Madow', madow['S'])):
    z = (Sx.mean(0) - p) / np.sqrt(p * (1 - p) / R)
    print('%-12s max |z| of the unit-level upload frequencies: %.2f' % (name, np.abs(z).max()))
    assert np.abs(z).max() < 5.0
f = np.full(n, y.mean())                                          # fixed predictor
th = y.mean()
for name, Sx in (('EXCL', excl['S']), ('Madow', madow['S'])):
    est = f.mean() + (Sx * (y - f) / p).sum(1) / n
    se = est.std(ddof=1) / np.sqrt(R)
    print('%-12s bias of the estimate: %.5f (standard error %.5f)' % (name, est.mean() - th, se)); assert abs(est.mean() - th) < 4.5 * se
rec_unit = np.cumsum(np.where(np.arange(n) == 0, True, gaps)) - 1
xr = excl['S'] * (y - f)[None, :] / p[None, :]
v1, v2, nb = E.block_v1v2(xr, excl['NEW'], rec_unit, n)
vhat = (excl['S'] * (y - f)[None, :] ** 2 / p[None, :] ** 2).sum(1) / n ** 2
assert np.allclose(v1, vhat, rtol=1e-10, atol=0)
est_e = f.mean() + (excl['S'] * (y - f) / p).sum(1) / n
print('V1 equals V-hat; mean V1 / empirical variance = %.3f, mean V2 / empirical variance = %.3f' % (v1.mean() / est_e.var(ddof=1), v2.mean() / est_e.var(ddof=1)))
print('all checks passed')
