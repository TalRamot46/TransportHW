# 02 — Evaluating the Closed Form

**The Q2 integral has a closed form, and almost every line of `exact.py` exists to stop that
closed form overflowing.**

Report §2 derives the integral and its `erfi` form. This is why the code does not use it.

## Dawson instead of `erfi`

The report's antiderivative is

    int_0^w0 e^w sqrt(w+b) dw = e^{w0} sqrt(w0+b) - sqrt(b) - (sqrt(pi)/2) e^{-b} [erfi(sqrt(w0+b)) - erfi(sqrt(b))]

Every term grows like `e^{w0}`, and the physical flux carries a compensating `e^{-t}`. Formed
naively, `erfi` **overflows past `w0 ~ 20`** — and `t = 15` is in the figure set, so this is
not a hypothetical.

`_dawson_antiderivative` uses Dawson's function instead, `D(z) = (sqrt(pi)/2) e^{-z^2} erfi(z)`,
which is bounded (it peaks near 0.54 and decays). The exponential is cancelled *analytically*
before evaluation rather than numerically afterwards, so `collided_integral` returns
`e^{-w0} int_0^w0 ...` directly and never forms a large intermediate:

    return _dawson_antiderivative(w0 + B_INTERP) - exp(-w0) * _dawson_antiderivative(B_INTERP)

That is the same identity as the report's, rearranged so that nothing in it is ever larger than
about `sqrt(w0 + b)`.

## The series is summed in log space

`G(w, form="series")` and the series branch of `collided_integral` build their terms as

    exp(gammaln(...) - gammaln(...) - gammaln(...) + N * log(w))

rather than forming ratios of factorials. `Gamma(3N/4 + 3/2)/Gamma(3N/4)/N!` at `N = 200`
is a ratio of numbers far outside double range; in log space it is a subtraction of three
`gammaln` calls. `SERIES_TERMS = 200` is converged to machine precision for `w0 <= 30`.

**`_safe_log` exists for one reason:** `w = 0` at the causal front, and `log(0)` would emit a
divide warning even though the term is later discarded. It substitutes 1 inside the log and the
result is masked afterwards by the `np.where(w0 > 0.0, total, 0.0)`. The mask is what is
correct; the substitution only keeps the log quiet.

## The front is handled twice, deliberately

`phi_c1` guards the causal front `|x| < t` in two separate places:

```
w0 = front_distance(np.where(inside, x, 0.0), t)     # feed the integral a safe value
return np.where(inside, (np.exp(-t) + collided) / (2.0 * t), 0.0)   # mask the result
```

The first `np.where` is not redundant with the second. NumPy evaluates the whole array before
masking, so without it `collided_integral` would be called on points outside the front where
`front_distance` has already clamped to zero — harmless in value, but it is the clamp in
`front_distance` (`np.maximum(1 - (x/t)^2, 0)`) that prevents a fractional power of a negative
number. Removing any one of the three guards produces `nan`s that only appear at some `t`.

## Where the units went

`Sigma_t = v = 1` is baked in, not carried. Paasschens' argument `a = v Sigma_t t` reduces to
`t`, which is why `front_distance(x, t)` takes no cross section, and the report's `Sigma_t`
factors have cancelled analytically (report §2). Restoring general units means revisiting
`front_distance` and `phi_c1` together, not just multiplying at the end.
