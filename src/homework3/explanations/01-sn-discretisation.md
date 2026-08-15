# 01 — The S_N Sweep in a Slab and in a Sphere

**One cell-balance equation serves both geometries; the sphere only adds a second outgoing
face, in angle, whose coefficient is fixed by requiring that a flat flux solve the discrete
equation exactly.**

## The slab

With Gauss–Legendre ordinates `(mu_m, w_m)` on `[-1, 1]` summing to 2, the discrete equation is

    mu_m dpsi_m/dx + Sigma_t psi_m = S/2,     phi = sum_m w_m psi_m

Integrating over cell `i` and closing with the diamond relation
`psi_i = (psi_out + psi_in)/2` gives the sweep

    psi_i = [ |mu_m| (A_out + A_in) psi_in + (S_i/2) V_i ] / [ Sigma_t V_i + 2 |mu_m| A_out ]
    psi_out = 2 psi_i - psi_in

with `A = 1` and `V = dx` in the slab. Reflection at `x = 0` and vacuum at `x = a/2` mean the
`mu < 0` ordinates are swept first, from the outer face inwards; each one leaves behind the
flux that its mirror ordinate `-mu` takes as its incoming value. That is why the ordinates are
kept in ascending `mu`: every inward ordinate is then swept before the outward one it feeds.

## The sphere

The spherical equation carries an angular derivative,

    mu dpsi/dr + (1 - mu^2)/r dpsi/dmu + Sigma_t psi = S/2

Its conservative form integrates over the cell to

    mu_m (A_p psi_p - A_m psi_m)
      + (A_p - A_m)/w_m [ alpha_{m+1/2} psi_{m+1/2} - alpha_{m-1/2} psi_{m-1/2} ]
      + Sigma_t V_i psi_i = (S_i/2) V_i

with `A = 4 pi r^2` and `V = (4 pi/3)(r_p^3 - r_m^3)`. This is the slab balance with a second
outgoing face — the angular one — so `cell_flux` in `sn.py` handles both by taking a list of
`(out coefficient, in coefficient, incoming flux)` links and closing each with the same
diamond relation, in space and in angle alike.

## Why `alpha_{m+1/2} = alpha_{m-1/2} - w_m mu_m`

Put `psi = 1` everywhere and `S/2 = Sigma_t`. The removal and source terms cancel, the
streaming term leaves `mu_m (A_p - A_m)`, and the angular term leaves
`(A_p - A_m)(alpha_{m+1/2} - alpha_{m-1/2})/w_m`. Requiring the two to cancel gives the
recursion, and `alpha_{1/2} = 0` starts it. Because `sum_m w_m mu_m = 0`, it also ends at
`alpha_{N+1/2} = 0`, so no current leaks out at `mu = +1`.

The exact coefficient is `alpha_{m+1/2} = 1 - mu_{m+1/2}^2`, which is what the angular
integral actually produces; the recursion is used instead because the quadrature
`integral f dmu ~ w_m f(mu_m)` is not exact, and consistency with a flat flux matters more
than consistency with the integral. **A factor of two here is not cosmetic:** the pair
`{coefficient (A_p - A_m)/w_m, recursion -w_m mu_m}` and the pair
`{coefficient (A_p - A_m)/(2 w_m), recursion -2 w_m mu_m}` are the same scheme, and mixing
them moved the `c = 1.5` critical radius from 1.686 to 1.607 mean free paths — a 5 % error
that no mesh refinement removes.

## The starting direction

The recursion needs `psi_{i, 1/2}`, the half-angle flux at `mu_{1/2} = -1`. There
`1 - mu^2 = 0`, the angular derivative drops out, and what is left,

    -dpsi/dr + Sigma_t psi = S/2

is the plain slab sweep at `|mu| = 1`, run inwards from the vacuum boundary. `_starting_direction`
solves exactly that, on the same mesh and with the same diamond closure.

## The centre of the sphere

`A_{1/2} = 0`, so the innermost face carries no current and the balance does not determine the
flux there — the diamond closure does, `psi_{1/2} = 2 psi_1 - psi_{3/2}`, and that value is
what the reflection `psi(0, mu) = psi(0, -mu)` hands to the outward sweep. The uncollided test
in `04` confirms the treatment: a pure absorber with a uniform source reproduces the analytic
centre flux `(1 - e^{-Sigma R})/Sigma` to five digits.
