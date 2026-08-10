# The planar–spherical Green's function relation (Assignment 2, Q2)

## Two different relations

These are easy to conflate, and only the second one is what Q2 asks for.

### 1. The operator identity — relates *solutions of the source-free equation*

$$\nabla^2 \phi = \frac{1}{r}\frac{\partial^2 (r\phi)}{\partial r^2}$$

So if $\phi(r)$ solves the homogeneous equation in spherical geometry, then $u = r\phi$ solves the planar one. This is what makes the sphere and the slab share flux shapes: $\dfrac{\sin(Br)}{r}$ and $\cos(Bx)$ in the criticality problem of Assignment 1, Q3.

It says nothing about how a *source* maps between geometries.

### 2. The source superposition — relates *Green's functions*

This is the one needed here. An infinite plane source is a continuum of point sources spread over the plane, so the planar Green's function is the integral of the point (spherical) one over that plane.

## Deriving it

Put the source plane at $x = 0$ with unit strength per unit area, and take a field point at perpendicular distance $x$. A source element at radius $\rho$ from the foot of the perpendicular sits a distance

$$r = \sqrt{x^2 + \rho^2}$$

away, and the area element is $dA = 2\pi\rho \, d\rho$. Superposing,

$$\phi_{\text{pl}}(x,t) = \int_0^\infty \phi_{\text{pt}}\!\left(\sqrt{x^2+\rho^2},\, t\right) 2\pi\rho \, d\rho$$

Substituting $r^2 = x^2 + \rho^2$, hence $r\,dr = \rho\,d\rho$, with $\rho = 0 \mapsto r = |x|$ and $\rho \to \infty \mapsto r \to \infty$:

$$\boxed{\;\phi_{\text{pl}}(x,t) = 2\pi \int_{|x|}^{\infty} r\, \phi_{\text{pt}}(r,t)\, dr\;}$$

Differentiating with respect to $x > 0$ inverts it:

$$\frac{\partial \phi_{\text{pl}}}{\partial x} = -2\pi x\, \phi_{\text{pt}}(x,t)
\qquad\Longrightarrow\qquad
\phi_{\text{pt}}(r,t) = -\frac{1}{2\pi r} \left.\frac{\partial \phi_{\text{pl}}}{\partial x}\right|_{x=r}$$

Two things to note:

- **Time is a spectator.** This is a purely spatial superposition, applied independently at each $t$.
- **It is a relation between scalar fluxes.** The angular flux does not transform this simply, because each point source on the plane contributes along a *different* direction at the field point.

## Why "just multiply by $r$" feels right

Check the steady-state diffusion Green's functions:

$$\phi_{\text{pt}}(r) = \frac{e^{-\kappa r}}{4\pi D r},
\qquad
\phi_{\text{pl}}(x) = \frac{e^{-\kappa |x|}}{2 D \kappa}$$

Then $r\,\phi_{\text{pt}} = \dfrac{e^{-\kappa r}}{4\pi D}$, which *is* proportional to $\phi_{\text{pl}}$ — so the guess appears to work.

But it is an accident of the pure exponential: integrating $e^{-\kappa r}$ simply reproduces $e^{-\kappa r}$. Confirming with the actual relation,

$$2\pi \int_{|x|}^{\infty} r \cdot \frac{e^{-\kappa r}}{4\pi D r}\, dr
= \frac{1}{2D}\int_{|x|}^{\infty} e^{-\kappa r} dr
= \frac{e^{-\kappa |x|}}{2 D \kappa} \quad\checkmark$$

The Paasschens shapes are nothing like exponentials in $r$, so the proportionality breaks and the integral has to be done properly.

## A check worth doing before the messy part

Apply the relation to the uncollided term of Paasschens, with $v = \Sigma_t = 1$:

$$2\pi \int_{|x|}^{\infty} r \cdot \frac{e^{-t}}{4\pi r^2}\,\delta(r - t)\, dr
= \frac{e^{-t}}{2t}\,\Theta(t - |x|)$$

Now derive the same thing independently. Pure streaming from a plane isotropic pulse gives the uncollided angular flux $\psi = \tfrac{1}{2}e^{-t}\delta(x - \mu t)$, so

$$\phi_{\text{unc}}(x,t) = \int_{-1}^{1} \tfrac{1}{2} e^{-t} \delta(x - \mu t)\, d\mu
= \frac{e^{-t}}{2t}, \qquad |x| < t$$

They agree, which validates both the relation and the normalization before the collided term is touched.

## Applying it to the collided term

$$\phi_{\text{pl}}(x,t) = \underbrace{\frac{e^{-t}}{2t}\Theta(t-|x|)}_{\text{analytic}}
\;+\; \underbrace{2\pi \int_{|x|}^{t} r\, \phi_{\text{coll}}(r,t)\, dr}_{\text{numerical}}$$

The upper limit is $t$, not $\infty$, because $\Theta(vt - r)$ cuts the collided term off at the causal front.

**Watch the endpoint.** Writing $u = 1 - r^2/t^2$, the collided term behaves as

$$\phi_{\text{coll}} \propto u^{1/8}\, G\!\left(t\, u^{3/4}\right),
\qquad G(\xi) \sim \xi^{-1/2} \ \ (\xi \to 0)$$

so $G \sim t^{-1/2} u^{-3/8}$ and the integrand goes as $u^{1/8 - 3/8} = u^{-1/4}$. Since $u \approx 2(t-r)/t$ near the front, this is an integrable $(t-r)^{-1/4}$ singularity — finite, but fatal to naive Gauss–Legendre or Simpson. Substitute it away (e.g. $r = t\sqrt{1-w^{4}}$) or use a quadrature built for endpoint singularities, such as tanh–sinh.

## Then map to general $c$

Question 1's scaling relation carries the $c = 1$ result to every other case:

$$\phi(x,t;c) = c\, e^{-(1-c)t}\, \phi(cx,\, ct;\, 1)$$

so the planar Green's function only ever has to be constructed once.
