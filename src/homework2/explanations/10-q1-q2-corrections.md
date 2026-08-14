# 10 — Corrections to Q1 and Q2

**Q1's scaling identity is four lines once three slips are fixed, and the factor `c` does not
come from the PDE at all; Q2 is not stuck at "evaluate numerically" — it has a closed form in
the Dawson function, and the `Sigma_t`'s cancel.**

Short version:

- **Q1** is provable in four lines, but three separate slips block it: a missing $\Sigma_t$ on the scattering term, a chain rule applied in the wrong direction, and then a *second* chain rule applied to a function that had already been differentiated. There is also a missing ingredient: the prefactor $c$ does not come from the PDE at all.
- **Q2** is *not* stuck at "evaluate numerically." The integral collapses in two substitutions to something elementary, and with Paasschens' $G$ it has a **closed form** in terms of the Dawson/`erfi` function. The $\Sigma_t$'s all cancel. There are also three small bookkeeping errors (a dropped $2\pi$, a `sign` that should be a $\Theta$, and a mislabeled left-hand side).

---

## Question 1

### 1.1 The two equations you are relating

This has to be pinned down first, because the identity is a statement *about a specific pair of equations*. With scattering ratio $c$ (isotropic, one speed), the transport equation is

$$\frac{1}{v}\frac{\partial \psi}{\partial t} + \mu \frac{\partial \psi}{\partial x} + \Sigma_t \psi
= \frac{c\,\Sigma_t}{2}\int_{-1}^{1}\psi\,d\mu' \tag{Q1.1}$$

and the $c=1$ equation, written in its own variables $(X,\mu,T)$, is

$$\frac{1}{v}\frac{\partial \Psi}{\partial T} + \mu \frac{\partial \Psi}{\partial X} + \Sigma_t \Psi
= \frac{\Sigma_t}{2}\int_{-1}^{1}\Psi\,d\mu' \tag{Q1.2}$$

> **Error 1.** Your Eq. (2) has $\tfrac12\int\tilde\psi\,d\mu'$ on the right, with no $\Sigma_t$. The scattering source is $c\Sigma_s^{\rm tot}$-weighted; for $c = 1$ it is $\tfrac{\Sigma_t}{2}\int\Psi\,d\mu'$. Dropping $\Sigma_t$ is what later leaves you with a leftover term $\tfrac{c-1}{v}\tilde\psi$ that cannot cancel against anything — note it has units of $\psi/(\text{length})\cdot(\text{time}/\text{length})$, while every other term carries $\Sigma_t\psi$. The mismatch is the tell.

### 1.2 The chain rule runs the other way

Define

$$\tilde\psi(x,\mu,t) \;\equiv\; \Psi(cx,\,\mu,\,ct), \qquad X = cx,\quad T = ct.$$

Then

$$\frac{\partial \tilde\psi}{\partial t} = c\,\Psi_T\Big|_{(cx,ct)}
\qquad\Longrightarrow\qquad
\Psi_T = \frac{1}{c}\frac{\partial \tilde\psi}{\partial t},$$

and likewise $\Psi_X = \tfrac{1}{c}\,\partial_x \tilde\psi$. Substituting into (Q1.2) and multiplying through by $c$:

$$\boxed{\;\frac{1}{v}\frac{\partial \tilde\psi}{\partial t} + \mu\frac{\partial \tilde\psi}{\partial x} + c\,\Sigma_t\,\tilde\psi
= \frac{c\,\Sigma_t}{2}\int_{-1}^{1}\tilde\psi\,d\mu' \;} \tag{Q1.3}$$

> **Error 2.** Your Eq. (3) reads $\frac{c}{v}\partial_t\tilde\psi + c\mu\,\partial_x\tilde\psi + \Sigma_t\tilde\psi = \tfrac12\int\tilde\psi$. The $c$'s are on the wrong terms. That equation is the one satisfied by $\Psi(x/c,\,t/c)$, i.e. the *inverse* scaling. The correct version, (Q1.3), moves the $c$ onto the **removal** term instead — and that is the whole mechanism of the identity, so it matters.

Read (Q1.3) against (Q1.1): they are identical *except* that $\tilde\psi$ has removal $c\Sigma_t$ where $\psi$ has removal $\Sigma_t$. The two right-hand sides already match. So the only job left for the exponential factor is to make up the removal deficit $\Sigma_t - c\Sigma_t = (1-c)\Sigma_t$. That immediately tells you what the exponent must be.

### 1.3 Units of the exponent

$(1-c)t$ is not dimensionless. The factor is

$$e^{-(1-c)\,v\Sigma_t t},$$

which reduces to $e^{-(1-c)t}$ in Paasschens' units $v = \Sigma_t = 1$ (equivalently, $t$ measured in mean free times). Since the rest of the assignment keeps $v$ and $\Sigma_t$ explicit, keep them here too.

### 1.4 The proof

Let $a \equiv v\Sigma_t$ and define the candidate

$$\Phi(x,\mu,t) \;=\; c\,e^{-(1-c)a t}\,\tilde\psi(x,\mu,t).$$

Insert it into the left-hand side of (Q1.1). Note $\tilde\psi$ is *already* a function of $(x,\mu,t)$, so $\partial_t$ and $\partial_x$ act on it directly — no second chain rule:

$$\frac{1}{v}\partial_t \Phi = c\,e^{-(1-c)at}\left[-(1-c)\Sigma_t\,\tilde\psi + \frac{1}{v}\partial_t\tilde\psi\right],
\qquad
\mu\,\partial_x \Phi = c\,e^{-(1-c)at}\,\mu\,\partial_x\tilde\psi,$$

$$\Sigma_t \Phi = c\,e^{-(1-c)at}\,\Sigma_t\,\tilde\psi.$$

Adding, the $\Sigma_t$ terms combine as $\Sigma_t - (1-c)\Sigma_t = c\Sigma_t$:

$$\text{LHS} = c\,e^{-(1-c)at}\left[\frac{1}{v}\partial_t\tilde\psi + \mu\,\partial_x\tilde\psi + c\Sigma_t\tilde\psi\right]
\overset{(Q1.3)}{=} c\,e^{-(1-c)at}\,\frac{c\Sigma_t}{2}\int_{-1}^{1}\tilde\psi\,d\mu'
= \frac{c\Sigma_t}{2}\int_{-1}^{1}\Phi\,d\mu',$$

which is exactly (Q1.1). $\blacksquare$

> **Error 3.** In your Eq. (38) the time-derivative term came out as $\frac{c^2}{v}e^{-(1-c)t}\partial_t\tilde\psi$ and the streaming term as $\mu c^2 e^{-(1-c)t}\partial_x\tilde\psi$. The extra factor of $c$ is a chain rule applied twice: once when you defined $\tilde\psi$, and again here. Once $\tilde\psi(x,\mu,t)$ is defined as a function of the unscaled variables, $\partial_t \tilde\psi$ is just $\partial_t\tilde\psi$.

### 1.5 The missing half of the proof — where the factor $c$ comes from

The PDE is linear and homogeneous, so if $\Phi$ solves it then so does $\lambda\Phi$ for any constant $\lambda$. **The differential equation alone can never produce the prefactor $c$.** It is fixed by the source/initial condition, and a complete proof has to say so.

For the plane-pulse Green's function the initial condition is $\psi(x,\mu,0^+) = \frac{v}{2}\delta(x)$ for every $c$. Evaluating the right-hand side of the identity at $t\to 0^+$, and using $\delta(cx) = \delta(x)/|c|$ for $c>0$:

$$c\,e^{0}\,\psi(cx,\mu,0^+;1) = c\cdot\frac{v}{2}\delta(cx) = c\cdot\frac{v}{2}\cdot\frac{\delta(x)}{c} = \frac{v}{2}\delta(x). \;\checkmark$$

So the prefactor $c$ is precisely the Jacobian of the spatial rescaling $x \mapsto cx$ acting on the source delta. This is worth one sentence in the writeup — it is the part that makes the identity a *theorem* rather than a statement about the PDE up to an unknown constant.

---

## Question 2

### 2.1 Bookkeeping errors first

Using $a \equiv v\Sigma_t t$ throughout.

1. **Left-hand side is mislabeled.** Lines 49–51 all read $\phi_{\rm pt}(r,t) = \dots$, but the object being computed is $\phi_{\rm pl}(x,t)$. The right-hand side has no $r$ in it after integration, so as written the equation is inconsistent.

2. **A dropped $2\pi$.** In line 49 the first (uncollided) integral carries $2\pi$, the second (collided) one does not. The $2\pi$ belongs to the relation $\phi_{\rm pl}(x) = 2\pi\int_{|x|}^{\infty} \phi_{\rm pt}(r)\,r\,dr$ itself, so it multiplies *both* terms. It is missing from lines 50 and 51 as well.

3. **`sign` should be $\Theta$.** You write $\text{sign}(vt-|x|)$ in front of the collided term. When $|x| > vt$ the factor $\Theta(vt-r)$ in the integrand kills the integrand for *every* $r \ge |x| > vt$, so the term is **zero**, not negative. The $\Theta(vt-r)$ converts an integral over $[|x|,\infty)$ into one over $[|x|, vt]$ *only when that interval is non-empty*; the correct bookkeeping is

$$\int_{|x|}^{\infty}\!\!(\cdots)\Theta(vt-r)\,dr \;=\; \Theta(vt-|x|)\int_{|x|}^{vt}\!\!(\cdots)\,dr.$$

4. **Density vs. flux.** Paasschens' $P(\mathbf r,t)$ is a *normalized probability density*: his Eq. (34) gives $\int d\mathbf r\, P_N = \frac{1}{N!}(ct/l)^N e^{-ct/l}$, which sums over $N$ to exactly 1. So $P$ is a number density, and the scalar flux is $\phi = vP$. Decide which one you are reporting and say so; if you want flux, there is a factor $v$ multiplying everything below.

Everything else in the setup is right. In particular your restoration of units in the Paasschens formula — $x_{\rm paper} = \Sigma_t r$, $t_{\rm paper} = v\Sigma_t t$, $(4\pi l ct/3)^{3/2} \to (4\pi v t/(3\Sigma_t))^{3/2}$, and the overall $\Sigma_t^3$ from the density normalization — is correct, and so is the $\xi = r/vt$ substitution in line 51 including the $(vt)^{1/2}$ prefactor. The plane/point superposition relation you use is also right.

### 2.2 The integral is not hard — it is elementary

This is the main point. Your line 51 leaves

$$I \;=\; \int_{|x|/vt}^{1}(1-\xi^2)^{1/8}\,G\!\left(a(1-\xi^2)^{3/4}\right)\xi\,d\xi.$$

Do **not** integrate this in $\xi$. Two substitutions kill it.

**Step 1: $u = 1-\xi^2$**, so $\xi\,d\xi = -\tfrac12 du$, and the limits flip to $u: 0 \to u_0 \equiv 1 - \dfrac{x^2}{(vt)^2}$:

$$I = \frac{1}{2}\int_0^{u_0} u^{1/8}\,G\!\left(a\,u^{3/4}\right)du.$$

**Step 2: $w = a\,u^{3/4}$**, so $u = (w/a)^{4/3}$, $du = \tfrac43 a^{-4/3} w^{1/3}dw$, $u^{1/8} = a^{-1/6}w^{1/6}$. The exponents add as $\tfrac16 + \tfrac13 = \tfrac12$:

$$\boxed{\;I \;=\; \frac{2}{3\,a^{3/2}}\int_0^{w_0}\sqrt{w}\;G(w)\,dw\;},
\qquad w_0 \equiv a\left(1-\frac{x^2}{(vt)^2}\right)^{3/4}.$$

The awkward fractional powers $\tfrac18$ and $\tfrac34$ were engineered by Paasschens to interpolate between $d=2$ and $d=4$; under $w = au^{3/4}$ they conspire into a plain $\sqrt{w}$. **This also completely removes the endpoint singularity** flagged in [09](09-planar-spherical-relation.md) — no tanh–sinh quadrature, no $r = t\sqrt{1-w^4}$ trick needed. The note's warning about the $(t-r)^{-1/4}$ behavior is correct as far as it goes, but it is an artifact of the variable, not of the integral.

### 2.3 Assembling it — the $\Sigma_t$'s cancel

Putting $I$ back through the prefactor of line 51 (with the missing $2\pi$ restored):

$$2\pi\,\frac{e^{-a}(vt)^{1/2}}{\left(4\pi/(3\Sigma_t)\right)^{3/2}}\cdot\frac{2}{3a^{3/2}}
= \frac{e^{-a}}{2vt}\sqrt{\frac{3}{\pi}},$$

because $a^{3/2} = (v\Sigma_t t)^{3/2}$ cancels the $\Sigma_t^{3/2}$ entirely and $(vt)^{1/2}/(vt)^{3/2} = 1/(vt)$. So both terms share the same prefactor and the answer is

$$\boxed{\;n_{\rm pl}(x,t) \;=\; \frac{e^{-v\Sigma_t t}}{2vt}\left[\,1 \;+\; \sqrt{\frac{3}{\pi}}\int_0^{w_0}\sqrt{w}\,G(w)\,dw\,\right]\Theta(vt-|x|)\;}$$

with $a = v\Sigma_t t$, $\;w_0 = a\left(1 - x^2/(vt)^2\right)^{3/4}$, and $\phi_{\rm pl} = v\,n_{\rm pl}$.

The uncollided term reduces to the "1", which is a good structural check: the collided part is measured in units of the uncollided plateau.

### 2.4 Closed form

Paasschens' interpolation (his Eq. 36b) is

$$G(x) \;\approx\; e^{x}\sqrt{1 + \frac{b}{x}}, \qquad b = 2.026,$$

and the miracle is that the $\sqrt{w}$ from Step 2 is exactly what rationalizes it:

$$\sqrt{w}\,G(w) = \sqrt{w}\;e^{w}\sqrt{1+\tfrac{b}{w}} = e^{w}\sqrt{w+b}.$$

So the remaining integral is elementary. With $s = w+b$, $\int e^s\sqrt{s}\,ds = e^s\sqrt s - \tfrac12\int e^s s^{-1/2}ds$ and $\int e^s s^{-1/2}ds = \sqrt\pi\,\mathrm{erfi}(\sqrt s)$:

$$\int_0^{w_0} e^{w}\sqrt{w+b}\;dw
= e^{w_0}\sqrt{w_0+b} - \sqrt{b} - \frac{\sqrt\pi}{2}e^{-b}\Big[\mathrm{erfi}\big(\sqrt{w_0+b}\big) - \mathrm{erfi}\big(\sqrt{b}\big)\Big].$$

**Use the Dawson form instead for numerics.** With $\mathrm{erfi}(z) = \tfrac{2}{\sqrt\pi}e^{z^2}D(z)$, where $D$ is the Dawson function (`scipy.special.dawsn`), the $e^{w_0}$ factors out cleanly:

$$\int_0^{w_0} e^{w}\sqrt{w+b}\;dw
= e^{w_0}\Big[\sqrt{w_0+b} - D\big(\sqrt{w_0+b}\big)\Big] - \Big[\sqrt{b} - D\big(\sqrt{b}\big)\Big].$$

Substituting into the boxed result, and noting $e^{-a}e^{w_0} = e^{-(a-w_0)} \le 1$ always:

$$n_{\rm pl}(x,t) = \frac{\Theta(vt-|x|)}{2vt}\left\{ e^{-a} + \sqrt{\tfrac{3}{\pi}}\left( e^{-(a-w_0)}\Big[\sqrt{w_0+b}-D(\sqrt{w_0+b})\Big] - e^{-a}\Big[\sqrt{b}-D(\sqrt{b})\Big]\right)\right\}$$

Every factor here is bounded, so this evaluates without overflow out to $t \sim 10^2$ mean free times and beyond. **No quadrature anywhere.**

### 2.5 If you want the exact $G$ rather than the interpolation

The same $\sqrt{w}$ reduction also makes the *exact* series integrable term by term. Paasschens' Eq. (36b) is

$$G(x) = 8(3x)^{-3/2}\sum_{N=1}^{\infty}\frac{\Gamma\!\left(\tfrac34 N + \tfrac32\right)}{\Gamma\!\left(\tfrac34 N\right)}\frac{x^N}{N!},$$

so $\sqrt{w}\,G(w) \propto \sum_N c_N w^{N-1}$ and

$$\int_0^{w_0}\sqrt{w}\,G(w)\,dw = \frac{8}{3^{3/2}}\sum_{N=1}^{\infty}\frac{\Gamma\!\left(\tfrac34 N+\tfrac32\right)}{\Gamma\!\left(\tfrac34 N\right)}\frac{w_0^{\,N}}{N\cdot N!}.$$

(Sum in log-space via `gammaln` to avoid overflow. Note the argument is $\tfrac34 N$, not $\tfrac14 N$ — easy to misread in the scanned PDF; it is forced by the exponent $\tfrac34 N - 1$ in his Eq. 35.)

### 2.6 Verification

All checks run against the formulas above.

| check | result |
|---|---|
| Dawson closed form vs. brute-force `quad` on the original $r$-integral | agrees to $\sim10^{-13}$ relative, for $t \in [0.5, 200]$ and $\|x\|/vt \in [0, 0.999]$ |
| $\int_{-vt}^{vt} n_{\rm pl}\,dx$ with the **exact series** $G$ | $= 1.00000000$ at $t = 0.3, 1, 3, 10, 30$ — exact particle conservation for $c=1$ |
| same, with the interpolation $G$ | $1.0008$–$1.0147$, i.e. the $\lesssim 2\%$ error Paasschens quotes for Eq. (36) |
| behavior at the front $\|x\|\to vt$ | $w_0 \to 0$ and the bracket $\to \sqrt b - \sqrt b = 0$, so $n_{\rm pl}\to e^{-a}/(2vt)$, the uncollided value alone $\checkmark$ |
| diffusion limit, $n_{\rm pl}(0,t)$ vs. $(4\pi Dt)^{-1/2}$, $D = v/(3\Sigma_t)$ | $t=10$: 0.1621 vs 0.1545; $t=100$: 0.04911 vs 0.04886; $t=300$: 0.028258 vs 0.028209 $\checkmark$ |

The exact-normalization result is the strong one: it pins down the overall prefactor $\sqrt{3/\pi}\,/(2vt)$ and the cancellation of $\Sigma_t$ independently of any numerics.

### 2.7 So: is this "the best form you can get to"?

Yes for the collided integral, and it is better than you thought — you were one substitution away from a closed form, not stuck at a numerical evaluation. What *cannot* be improved:

- $G$ itself is an interpolation between the exactly-invertible $d=2$ and $d=4$ results; the true $d=3$ Green's function has no closed form (Paasschens inverts it numerically in his Appendix). So "closed form" here means closed form *given* Eq. (36), which is accurate to a few percent.
- The uncollided term is already exact.

And your closing remark is right: with Q1's scaling relation, $\;n(x,t;c) = c\,e^{-(1-c)v\Sigma_t t}\,n(cx,ct;1)$, the $c=1$ result above is all you ever need to build.
