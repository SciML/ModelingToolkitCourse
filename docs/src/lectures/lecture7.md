# Lecture 7: Pantelides Algorithm, Dummy Derivatives, and Demonstrations

### The Mathematical Structure of DAEs

There are both differential equations and algebraic equations in acasual models.
Thus, a generic acasual model is a system of differential-algebraic equations
(DAEs). In general, additional processing steps are required to simulate DAEs.
To see this more clearly, we will analyze the implicit Euler algorithm which is
the most basic form of both Runge-Kutta and linear multistep methods for DAEs.

Given the DAE of the form
```math
0 = F(u'(t), u(t), p, t),
```
where ``F: (\mathbb{R}^n, \mathbb{R}^n, \mathbb{R}^m, \mathbb{R}) \rightarrow
\mathbb{R}^n``. The implicit Euler solves for ``u(t+h)`` from
```math
0 = F\left(\frac{\hat{u}(t+h) - u(t)}{h},\; \hat{u}(t+h),\; p,\; t+h\right),
```
with fixed ``t, h, p``, and ``u(t)``.

Numerically, we use Newton's method to solve potentially nonlinear equations by
solving the best approximating linear equations (i.e. the Jacobian of the
nonlinear function with respect to the unknowns) iteratively to refine an
initial guess. By the chain rule, we have
```math
\frac{\partial F}{\partial u} = \frac{1}{h}F_{u'} + F_{u}.
```
The Newton iteration is then
```math
\begin{align}
\frac{\partial F}{\partial u} \Delta^{[i]} &= F\left(\frac{\hat{u}^{[i]}(t+h) - u(t)}{h}, \hat{u}^{[i]}(t+h), p, t+h\right) \\
\hat{u}^{[i+1]} &= \hat{u}^{[i]} - \Delta^{[i]},
\end{align}
```
where ``\{\cdot\}^{[i]}`` denotes the iteration variable at the ``i``-th
iteration. Thus, for the implicit Euler algorithm to work, ``\lambda F_{u'} +
F_{u}`` has to be non-singular for sufficiently small ``\lambda\in\mathbb{R}``.
In fact, this conclusion holds for all linear multistep methods and Runge-Kutta
methods, that is, they all need to solve an iteration matrix of the form
``\lambda F_{u'} + F_{u}``. Curious readers can check this
[development documentation](https://github.com/SciML/DiffEqDevMaterials/blob/master/newton/output/main.pdf)
on how the Newton iteration is set up for all the other cases.

The above criterion is too hard to verify at compile time, since we do not yet
know the exact values to solve the DAE. Formally, the above criterion is
```math
\forall u, u', \exists \lambda > 0, \det(\lambda F_{u'} + F_{u}) \ne 0.
```
To make it computationally easier to check, we can change the ``\forall`` to
``\exists``, which is
```math
\exists u, u', \lambda > 0, \det(\lambda F_{u'} + F_{u}) \ne 0,
```
so that we only need to validate a single instance. To make the criterion even
easier to check, we introduce the following definitions.

!!! definition "Sparse Matrix"

    A sparse matrix is a matrix that could contain structural zeros. Structural
    zeros are enties that are zero by construction denoted by ``\hat{0}``. We
    also define ``\mathbb{F}:=\mathbb{R}\cup\hat{0}``.

!!! definition "Incidence Matrix"

    A incidence matrix of the symbolically defined function ``f: \mathbb{R}^n
    \mapsto \mathbb{R}^m`` with respect to the indexed set of variables ``\{x_{j}
    \} `` is a sparse matrix ``M\in\mathbb{F}^{m \times n}`` defined by
    ```math
    M_{ij} := \begin{cases}
    1, & \text{if } x_{j} \text{ appears in the expression for computing the
    $i$-th output of $f$, i.e. } f(\{x_{j}\})_i \\
    \hat{0}, & \text{else}
    \end{cases}.
    ```

!!! note "Structural Zeros"

    For convince, we will simply use ``0`` and ``\mathbb{R}`` when structural
    zeros are obvious in context.

!!! definition "Structurally Non-singular"

    A matrix ``M`` is structurally non-singular if and only if there exists a
    set of elements for its structurally nonzero entries such that ``M``
