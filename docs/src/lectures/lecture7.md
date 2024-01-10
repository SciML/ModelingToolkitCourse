# Lecture 7: Numerical and Structural Characterizations for DAEs

## Numerical Integrability Criterion for DAEs

There are both differential equations and algebraic equations in acasual models.
Thus, a generic acasual model is a system of differential-algebraic equations
(DAEs). In general, additional processing steps are required to simulate DAEs.
To see this more clearly, we will analyze the implicit Euler algorithm which is
the most basic form of both Runge-Kutta and linear multistep methods for DAEs.

Given the DAE of the form
```math
\begin{equation}
0 = F(u'(t), u(t), p, t),
\end{equation}
```
where ``F: (\mathbb{R}^n, \mathbb{R}^n, \mathbb{R}^m, \mathbb{R}) \rightarrow
\mathbb{R}^n``. The implicit Euler solves for ``u(t+h)`` from
```math
\begin{equation}
0 = F\left(\frac{\hat{u}(t+h) - u(t)}{h},\; \hat{u}(t+h),\; p,\; t+h\right),
\end{equation}
```
with fixed ``t, h, p``, and ``u(t)``.

Numerically, we use Newton's method to solve potentially nonlinear equations by
solving the best approximating linear equations (i.e. the Jacobian of the
nonlinear function with respect to the unknowns) iteratively to refine an
initial guess. By the chain rule, we have
```math
\begin{equation}
\frac{\partial F}{\partial u} = \frac{1}{h}F_{u'} + F_{u}.
\end{equation}
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
Formally, the *numerical integrability criterion* is
```math
\begin{equation}
\forall u, u', \exists \lambda > 0, \det(\lambda F_{u'} + F_{u}) \ne 0.
\end{equation}
```

## Structural Analysis

The above criterion is too hard to verify at compile time, since we do not yet
know the exact values to solve the DAE. To make it computationally easier to
check, we can change the ``\forall`` to ``\exists``, which is
```math
\begin{equation}
\exists u, u', \lambda > 0, \det(\lambda F_{u'} + F_{u}) \ne 0,
\end{equation}
```
so that we only need to validate a single instance. To make the criterion even
easier to check, we introduce the following definitions.

!!! definition "Sparse Matrix"

    A sparse matrix is a matrix that could contain structural zeros. Structural
    zeros are entries that are zero by construction denoted by ``\hat{0}``. We
    call a entry structural nonzero if it is not a structural zero. We also
    define ``\mathbb{F}:=\mathbb{R}\cup\hat{0}``.

!!! definition "Incidence Matrix"

    A incidence matrix of the symbolically defined function ``f: \mathbb{R}^n
    \to \mathbb{R}^m`` with respect to the indexed set of variables ``\{x_{j}
    \} `` is a sparse matrix ``M\in\mathbb{F}^{m \times n}`` defined by
    ```math
    \begin{equation}
    M_{ij} := \begin{cases}
    1, & \text{if } x_{j} \text{ appears in the expression for computing the
    $i$-th output of $f$, i.e. } f(\{x_{j}\})_i \\
    \hat{0}, & \text{else}
    \end{cases}.
    \end{equation}
    ```
    We use ``\mathfrak{I}(f, \{x_j\}) := M`` to denote the incidence matrix of
    ``f`` with respect to ``\{x_j\}``.

!!! note "Structural Zeros"

    For convince, we will simply use ``0`` and ``\mathbb{R}`` when structural
    zeros are obvious in context.

!!! definition "Structurally Non-singular"

    A sparse matrix ``A`` is structurally non-singular if and only if there
    exists a set of real numbers when it replaces all the structural nonzero
    entries, the new matrix ``\mathfrak{N}(A)`` is numerically non-singular.

!!! example "Structurally Non-singular"
    - ``A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}`` is structurally
      non-singular but numerically singular because we can replace a nonzero such
      that we get ``\mathfrak{N}(A) = \begin{pmatrix} 1 & 10 \\ 1 & 1 \end{pmatrix}`` which is
      numerically non-singular.
    - ``A = \begin{pmatrix} 1 & 2 & 3\\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}`` is
      structurally non-singular but numerically singular because we can replace
      nonzeros such that we get
      ``\mathfrak{N}(A) = \begin{pmatrix} 10 & 2 & 3\\ 4 & 5 & 6 \\ 7 & 8 & 10 \end{pmatrix}``
      which is numerically non-singular.
    - ``\begin{pmatrix} 1 & 1 \\ \hat{0} & \hat{0} \end{pmatrix}`` is structurally
      singular because ``\begin{pmatrix} a & b \\ 0 & 0 \end{pmatrix}`` is
      numerically singular for all real ``a`` and ``b``.

Although structurally non-singular is a strictly weaker condition than
numerically non-singular, checking it seems to be more difficult on the surface
because we need to come up with an example that is numerically non-singular.
However, the following powerful theorem gives us a dramatically simpler way of
checking if a sparse matrix is structurally non-singular.

!!! info "Structurally Non-singularity Theorem"

    A square sparse matrix ``A`` is structurally non-singular if and only if
    there exists permutation matrices ``P`` and ``Q`` such that all the diagonal
    entries of ``PA`` and ``AQ`` are structural nonzeros.

    ##### Proof:
    - ``\Leftarrow``: Suppose ``A`` is a square sparse matrix such that all
      diagonal entries of ``PA`` and ``AQ`` are structural nonzeros, where ``P``
      and ``Q`` are permutation matrices. Let ``\mathfrak{N}(A) = P^{-1}I`` or
      ``\mathfrak{N}(A) = IQ^{-1}`` which are numerically non-singular. Thus,
      ``A`` is structurally non-singular.
    - ``\Rightarrow``: Suppose ``\hat{A}\in\mathbb{R}^{n\times n}`` is a square
      sparse matrix that is structurally non-singular. Let ``A = \mathfrak{N}
      (\hat{A})``, we have ``\det(A) \ne 0``. Expanding the definition of the
      determinant, we have
      ```math
      \begin{equation}
      \det(A)=\sum_{\sigma \in S_{n}}\operatorname{sgn}(\sigma)\prod_{i = 1}^n
      a_{i,\sigma (i)} \ne 0,
      \end{equation}
      ```
      where ``S_n`` denotes the set of all permutations of the set ``\{1, 2, ...,
      n\}`` (the symmetric group of order ``n``). For ``\det(A)\ne 0``, there
      must exist one ``\sigma\in S_n`` such that ``\prod_{i = 1}^n a_{i,\sigma
      (i)} \ne 0``. Note that we also have ``\prod_{i = 1}^n a_{i,\sigma
      (i)} = \prod_{i = 1}^n a_{\sigma^{-1}(i), i}``. Thus, the desired ``Q`` is
      the permutation matrix corresponding to ``\sigma``, and the desired ``P``
      is ``Q^{-1}``. ``\blacksquare``

The above condition is equivalent with checking the existence of a perfect
matching on the induced bipartite graph of the incidence matrix, and it can be
efficiently solved by using the augmenting path algorithm to find the maximum
cardinality matching.

!!! definition "Induced Bipartite Graph"

    Given an incidence matrix ``A\in\mathbb{R}^{m\times n}`` the induced
    bipartite graph ``(U, V, E)`` is a tuple of a set of source verties ``U =
    \{1, 2, ..., m\}``, a set of destination verties ``V\in \{1, 2, ..., n\}``,
    and edges between them ``E \subseteq U\times V`` defined by
    ```math
    \begin{equation}
    \forall (i, j) \in U\times V, (i, j) \in E \iff A_{i, j} = 1.
    \end{equation}
    ```
    Similarly, the induced bipartite graph of a sparse matrix is the induced
    bipartite graph of its induced incidence matrix.

!!! definition "Bipartite Matching"

    A matching of bipartite graph ``(U, V, E)`` is a set ``M \subseteq E`` where
    every vertex in ``U`` and ``V`` can appear at most once in ``M``. A matching
    ``M`` is perfect if ``|M| = |U| = |V|``. We call an edge in a matching
    matched, otherwise, free. It is often more convient to interprate matching
    as the function ``m: U \to (V \cup \emptyset)`` defined by
    ```math
    \begin{equation}
    x \mapsto \cup \{y | (x, y) \in M\}.
    \end{equation}
    ```

!!! info "Structural Non-singularity and Perfect Matching Equivalence Theorem"

    A sparse matrix ``A\in\mathbb{R}^{m\times n}`` is structurally non-singular
    if and only if its induced bipartite graph has a perfect matching.

!!! definition "Augmenting Path"

    Given a particular matching, an alternating path is a sequence of adjacent
    edges that alternate between being matched and free. In particular, an
    augmenting path is a alternating path that starts and ends with free edges.

!!! info "Augmenting Path Algorithm for Finding a Maximum Cardinality Matching"

    Input: bipartite graph ``G = (U, V, E)``.

    Output: matching ``M``.
    ```julia
    M = []
    for each u in U
        p ← find an augmenting path w.r.t. M that starts with u
        if p === nothing
            continue
        else
            add all free edges of p to M
            remove all matched edges of p from M
        end
    end
    ```
    Note that by the definition of augmenting paths, whenever ``p`` is not
    `nothing` in the above algorithm, we increase the cardinality of ``M`` by
    ``1``. We will assert without a proof that the above algorithm outputs a
    maximum cardinality matching, and in particular, if ``p`` is `nothing` for a
    source vertex ``u``, then no maximum cardinality matching contains an edge
    that starts with ``u``. More details of this algorithm including the search
    algorithm of an augmenting path are available in the original Pantelides
    paper [^Pantelides1988].

    Note that a perfect matching for ``G`` exists if and only if a maximum
    cardinality matching ``M`` satisifies ``|M| = |U| = |V|``.



[^Pantelides1988]: Pantelides, Constantinos C. "The consistent initialization of differential-algebraic systems." SIAM Journal on scientific and statistical computing 9.2 (1988): 213-231.

## Structural Integrability Criterion for DAEs

To utilize the structural analysis framework, we need weaken the integrability
criterion further from
```math
\begin{equation}
\exists u, u', \lambda > 0, \det(\lambda F_{u'} + F_{u}) \ne 0,
\end{equation}
```
to the *structural integrability criterion*
```math
\begin{equation}
\mathfrak{I}(F, \{u_i'\}) + \mathfrak{I}(F, \{u_i\}) \text{ is structurally
non-singular.}
\end{equation}
```

## Consistency Solvability Criterion for DAEs

Consider the DAE system
```math
\begin{equation}
F(\{x', y'\}, \{x, y\}, t) =
\begin{pmatrix}
f_1(x, t) \\
f_2(x', y', t)
\end{pmatrix} = 0,
\end{equation}
```
where ``f_1`` and ``f_2`` are some arbitrary smooth functions. We have
```math
\begin{equation}
\mathfrak{I}(F, \{x', y'\}) = \begin{pmatrix}
0 & 0 \\
1 & 1
\end{pmatrix},\quad
\mathfrak{I}(F, \{x, y\}) = \begin{pmatrix}
1 & 0 \\
0 & 0
\end{pmatrix}.
\end{equation}
```
Thus,
```math
\begin{equation}
\mathfrak{I}(F, \{u_i'\}) + \mathfrak{I}(F, \{u_i\}) = \begin{pmatrix}
1 & 0 \\
1 & 1
\end{pmatrix}
\end{equation}
```
is structurally non-singular, which means that the DAE system is structurally
integrable.

However, solving for a consistent initial condition ``u(t_0)`` and
``u'(t_0)`` is not as simple as simply solving for ``F(u'(t_0), u(t_0), t_0) =
0``. Because given a general DAE in the form of ``F(u', u, p, t) = 0``, all its
total time derivatives are also valid constraints, i.e.
```math
\begin{align}
F(u', u, p, t) &= 0 \\
F'(u'', u', u, p, t) &= 0  \\
F''(u''', u'', u', u, p, t) &= 0 \\
    &\vdots \nonumber
\end{align}
```
For the above example, we can differentiate the ``f_1`` equation once and get
```math
\begin{equation}
\begin{pmatrix}
f'_1(x', t) \\
f_2(x', y', t)
\end{pmatrix} = 0.
\end{equation}
```
Note that differentiating ``f'_1`` and ``f_2`` further is not necessary because
we will not get additional constraints for the states (``u(t)`` and ``u'(t)``)
of the DAE system.

We need a more systematic way of knowing when differentiating ``F`` does not add
new "information" into the system. First, let's develop a characterization on
the variables. Let ``z=\{z_i\}`` be the set of the highest order derivative variables, and
let ``\lambda = \{\lambda_i\}`` contain the rest of the variables. Note that ``z``
and ``\lambda`` must be disjoint. Therefore, DAEs can then be written as
```math
0 = F(z, \lambda, p, t) \\
```
Differentiating the above equation gives us
```math
F_z z' + F_\lambda \lambda' + F_t = 0
```
Note that ``z' = \{z'_i\}`` contains all the new terms generated by the
differentiation, as it contains variables with higher order derivatives than
before. Rearranging terms, we get
```math
F_z z' = -F_\lambda \lambda' - F_t.
```
When ``F_z`` is non-singular, we can use old terms to explicitly
solve for ``z'``, so we don't generate genuinely new equations, and this is the
*numerical consistency solvability criterion*. Therefore, we only add new
equations to the system if and only if ``F_z`` is singular. It's wasteful to
differentiate the entire system until the matrix is invertible; we can
differentiate a minimal subset of equations to make ``F_z`` non-singular.

Similarly, the *structural consistency solvability criterion* is then
```math
\mathfrak{I}(F, z) \text{ is structurally non-singular.}
```
For the above system, the differentiated system has
```math
\mathfrak{I}(F, \{x', y'\}) = \begin{pmatrix}
1 & 0 \\
1 & 1
\end{pmatrix},
```
which is structurally non-singular. Thus, the differentiated system satisfies
the structural consistency solvability criterion.

Note that the sparsity pattern of ``\mathfrak{I}(F, \{u_i'\}) + \mathfrak{I}(F,
\{u_i\})`` is always a subset of ``\mathfrak{I}(F, z)`` for arbitrary systems.
The structural consistency solvability criterion is stronger than the
structural integrability criterion, so we just check the structural consistency
solvability criterion.
