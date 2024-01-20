# Lecture 7: Numerical and Structural Characterizations for DAEs

Numerical solvers cannot solve all DAEs. Consider the DAE system
```math
F(\{x', y', z'\}, \{x, y, z\}, t) = \begin{pmatrix}
    x + y - \sin(t) \\
    z - \sin(t) \\
    z' - \cos(t)
\end{pmatrix} = 0.
```
Note that the second equation and the third equation are equivalent, and there
are not enough constraints to uniquely determine ``x`` and ``y``. Let's see how
numerical solvers and ModelingToolkit behave.
```@example l7
using DifferentialEquations, Sundials, ModelingToolkit, Plots, LinearAlgebra

function f!(out, du, u, p, t)
    # u[1]: x, du[1]: x'
    # u[2]: y, du[2]: y'
    # u[3]: z, du[3]: z'
    out[1] = u[1] + u[2] - sin(t)
    out[2] = u[3] - sin(t)
    out[3] = du[3] - cos(t)
end
prob = DAEProblem(f!, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], (0, 100.0), differential_vars=[false, false, true])
sol1 = solve(prob, IDA())
sol2 = solve(prob, DFBDF())
println(sol1.retcode, "\t", sol2.retcode)
```
As expected, numerical solvers exit early and they cannot take more than one
step. We call such systems non-integrable or singular. Let's try to use
ModelingToolkit to analyze this problem.

```@example l7
@variables t x(t) y(t) z(t)
D = Differential(t)
eqs = [
    x + y ~ sin(t)
    z ~ sin(t)
    D(z) ~ cos(t)
]
@named sys = ODESystem(eqs, t)
sys = complete(sys);
try structural_simplify(sys) catch e println(e) end
```
We can see that ModelingToolkit also errors, correctly identifying that the
system is singular.

Let's consider another DAE system
```math
F(\{x', y'\}, \{x, y\}, t) = \begin{pmatrix}
    x - \sin(t) \\
    x' + y' - \cos(t)
\end{pmatrix} = 0.
```
If we plug the first equation into the second equation, we have
```math
(\sin(t))' + y' - \cos(t) = y' = 0.
```
Thus, the solution for ``y`` is just a constant function. Let's try to solve
this simple DAE using a numerical solver.
```@example l7
function f!(out, du, u, p, t)
    # u[1]: x, du[1]: x'
    # u[2]: y, du[2]: y'
    out[1] = u[1] - sin(t)
    out[2] = du[1] + du[2] - cos(t)
end
prob = DAEProblem(f!, [1, 0.0], [0.0, 0.0], (0, 100.0), differential_vars=[true, true])
sol1 = solve(prob, IDA())
sol2 = solve(prob, DFBDF())
println("[sol1: ", sol1.retcode, ": y(100)=", sol1[2, end], " steps: ", length(sol1.t), "]",
"\n", "[sol2 ", sol2.retcode, ": y(100)=", sol2[2, end], " steps: ", length(sol2.t), "]")
```
Note that we set ``y(0) = 0``, so the analytic solution is ``y(t) = 0``. Also,
superficially, we can have the initial condition
```math
x'(0) = x(0) = y(0) = 0, y'(0) = 1.
```
However, if we differentiate the first equation once, we get the hidden
constraint
```math
x'(t) - \cos(t) = 0.
```
Thus, the true consistent initial condition is
```math
x'(0) = 1, x(0) = y(0) = y'(0) = 0.
```
Let's replace the first equation with the differentiated equation and solve it
numerically,
```@example l7
function g!(out, du, u, p, t)
    # u[1]: x, du[1]: x'
    # u[2]: y, du[2]: y'
    out[1] = du[1] - cos(t)
    out[2] = du[1] + du[2] - cos(t)
end
prob = DAEProblem(g!, [1, 0.0], [0.0, 0.0], (0, 100.0), differential_vars=[true, true])
sol1_diff = solve(prob, IDA())
sol2_diff = solve(prob, DFBDF())
println("[sol1_diff: ", sol1_diff.retcode, ": y(100)=", sol1_diff[2, end], " steps: ", length(sol1_diff.t), "]",
"\n", "[sol2_diff ", sol2_diff.retcode, ": y(100)=", sol2_diff[2, end], " steps: ", length(sol2_diff.t), "]")
```
We can see that it takes far fewer iterations to solve the system, and the
numerical solution is much closer to the analytical solution ``y(t) = 0``.
If we check the residual of the original constraint, we get
```@example l7
plot(sol1_diff.t, sol1_diff[1, :] - sin.(sol1_diff.t), lab = "IDA")
plot!(sol2_diff.t, sol2_diff[1, :] - sin.(sol2_diff.t), lab = "DFBDF")
```
We see a significant numerical drift from the original constraint for both DAE
solvers. Again, let's see how ModelingToolkit does.

```@example l7
@variables t x(t) y(t)
D = Differential(t)
eqs = [
    x ~ sin(t)
    D(x) + D(y) ~ cos(t)
]
@named sys = ODESystem(eqs, t)
sys = complete(sys)
model = structural_simplify(sys)
prob = ODEProblem(model, [x=>0.0, y=>0.0, D(x)=>1.0, D(y)=>0.0], (0, 100.0))
sol = solve(prob, Rodas5P())
println("[sol: ", sol.retcode, ": y(100)=", sol[y, end], " steps: ", length(sol.t), "]")
```
```@example l7
norm(sol[x, :] - sin.(sol.t))
```
We can see that ModelingToolkit handles this DAE system perfectly. ``y(t)`` is
completely accurate, the original constraints are satisfied without drift, and
the numerical solver needs to take minimal steps. In my following lectures, we
will study how ModelingToolkit handles this example system.

## Numerical Integrability Criterion for DAEs

There are both differential equations and algebraic equations in acausal models.
Thus, a generic acausal model is a system of differential-algebraic equations
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
iteratively solving the best approximating linear equations (i.e., the Jacobian
of the nonlinear function with respect to the unknowns) to refine an initial
guess. By the chain rule, we have
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
    call an entry structural nonzero if it is not a structural zero. We also
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

    For convenience, we will simply use ``0`` and ``\mathbb{R}`` when structural
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
    there exist permutation matrices ``P`` and ``Q`` such that all the diagonal
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
      a_{i,\sigma(i)} \ne 0,
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
    bipartite graph ``(U, V, E)`` is a tuple of a set of source vertices ``U =
    \{1, 2, ..., m\}``, a set of destination vertices ``V\in \{1, 2, ..., n\}``,
    and edges between them ``E \subseteq U\times V`` defined by
    ```math
    \begin{equation}
    E = \{(i, j): A_{i,j} \ne 0\}.
    \end{equation}
    ```
    Similarly, the induced bipartite graph of a sparse matrix is the induced
    bipartite graph of its induced incidence matrix.

!!! definition "Bipartite Matching"

    A matching of a bipartite graph ``(U, V, E)`` is a set ``M \subseteq E``
    where every vertex in ``U`` and ``V`` can appear at most once in ``M``. A
    matching ``M`` is perfect if ``|M| = |U| = |V|``. We call an edge in a
    matching matched, otherwise, free. It is often more convenient to interpret
    matching as the function ``M: U \to (V \cup \emptyset)`` defined as
    ```math
    \begin{equation}
    m(i) = \begin{cases}
        j, & \text{if } (i, j) \in E \\
        \emptyset, & \text{else}
    \end{cases}.
    \end{equation}
    ```

!!! info "Structural Non-singularity and Perfect Matching Equivalence Theorem"

    A sparse matrix ``A\in\mathbb{R}^{m\times n}`` is structurally non-singular
    if and only if its induced bipartite graph has a perfect matching.

    ##### Proof:
    - ``\Leftarrow``: Suppose ``A\in\mathbb{R}^{m\times n}`` is a structurally
      non-singular sparse matrix, then ``m=n`` and there exists a row
      permutation ``\sigma\in\S_n`` such that ``\forall i\in\{1, ..., n\},
      A_{\sigma(i), i} \ne 0``. Note that ``\sigma`` is a perfect matching.
    - ``\Rightarrow``: Suppose the induced bipartite graph ``(U, V, E)`` of
      ``A\in\mathbb{R}^{m\times n}`` has a perfect bipartite matching ``\sigma``,
      then ``|U| = |V| = |\sigma| = m = n`` and ``\sigma\in\S_n``. In particular,
      ``\forall i\in\{1, ..., n\}, A_{\sigma(i), i} \ne 0``. Hence, ``A`` is
      structurally non-singular. ``\blacksquare``

!!! definition "Augmenting Path"

    Given a particular matching, an alternating path is a sequence of adjacent
    edges that alternate between being matched and free. In particular, an
    augmenting path is a alternating path that starts and ends with free edges.

!!! info "Bipartite Graph Maximum Cardinality Matching Theorem"

    A matching ``M`` of a bipartite graph has maximum cardinality matching if
    and only if there is no augmenting path with respect to ``M``.

    ##### Proof:
    - ``\Rightarrow``: We will show this using proof by contrapositive. Suppose
      a bipartite graph has a matching ``M`` and an augmenting path ``A``. Let
      ``\hat{M} = M \bigtriangleup A := (M \setminus A)\union (A\setminus M)``,
      then ``\hat{M}`` is a matching and ``|\hat{M}| = |M| + 1``. Thus, ``M`` is
      not a maximum cardinality matching.
    - ``\Leftarrow``: We will should this using proof by contrapositive, again.
      Suppose a bipartite graph ``(U, V, E)`` has a non-maximum cardinality
      matching ``B``, we want to seek an augmenting path. Let ``A`` be a
      maximum cardinality matching. We claim ``P = A\bigtriangleup B`` contains
      at least one augmenting path by the following arguments
        - Since all edges of ``P`` come from two matchings, by the definition of
          a matching, each vertex can have at most two edges. Therefore, ``P``
          contains either paths or cycles, and such segments are alternating
          between ``A`` and ``B``.
        - ``|P\cap A| > |P\cap B|``. Note that
          ``P\cap A = ((A\setminus B) \cup (B\setminus A)) \cap A = ((A\setminus B)\cap A)\cup ((B\setminus A) \cap A) = (A\setminus B)``,
          and similarly ``|P\cap B| = |B\setminus A|``. Thus,
          ``|P\cap A| = |A| - |A \cap B|`` and ``|P\cap B| = |B| - |A \cap B|``.
          By the maximality of ``A``, we know ``|A|>|B|``. Therefore, ``|P\cap
          A| > |P\cap B|``.
        - By the previous argument, there must be at least one connected
          component such that it contains more edges in ``A`` than ``B``. Since
          all cycles in ``P`` must be even length and alternating, such
          segment can only be a path, and in particular, an augmenting path
          with respect to ``B``.
      ``\blacksquare``

!!! info "Augmenting Path Algorithm for Finding a Maximum Cardinality Matching"

    Input: bipartite graph ``G = (U, V, E)``.

    #TODO: make this better
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
    cardinality matching ``M`` satisfies ``|M| = |U| = |V|``.



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
\mathfrak{I}(\mathfrak{I}(F, \{u_i'\}) + \mathfrak{I}(F, \{u_i\})) \text{ is
structurally non-singular.}
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
\mathfrak{I}(\mathfrak{I}(F, \{u_i'\}) + \mathfrak{I}(F, \{u_i\})) =
\begin{pmatrix}
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

Note that the sparsity pattern of ``\mathfrak{I}(\mathfrak{I}(F, \{u_i'\}) +
\mathfrak{I}(F, \{u_i\}))`` is always a subset of ``\mathfrak{I}(F, z)`` for
arbitrary systems.

!!! info "Structural Consistency Solvability Theorem"

    Structural consistency solvability implies structural integrability.
    Moreover, if a DAE system is not integrable, then structural consistency
    solvability cannot be achieved by differentiating any equation any number of
    times. [^2]


[^2]: Curious readers can read the original Pantelides paper for ideas to prove
    this.

## Pantelides Algorithm

We can use the Pantelides algorithm to efficiently convert a DAE system that has
structural integrability to a system with structural consistency solvability,
even if it initially lacks this property..
