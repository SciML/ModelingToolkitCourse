# Numerical and Structural Characterizations for DAEs

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
    x' - y
\end{pmatrix} = 0.
```
If we plug the first equation into the second equation, we have
```math
y(t) = (\sin(t))' = \cos(t).
```
Superficially, we can have the initial condition
```math
x'(0) = x(0) = y(0) = y'(0) = 0.
```
However, if we differentiate the first equation once, we get the hidden
constraint
```math
x'(t) - \cos(t) = 0.
```
Thus, the true consistent initial condition is
```math
x'(0) = y(0) = 1, x(0) = y'(0) = 0.
```

Let's try to solve this simple DAE using a numerical solver.
```@example l7
function f!(out, du, u, p, t)
    # u[1]: x, du[1]: x'
    # u[2]: y, du[2]: y'
    out[1] = u[1] - sin(t)
    out[2] = du[1] - u[2]
end
prob = DAEProblem(f!, [1, 0.0], [0.0, 1.0], (0, 100pi), differential_vars=[true, false])
sol1 = solve(prob, IDA())
sol2 = solve(prob, DFBDF())
println("[sol1: ", sol1.retcode, "]",
"\n", "[sol2 ", sol2.retcode, ": y(100pi)=", sol2[2, end], " steps: ", length(sol2.t), "]")
```

To better understand the numerical behavior, let's analyze the variable step
size behavior of the implicit Euler method of the original system. The local
truncation error for ``y`` is
```math
\begin{align}
&\text{lte} = \frac{\frac{y(t_n) - y(t_{n-1})}{h_n} - \frac{y(t_{n-1}) - y(t_{n-2})}{h_{n-1}}}
{h_{n} + h_{n-1}} (h_{n} + h_{n-1}) h_n \\
=& \left(\frac{y(t_n) - y(t_{n-1})}{h_n} - \frac{y(t_{n-1}) - y(t_{n-2})}{h_{n-1}}\right) h_n
\\
=& h_n\left(\frac{\frac{\sin(t_{n}) - \sin(t_{n-1})}{h_n} - \frac{\sin(t_{n-1}) - \sin(t_{n-2})}{h_{n-1}}}{h_n} -
\frac{\frac{\sin(t_{n-1}) - \sin(t_{n-2})}{h_{n-1}} - \frac{\sin(t_{n-2}) - \sin(t_{n-3})}{h_{n-2}}}{h_{n-1}}\right)
\\
=&
\frac{\sin(t_{n}) - \sin(t_{n-1})}{h_n} - \frac{\sin(t_{n-1}) - \sin(t_{n-2})}{h_{n-1}} -
h_n
\frac{\frac{\sin(t_{n-1}) - \sin(t_{n-2})}{h_{n-1}} - \frac{\sin(t_{n-2}) - \sin(t_{n-3})}{h_{n-2}}}{h_{n-1}}
\end{align}
```
Note that when ``h_{n} \to 0``, the local truncation error becomes
```math
\lim_{h_{n}\to 0} \text{lte} = \frac{\sin(t_{n}) - \sin(t_{n-1})}{h_n} - \frac{\sin(t_{n-1}) - \sin(t_{n-2})}{h_{n-1}}
= \cos(t_{n}) - \frac{\sin(t_{n-1}) - \sin(t_{n-2})}{h_{n-1}}
```
Note that this in general does not converge to ``0``. Thus, numerical solvers
could have difficulties in solving this system. Let's also confirm this
numerical behavior by setting the maximum order to ``1``.
```@example l7
solve(prob, DFBDF(max_order=Val(1)))
```

Let's replace the first equation with the differentiated equation and solve it
numerically,
```@example l7
function g!(out, du, u, p, t)
    # u[1]: x, du[1]: x'
    # u[2]: y, du[2]: y'
    out[1] = du[1] - cos(t)
    out[2] = du[1] - u[2]
end
prob = DAEProblem(g!, [1, 0.0], [0.0, 1.0], (0, 100pi), differential_vars=[true, false])
sol1_diff = solve(prob, IDA())
sol2_diff = solve(prob, DFBDF())
println("[sol1_diff: ", sol1_diff.retcode, ": y(100pi)=", sol1_diff[2, end], " steps: ", length(sol1_diff.t), "]",
"\n", "[sol2_diff ", sol2_diff.retcode, ": y(100pi)=", sol2_diff[2, end], " steps: ", length(sol2_diff.t), "]")
```
We can see that it takes far fewer iterations to solve the system, and the
numerical solution is much closer to the analytical solution ``y(100\pi) = 1``.
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
    D(x) ~ y
]
@named sys = ODESystem(eqs, t)
sys = complete(sys)
model = structural_simplify(sys)
prob = ODEProblem(model, [x=>0.0, y=>1.0, D(x)=>0.0, D(y)=>1.0], (0, 100pi))
sol = solve(prob, Rodas5P())
println("[sol: ", sol.retcode, ": y(100pi)=", sol[y, end], " steps: ", length(sol.t), "]")
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
      permutation ``\sigma\in S_n`` such that ``\forall i\in\{1, ..., n\},
      A_{\sigma(i), i} \ne 0``. Note that ``\sigma`` is a perfect matching.
    - ``\Rightarrow``: Suppose the induced bipartite graph ``(U, V, E)`` of
      ``A\in\mathbb{R}^{m\times n}`` has a perfect bipartite matching ``\sigma``,
      then ``|U| = |V| = |\sigma| = m = n`` and ``\sigma\in S_n``. In particular,
      ``\forall i\in\{1, ..., n\}, A_{\sigma(i), i} \ne 0``. Hence, ``A`` is
      structurally non-singular. ``\blacksquare``

!!! definition "Augmenting Path"

    Given a particular matching, an alternating path is a sequence of adjacent
    edges that alternate between being matched and free. In particular, an
    augmenting path is a alternating path that starts and ends with free vertices.

!!! info "Bipartite Graph Maximum Cardinality Matching Theorem"

    A matching ``M`` of a bipartite graph has maximum cardinality matching if
    and only if there is no augmenting path with respect to ``M``.

    ##### Proof:
    - ``\Rightarrow``: We will show this using proof by contrapositive. Suppose
      a bipartite graph has a matching ``M`` and an augmenting path ``A``. Let
      ``\hat{M} = M \bigtriangleup A := (M \setminus A)\cup (A\setminus M)``,
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

!!! info "Augmenting Path Algorithm"

    Input: bipartite graph ``g = (U, V, E)``, a vertex ``\text{vsrc} \in U``,
    and a partial matching.

    Output: return a boolean indicating the existence of an augmenting path,
    and if one is present, use the augmenting path to increase the cardinality
    of the partial matching by exactly one.

    In ModelingToolkit, there are the `𝑠neighbors(g, i)` function that returns a
    sorted list containing ``\{j: (i, j) \in E\}``, and the `𝑑neighbors(g, j)`
    function that returns a sorted list containing ``\{i: (i, j) \in E\}``.
    ModelingToolkit also encodes matching `M` using the `m::Matching` structure,
    let `j = m[i]`, it holds that `j::Int` if and only if ``(i, j) \in M`` and
    `j::Unassigned` if and only if ``(i, j) \not\in M``. It following code comes
    directly from ModelingToolkit.
    ```julia
    function construct_augmenting_path!(matching::Matching, g::BipartiteGraph, vsrc, dstfilter,
            dcolor = falses(ndsts(g)), scolor = nothing)
        scolor === nothing || (scolor[vsrc] = true)

        # if a `vdst` is unassigned and the edge `vsrc <=> vdst` exists
        for vdst in 𝑠neighbors(g, vsrc)
            if dstfilter(vdst) && matching[vdst] === unassigned
                matching[vdst] = vsrc
                return true
            end
        end

        # for every `vsrc` such that edge `vsrc <=> vdst` exists and `vdst` is uncolored
        for vdst in 𝑠neighbors(g, vsrc)
            (dstfilter(vdst) && !dcolor[vdst]) || continue
            dcolor[vdst] = true
            if construct_augmenting_path!(matching, g, matching[vdst], dstfilter, dcolor,
                scolor)
                matching[vdst] = vsrc
                return true
            end
        end
        return false
    end
    ```
    Note that the augmenting path algorithm never removes any matched vertices
    in ``U``.

!!! info "Augmenting Path Algorithm for Finding a Maximum Cardinality Matching"

    Input: bipartite graph ``g = (U, V, E)``.

    Output: matching ``M``.

    The following code comes directly from ModelingToolkit. Note that
    `𝑠vertices(g)` returns `1:n` where ``n=|U|``.
    ```julia
    function maximal_matching(g::BipartiteGraph, srcfilter = vsrc -> true,
            dstfilter = vdst -> true, ::Type{U} = Unassigned) where {U}
        matching = Matching{U}(ndsts(g))
        foreach(Iterators.filter(srcfilter, 𝑠vertices(g))) do vsrc
            construct_augmenting_path!(matching, g, vsrc, dstfilter)
        end
        return matching
    end
    ```
    Given that the augmenting path algorithm never removes any matched vertices
    in ``U``, and if no augmenting path starts in vertex ``i\in U``, then ``i``
    will never be matched using the augmenting path algorithm. It is sufficient
    to run the augmenting path algorithm for all vertices in ``U`` by the
    Bipartite Graph Maximum Cardinality Matching Theorem.

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

We can use the Pantelides algorithm [^Pantelides1988] to efficiently convert a
DAE system that has structural integrability to a system with structural
consistency solvability, even if it initially lacks this property.

The gist of the Pantelides algorithm is that, we can try to find an augmenting
path for all equation (source) vertices on the sub-graph that only contains
highest differentiated variable (destination) vertices, and if there is no
augmenting path starting at an equation vertex, then we can differentiate all
the equations and variables reached in the augmenting path search until there is
an augmenting path starting at the differentiated equation vertex. Note that if
we differentiate an equation in the form of
``f(x, y, ...)`` we get
```math
\frac{d}{dt}f(x, y, ...) = \frac{\partial f}{\partial x}x' + \frac{\partial f}{\partial y}y' + ...
```
Thus, the incidence of the differentiated equation is trivial to compute.

The following code comes directly from ModelingToolkit.
```julia
function pantelides!(state::TransformationState; finalize = true, maxiters = 8000)
    @unpack graph, solvable_graph, var_to_diff, eq_to_diff = state.structure
    neqs = nsrcs(graph)
    nvars = nv(var_to_diff)
    vcolor = falses(nvars)
    ecolor = falses(neqs)
    var_eq_matching = Matching(nvars)
    neqs′ = neqs
    nnonemptyeqs = count(eq -> !isempty(𝑠neighbors(graph, eq)) && eq_to_diff[eq] === nothing,
        1:neqs′)

    varwhitelist = computed_highest_diff_variables(state.structure)

    if nnonemptyeqs > count(varwhitelist)
        throw(InvalidSystemException("System is structurally singular"))
    end

    for k in 1:neqs′
        eq′ = k
        eq_to_diff[eq′] === nothing || continue
        isempty(𝑠neighbors(graph, eq′)) && continue
        pathfound = false
        # In practice, `maxiters=8000` should never be reached, otherwise, the
        # index would be on the order of thousands.
        for iii in 1:maxiters
            # run matching on (dx, y) variables
            #
            # the derivatives and algebraic variables are zeros in the variable
            # association list
            resize!(vcolor, nvars)
            fill!(vcolor, false)
            resize!(ecolor, neqs)
            fill!(ecolor, false)
            pathfound = construct_augmenting_path!(var_eq_matching, graph, eq′,
                v -> varwhitelist[v], vcolor, ecolor)
            pathfound && break # terminating condition
            if is_only_discrete(state.structure)
                error("The discrete system has high structural index. This is not supported.")
            end
            for var in eachindex(vcolor)
                vcolor[var] || continue
                if var_to_diff[var] === nothing
                    # introduce a new variable
                    nvars += 1
                    var_diff = var_derivative!(state, var)
                    push!(var_eq_matching, unassigned)
                    push!(varwhitelist, false)
                    @assert length(var_eq_matching) == var_diff
                end
                varwhitelist[var] = false
                varwhitelist[var_to_diff[var]] = true
            end

            for eq in eachindex(ecolor)
                ecolor[eq] || continue
                # introduce a new equation
                neqs += 1
                eq_derivative!(state, eq)
            end

            for var in eachindex(vcolor)
                vcolor[var] || continue
                # the newly introduced `var`s and `eq`s have the inherits
                # assignment
                var_eq_matching[var_to_diff[var]] = eq_to_diff[var_eq_matching[var]]
            end
            eq′ = eq_to_diff[eq′]
        end # for _ in 1:maxiters
        pathfound ||
            error("maxiters=$maxiters reached! File a bug report if your system has a reasonable index (<100), and you are using the default `maxiters`. Try to increase the maxiters by `pantelides(sys::ODESystem; maxiters=1_000_000)` if your system has an incredibly high index and it is truly extremely large.")
    end # for k in 1:neqs′

    finalize && for var in 1:ndsts(graph)
        varwhitelist[var] && continue
        var_eq_matching[var] = unassigned
    end
    return var_eq_matching
end
```

[^Pantelides1988]: Pantelides, Constantinos C. "The consistent initialization of
    differential-algebraic systems." SIAM Journal on scientific and statistical
    computing 9.2 (1988): 213-231.
