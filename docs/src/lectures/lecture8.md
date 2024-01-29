# Tearing

!!! definition "Induced Directed Graphs"

    The induced directed graph ``G = (V, E_v)`` from the destination vertices of
    a bipartite graph ``(U, V, E)`` and a perfect matching ``M`` is defined as
    ```math
    E_v = \{(i, j): (M(i), j) \in E\}.
    ```
    Similarly, the induced directed graph ``G = (V, E_u)`` from the source
    vertices is defined as
    ```math
    E_u = \{(i, j): (i, M^{-1}(j)) \in E\}.
    ```

!!! definition "Strongly Connected Component"

    A strongly connected component ``c\subseteq E`` of a directed graph ``G = (V,
    E)`` is a maximum cardinality set of vertices such that any pair ``i \in c,
    j \in c``, there exists a path ``i \rightsquigarrow j`` in ``G``.

!!! note "Strongly Connected Components Uniqueness Theorem"

    The strongly connected components are unique for induced directed graphs
    from bipartite graphs with a perfect matching.

    #### Proof:
    See [^DulmageMendelsohn1958].

[^DulmageMendelsohn1958]: Dulmage, Andrew L., and Nathan S. Mendelsohn.
    "Coverings of bipartite graphs." Canadian Journal of Mathematics 10 (1958):
    517-534.


#### Example

Consider a nonlinear system represented by
```math
\begin{align}
f_1(v_1, v_3) &= 0 \\
f_2(v_1, v_3) &= 0 \\
f_3(v_1, v_2) &= 0
\end{align}
```
The incidence matrix is
```math
\begin{pmatrix}
1 & 0 & 1 \\
1 & 0 & 1 \\
1 & 1 & 0
\end{pmatrix}
```
where a perfect matching ``m`` is then defined as
```math
1 \mapsto 1, 2 \mapsto 3, 3 \mapsto 2.
```
The permuted matrix is then
```math
\begin{pmatrix}
1 & 0 & 1 \\
1 & 1 & 0 \\
1 & 0 & 1
\end{pmatrix}.
```
The matching can be interpreted as a solvability assignment, i.e.
```math
v_1 = \hat{f}_1(v_3)
v_2 = \hat{f}_3(v_1)
v_3 = \hat{f}_2(v_2)
```
Even if all ``\{\hat{f}_i\}`` are symbolically solvable, the above assignment
will not work because the interdependence of variables. The strongly connected
component definition captures this idea well. Variables in a non-trivial
strongly connected component are the largest set of variables that are
interdependent. The strongly connected components of the above system are
``\{\{1, 3\}, \{2\}\}``. Thus, by the previous matching, we should reorder the
equations as ``e_1, e_2, e_3`` and variables as ``v_1, v_3, v_2`` to isolate the
interdependent part. The resulting system is then
```math
\begin{pmatrix}
1 & 1 & 0 \\
1 & 1 & 0 \\
1 & 0 & 1
\end{pmatrix}.
```
Note that resulting matrix is block lower triangular and this is not a
coincidence. We can always reorder the system to be block lower triangular
granted by the following theorem.

!!! definition "Condensation Graph"

    A condensation graph ``G_c = (V_c, E_c)`` of a directed graph
    ``G = (V, E)`` is a directed graph that has vertices
    ```math
    V_c = \{\text{strongly connected components of } G\},
    ```
    and edges
    ```math
    E_c = \{(i, j): \exists i_e \in i, j_e \in j, i_e \ne j_e \land (i_e, j_e)
    \in E\}.
    ```

!!! note "Condensation Graphs are Acyclic Theorem"

    The condensation graph ``G_c = (V_c, E_c)`` induced from the directed graph
    ``G = (V, E)`` is acyclic.

    #### Proof:
    Suppose ``G_c`` is cyclic with a cycle consisting of vertices ``s = \{v_1,
    v_2, ...\} \subseteq V_c``. Then, any original vertices in ``v_i`` has a
    path to any vertices in ``v_j`` for all ``i, j``. Thus, ``G_c`` must have
    only one vertex. By the definition of a condensation graph, ``G_c`` must has
    no edges, and therefore, no cycles.

Since the condensation graph has no cycles, we can topologically sort strongly
connected components so that the resulting system is always block lower
triangular. Further, each block on the diagonal must be square, because the
perfect matching will map all variables in each strongly connected components to
distinct equations.

## Demo

```julia
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T

@parameters t
D = Differential(t)

function System(use_input, f; name)
    @parameters t

    pars = @parameters begin
        p_s = 200e5
        p_r = 5e5

        A_1 = 360e-4
        A_2 = 360e-4

        p_1 = 45e5
        p_2 = 45e5

        l_1 = 0.01
        l_2 = 0.05
        m_f = 250
        g = 0

        d = 100e-3

        Cd = 0.01

        m_piston = 880
    end

    vars = @variables begin
        ddx(t) = 0
    end

    systems = @named begin
        src = IC.FixedPressure(; p = p_s)
        valve = IC.SpoolValve2Way(; p_s_int = p_s, p_a_int = p_1, p_b_int = p_2,
            p_r_int = p_r, g, m = m_f, x_int = 0, d, Cd)
        piston = IC.Actuator(5;
            p_a_int = p_1,
            p_b_int = p_2,
            area_a = A_1,
            area_b = A_2,
            length_a_int = l_1,
            length_b_int = l_2,
            m = m_piston,
            g = 0,
            x_int = 0,
            minimum_volume_a = A_1 * 1e-3,
            minimum_volume_b = A_2 * 1e-3,
            damping_volume_a = A_1 * 5e-3,
            damping_volume_b = A_2 * 5e-3)
        body = T.Mass(; m = 1500)
        pipe = IC.Tube(5; p_int = p_2, area = A_2, length = 2.0)
        snk = IC.FixedPressure(; p = p_r)
        pos = T.Position()

        m1 = IC.FlowDivider(; p_int = p_2, n = 3)
        m2 = IC.FlowDivider(; p_int = p_2, n = 3)

        fluid = IC.HydraulicFluid()
    end

    if use_input
        @named input = B.SampledData(Float64)
    else
        @named input = B.TimeVaryingFunction(f)
    end

    push!(systems, input)

    eqs = [connect(input.output, pos.s)
        connect(valve.flange, pos.flange)
        connect(valve.port_a, piston.port_a)
        connect(piston.flange, body.flange)
        connect(piston.port_b, m1.port_a)
        connect(m1.port_b, pipe.port_b)
        connect(pipe.port_a, m2.port_b)
        connect(m2.port_a, valve.port_b)
        connect(src.port, valve.port_s)
        connect(snk.port, valve.port_r)
        connect(fluid, src.port, snk.port)
        D(body.v) ~ ddx]

    ODESystem(eqs, t, vars, pars; name, systems)
end

@named system = System(true, nothing)

# sys = structural_simplify(system)
using ModelingToolkit.StructuralTransformations, ModelingToolkit.BipartiteGraphs,
    Graphs
ts = TearingState(ModelingToolkit.expand_connections(system))
m = BipartiteGraphs.maximal_matching(ts.structure.graph, _->true, x->ts.structure.var_to_diff[x] === nothing);
count(x->x isa Int, m)
count(x->x===nothing, ts.structure.eq_to_diff)
ModelingToolkit.pantelides!(ts)
m = BipartiteGraphs.maximal_matching(ts.structure.graph, x->ts.structure.eq_to_diff[x]===nothing, x->ts.structure.var_to_diff[x] === nothing);
count(x->x isa Int, m)
count(x->x===nothing, ts.structure.eq_to_diff)
M = incidence_matrix(ts.structure.graph)
A = M[Int[m[i] for i in eachindex(m) if m[i] isa Int], Int[i for i in eachindex(m) if m[i] isa Int]]
all(isequal(1), diag(A))
g = BipartiteGraphs.DiCMOBiGraph{true}(complete(ts.structure.graph), complete(m));
scc = strongly_connected_components(g);
M[Int[m[i] for i in reduce(vcat, scc) if m[i] isa Int], Int[i for i in reduce(vcat, scc) if m[i] isa Int]]
for c in scc
    length(c) > 1 || continue
    B = M[Int[m[i] for i in c if m[i] isa Int], Int[i for i in c if m[i] isa Int]]
    display(B)
end
```
