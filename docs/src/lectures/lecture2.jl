using ModelingToolkit
using DifferentialEquations
using Symbolics
using Plots

@parameters t
D = Differential(t)

# parameters -------
pars = @parameters begin
    r₀ = 876 #kg/s
    β = 1.2e9 #Pa
    A = 0.01 #m²
    x₀ = 1.0 #m
    M = 10_000 #kg
    g = 9.807 #m/s²
    amp = 5e-2 #m
    f = 15 #Hz    
end

dt = 1e-4 #s
t_end = 0.2 #s
time = 0:dt:t_end

x_fun(t,amp,f) = amp*sin(2π*t*f) + x₀
ẋ_fun = build_function(expand_derivatives(D(x_fun(t,amp,f))), t, amp, f; expression=false)


vars = @variables begin
    x(t) # = x₀
    ẋ(t)
    ẍ(t)
    p(t) # = m*g/A #Pa
    ṁ(t)
    r(t)
    ṙ(t)
end 

function get_base_equations(density_type) 
    
    eqs = [
        D(x) ~ ẋ 
        D(ẋ) ~ ẍ
        D(r) ~ ṙ

        r ~ r₀*(1 + p/β)

        ṁ ~ ṙ*x*A + (density_type)*ẋ*A
        M*ẍ ~ p*A - M*g
    ]

    return eqs
end

eqs_ṁ1 = [
    get_base_equations(r₀)...
    ṁ ~ ẋ_fun(t,amp,f)*A*r
]

eqs_ṁ2 = [
    get_base_equations(r)...
    ṁ ~ ẋ_fun(t,amp,f)*A*r
]

eqs_x = [
    get_base_equations(r)...
    x ~ x_fun(t,amp,f)
]

@mtkbuild odesys_x = ODESystem(eqs_x, t, vars, pars)
prob_x = ODEProblem(sys_x, [], (0, t_end))
sol_x = solve(prob_x; saveat=time)

@mtkbuild odesys_ṁ1 = ODESystem(eqs_ṁ1, t, vars, pars)
u0 = [sol_x[s][1] for s in states(sys_ṁ1)]
prob_ṁ1 = ODEProblem(sys_ṁ1, u0, (0, t_end))
sol_ṁ1 = solve(prob_ṁ1)

plot(sol_ṁ1; idxs=ṁ, label="guess", ylabel="ṁ")
plot!(sol_x; idxs=ṁ, label="solution")

@mtkbuild odesys_ṁ2 = ODESystem(eqs_ṁ2, t, vars, pars)
prob_ṁ2 = ODEProblem(sys_ṁ2, u0, (0, t_end))
sol_ṁ2 = solve(prob_ṁ2)

plot(sol_x; idxs=x, label="solution", ylabel="x")
plot!(sol_ṁ1; idxs=x, label="case 1")
plot!(sol_ṁ2; idxs=x, label="case 2")

plot(time, (sol_ṁ1(time)[x] .- sol_ṁ2(time)[x])/1e-3, ylabel="x error [mm]", xlabel="t [s]")

# -----------------------------------------------------

using DataInterpolations
mass_flow_fun = LinearInterpolation(sol_x[ṁ], sol_x.t)

open("dm.jl","w") do io
    print(io, "u = [")
    join(io, sol_x[ṁ], ',')
    print(io, "]")
end


plot(sol_x; idxs=ṁ)
plot!(time, -263.9489641054512*cos.(2π*time*15))

import ModelingToolkitStandardLibrary.Mechanical.Translational as T
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B

function MassVolume(; name, dx, drho, dm)

    pars = @parameters begin
        A = 0.01 #m²
        x₀ = 1.0 #m
        M = 10_000 #kg
        g = 9.807 #m/s²
        amp = 5e-2 #m
        f = 15 #Hz   
        p_int=M*g/A
        dx=dx
        drho=drho
        dm=dm
    end
    vars = []
    systems = @named begin
        fluid = IC.HydraulicFluid(; density = 876, bulk_modulus = 1.2e9)
        mass = T.Mass(;v=dx,m=M,g=-g)
        vol = IC.Volume(;area=A, x=x₀, p=p_int, dx, drho, dm)
        mass_flow = IC.MassFlow(;p_int)
        mass_flow_input = B.TimeVaryingFunction(;f = mass_flow_fun)
    end

    eqs = [
        connect(mass.flange, vol.flange)
        connect(vol.port, mass_flow.port)
        connect(mass_flow.dm, mass_flow_input.output)
        connect(mass_flow.port, fluid)
    ]

    return ODESystem(eqs, t, vars, pars; systems, name)
end

dx = sol_x[ẋ][1]
drho = sol_x[ṙ][1]
dm = sol_x[ṁ][1]

@named odesys = MassVolume(; dx, drho, dm)
# using DAE2AE
sys = structural_simplify(odesys)
# sys = no_simplify(expand_connections(odesys)) |> complete
prob = ODEProblem(sys, [], (0, t_end))
sol=solve(prob)
plot(sol; idxs=sys.vol.x, linewidth=2)
plot!(sol_x; idxs=x)

plot(sol; idxs=sys.vol.dx, linewidth=2)
plot!(sol_x; idxs=ẋ)



du0 = [
    0 # mass₊v(t)
    sol_x[ẋ][1] # vol₊x(t)
    sol_x[ẋ][1] # vol₊moving_volume₊x(t)
    sol_x[ṙ][1] # vol₊moving_volume₊rho(t)
    sol_x[ṙ*β/r₀][1] # vol₊moving_volume₊port₊p(t)
    sol_x[ṙ*β/r₀][1] # vol₊damper₊port_b₊p(t)
    0 # vol₊moving_volume₊drho(t)
]
prob = DAEProblem(sys, du0, [], (0, t_end))
sol = solve(prob, DImplicitEuler(nlsolve=NLNewton(check_div=false)))

prob = ODEProblem(sys, [], (0, t_end))
sol = solve(prob, ImplicitEuler(nlsolve=NLNewton(check_div=false, max_iter=100, relax=4//10)); dt, adaptive=false)

using SimpleEuler
sol = solve(prob, BackwardEuler(); dt, adaptive=false)

plot(sol; idxs=sys.vol.x, linewidth=2)
plot!(sol_x; idxs=x)

plot(sol; idxs=sys.mass_flow.dm.u)

plot(sol; idxs=sys.vol.moving_volume.port.p)
