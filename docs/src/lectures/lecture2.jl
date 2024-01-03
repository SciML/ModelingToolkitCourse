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
        m*ẍ ~ p*A - m*g
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
    x ~ x_fun(t)
]

@named odesys_x = ODESystem(eqs_x, t, vars, pars)
sys_x = structural_simplify(odesys_x)
prob_x = ODEProblem(sys_x, [], (0, t_end))
sol_x = solve(prob_x; saveat=time)

@named odesys_ṁ1 = ODESystem(eqs_ṁ1, t, vars, pars)
sys_ṁ1 = structural_simplify(odesys_ṁ1)
u0 = [sol_x[s][1] for s in states(sys_ṁ1)]
prob_ṁ1 = ODEProblem(sys_ṁ1, u0, (0, t_end))
sol_ṁ1 = solve(prob_ṁ1)

plot(sol_ṁ1; idxs=ṁ, label="guess", ylabel="ṁ")
plot!(sol_x; idxs=ṁ, label="solution")

@named odesys_ṁ2 = ODESystem(eqs_ṁ2, t, vars, pars)
sys_ṁ2 = structural_simplify(odesys_ṁ2)
prob_ṁ2 = ODEProblem(sys_ṁ2, u0, (0, t_end))
sol_ṁ2 = solve(prob_ṁ2)

plot(sol_x; idxs=x, label="solution", ylabel="x")
plot!(sol_ṁ1; idxs=x, label="case 1")
plot!(sol_ṁ2; idxs=x, label="case 2")

plot(time, (sol_ṁ1(time)[x] .- sol_ṁ2(time)[x])/1e-3, ylabel="x error [mm]", xlabel="t [s]")

# -----------------------------------------------------

using DataInterpolations
mass_flow_fun = LinearInterpolation(sol_x[ṁ], sol_x.t)


using ModelingToolkitStandardLibrary.Mechanical.Translational
using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible

@mtkmodel begin
    @parameters begin
        A = 0.01 #m²
        x₀ = 1.0 #m
        M = 10_000 #kg
        g = 9.807 #m/s²
        amp = 5e-2 #m
        f = 15 #Hz   
        p_int=M*g/A
    end
    @components begin
        mass = Mass(;m=M)
        vol = DynamicVolume(0, false; p_int, area=A, x_int=x₀, x_max=10*x₀)
        mass_flow = MassFlow(;p_int)
        mass_flow_input = TimeVaryingFunction(;f = mass_flow_fun)
    end
    @equations begin
        connect()
    end
    
end

