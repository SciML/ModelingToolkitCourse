using ModelingToolkit
using DifferentialEquations
using Symbolics

@parameters t
D = Differential(t)

# parameters -------
pars = @parameters begin
    r₀ = 876 #kg/s
    β = 1.2e9 #Pa
    A = 0.01 #m²
    x₀ = 1.0 #m
    m = 1000 #kg
    g = 9.807 #m/s²
    amp = 1e-3 #m
    f = 10 #Hz    
end

dt = 1e-4 #s
t_end = 0.2 #s
time = 0:dt:t_end

x_fun(t,amp,f) = amp*sin(2π*t*f) + x₀
ẋ_fun = build_function(expand_derivatives(D(x_fun(t,amp,f))), t, amp, f; expression=false)
ẍ_fun = build_function(expand_derivatives(D(ẋ_fun(t,amp,f))), t, amp, f; expression=false)


vars = @variables begin
    x(t) = x₀
    ẋ(t) = ẋ_fun(0.0, amp, f)
    ẍ(t) = ẍ_fun(0.0, amp, f)
    p(t) = m*g/A #Pa
    ṁ(t)
    r(t)=r₀*(1 + p/β)
    ṙ(t)=0.0
end 

# let ------
V = x*A
V̇ = ẋ*A

eqs = [
    D(x) ~ ẋ 
    D(ẋ) ~ ẍ
    D(r) ~ ṙ

    r ~ r₀*(1 + p/β)

    ṁ ~ ṙ*V + r*V̇
    m*ẍ ~ p*A - m*g

    ṁ ~ ẋ_fun(t,amp,f)*A*r
]

@named odesys = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(odesys)
prob = ODEProblem(sys, [], (0, t_end))

sol_ṁ = solve(prob)

eqs = [
    D(x) ~ ẋ 
    D(ẋ) ~ ẍ
    D(r) ~ ṙ

    r ~ r₀*(1 + p/β)

    ṁ ~ ṙ*V + r*V̇
    m*ẍ ~ p*A - m*g

    x ~ x_fun(t)
]

@named odesys = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(odesys)
prob = ODEProblem(sys, [], (0, t_end))
sol_x = solve(prob)


using Plots
plot(sol_ṁ; idxs=r)
plot!(time, sol_x(time)[r])

plot(sol_ṁ; idxs=p/1e5)
plot!(time, sol_x(time)[p]/1e5)

plot(sol_ṁ; idxs=x)
plot!(time, sol_x(time)[x])

plot(sol_ṁ; idxs=ṁ)
plot!(time, sol_x(time)[ṁ])


















