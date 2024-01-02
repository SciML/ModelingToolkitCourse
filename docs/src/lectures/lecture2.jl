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
    m = 1000 #kg
    g = 9.807 #m/s²
    amp = 1e-3 #m
    f = 10 #Hz    
    ξ = 1.0
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

eqs = [
    D(x) ~ ẋ 
    D(ẋ) ~ ẍ
    D(r) ~ ṙ

    r ~ r₀*(1 + p/β)

    ṁ ~ ṙ*x*A + r*ẋ*A
    m*ẍ ~ p*A - m*g
]



eqs_ṁ = [
    eqs...
    ṁ ~ ẋ_fun(t,amp,f)*A*r*ξ
]

eqs_x = [
    eqs...
    x ~ x_fun(t)
]

@named odesys_x = ODESystem(eqs_x, t, vars, pars)
sys_x = structural_simplify(odesys_x)
prob_x = ODEProblem(sys_x, [], (0, t_end))
sol_x = solve(prob_x; saveat=time)

@named odesys_ṁ = ODESystem(eqs_ṁ, t, vars, pars)
sys_ṁ = structural_simplify(odesys_ṁ)
u0 = [sol_x[s][1] for s in states(sys_ṁ)]
prob_ṁ = ODEProblem(sys_ṁ, u0, (0, t_end))
sol_ṁ = solve(prob_ṁ)


plot(sol_ṁ; idxs=r)
plot!(sol_x; idxs=r)

plot(sol_ṁ; idxs=ṙ)
plot!(sol_x; idxs=ṙ)

plot(sol_ṁ; idxs=p/1e5)
plot!(sol_x, idxs=p/1e5)

plot(sol_ṁ; idxs=x)
plot!(sol_x; idxs=x)

plot(sol_ṁ; idxs=ẋ)
plot!(time, sol_x(time)[ẋ])

plot(sol_ṁ; idxs=ṁ)
plot!(time, sol_x(time)[ṁ])

plot(sol_ṁ; idxs=ṁ*0.67)
plot!(time, sol_x(time)[ṁ])

# =========================
sol_ṁ′ = solve(remake(prob_ṁ; p=[ξ=>0.67]))

plot(sol_ṁ′; idxs=r)
plot!(time, sol_x(time)[r])

plot(sol_ṁ′; idxs=ṙ)
plot!(time, sol_x(time)[ṙ])

plot(sol_ṁ′; idxs=p/1e5)
plot!(sol_ṁ; idxs=p/1e5)
plot!(time, sol_x(time)[p]/1e5)

plot(sol_ṁ′; idxs=x)
plot!(time, sol_x(time)[x])

plot(sol_ṁ′; idxs=ẋ)
plot!(time, sol_x(time)[ẋ])

plot(sol_ṁ′; idxs=ṁ)
plot!(time, sol_x(time)[ṁ])

plot(sol_ṁ′; idxs=ṁ)
plot!(time, sol_x(time)[ṁ])












