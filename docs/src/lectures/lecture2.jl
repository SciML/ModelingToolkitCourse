using ModelingToolkit
using OrdinaryDiffEq
using ForwardDiff

NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 1000, relax = 4 // 10)

@parameters t
D = Differential(t)

# parameters -------
r₀ = 876 #kg/s
β = 1.2e9 #Pa
A = 0.01 #m²
x₀ = 1.0 #m
m = 1000 #kg
g = 9.807 #m/s²
amp = 1e-3 #m
f = 10 #Hz
dt = 1e-4 #s
t_end = 0.2 #s

vars = @variables begin
    x(t) = x₀
    ẋ(t) = 2π*f*amp #m/s
    ẍ(t) = 0 #m/s²
    p(t) = m*g/A #Pa
    ṁ(t) = 0 #kg/s
    r(t) = r₀*(1 + p/β)
    ṙ(t) = 0
end 

# let ------
V = x*A
V̇ = ẋ*A

time = 0:dt:t_end

x_fun(t) = amp*sin(2π*t*f) + x₀
# ẋ_fun(t) = ForwardDiff.derivative(x_fun, t)

eqs = [
    D(x) ~ ẋ 
    D(ẋ) ~ ẍ
    D(r) ~ ṙ

    r ~ r₀*(1 + p/β)

    ṁ ~ ṙ*V + r*V̇
    m*ẍ ~ p*A - m*g

    x ~ x_fun(t)
]

@named odesys = ODESystem(eqs, t, vars, [])
sys = structural_simplify(odesys)
observed(sys)


prob = ODEProblem(sys, [], (0,t_end), [])

sol_i1 = solve(prob_i1, ImplicitEuler())