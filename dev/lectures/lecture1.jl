using ForwardDiff
using Plots

d=1
k=1000
Δt=1e-3
F = 100

function f(xᵢ, xᵢ₋₁)

    ẋᵢ = (xᵢ - xᵢ₋₁)/Δt
    lhs = d*ẋᵢ + k*xᵢ^1.5
    rhs = F

    return lhs - rhs
end

# Newton's Method
# first time step (i=2)
xᵢ₋₁ = 0.0
xᵢ = xᵢ₋₁ #<-- guess
g(xᵢ) = f(xᵢ, xᵢ₋₁)
xᵢ -= g(xᵢ)/ForwardDiff.derivative(g, xᵢ)
xᵢ -= g(xᵢ)/ForwardDiff.derivative(g, xᵢ)
xᵢ -= g(xᵢ)/ForwardDiff.derivative(g, xᵢ)


# ---------------------------------

tol = 1e-3
x = zeros(10)
for i=2:10
    g(ξ) = f(ξ, x[i-1])
    Δx = Inf
    while abs(Δx) > tol
        Δx = g(x[i])/ForwardDiff.derivative(g, x[i]) 
        x[i] -= Δx
    end
end

plot(x; ylabel="x [m]", xlabel="time step")

# ---------------------------------
using DifferentialEquations

prob = NonlinearProblem(f, 0.0, 0.0)
sol=solve(prob, NewtonRaphson(); abstol=tol)


fe(x,p,t) = error_tracker(f, x, p, t)

function step_julia(xₒ, t)

    x = [xₒ]
    p = [Δt, d, k, xₒ]
    t += Δt
    g(x) = fe(x,p,t)

    while norm(g(x)) > tol
        x = x .- ForwardDiff.jacobian(g, x)\g(x)
    end

    return x[1], t
end

function solver_julia()
    xs = zeros(10)
    t = 0.0
    for i=2:10
        xs[i], t = step_julia(xs[i-1], t)
    end
    return xs
end

@time x_julia=solver_julia();
errors_julia = copy(errors)

plot(x_julia)

plot(iteration_error(errors_julia); yscale=:log10)
hline!([tol])

plot(iteration_number(errors_julia))




# DifferentialEquations.jl
using DifferentialEquations

xₒ=0.0
p = [Δt, d, k, xₒ]
prob = SteadyStateProblem(fe, [xₒ], p)
sol = solve(prob, NewtonRaphson(); abstol=tol)

function step_diffeq(xₒ, t)
    
    prob′ = remake(prob, p=[Δt, d, k, xₒ], u0=[xₒ])
    sol = solve(prob′, NewtonRaphson(); abstol=tol)

    for (i,er) in enumerate(errors)
        if isinf(er[1])
            errors[i] = (t, er[2], er[3])
        end
    end

    return sol.u[]
end

empty!(errors)
function solver_diffeq()
    xs = zeros(10)
    t = 0.0
    for i=2:10
        t += Δt
        xs[i] = step_diffeq(xs[i-1], t)
    end
    return xs
end

@time x_diffeq=solver_diffeq();
plot(x_diffeq)
plot!(x_julia)

errors_diffeq = copy(errors)

plot(iteration_error(errors_julia); yscale=:log10)
plot!(iteration_error(errors_diffeq); yscale=:log10)
hline!([tol])

plot(iteration_number(errors_julia))
plot!(iteration_number(errors_diffeq))

fode = ODEFunction(fe, mass_matrix=zeros(1,1))
prob = ODEProblem(fode, [0.0], (0, Inf), [Δt, d, k, xₒ])
empty!(errors)
integrator = init(prob, ImplicitEuler(nlsolve=NLNewton(always_new=false, max_iter=100)); adaptive=false, dt=Δt, abstol=tol, initializealg=NoInit())

function solver_ie()
    xs = zeros(10)
    
    for i=2:10
        integrator.p[end] = xs[i-1]
        step!(integrator)
        xs[i] = integrator.u[1]
    end
    return xs
end

x_ie = solver_ie()

errors_ie = copy(errors)

plot(x_diffeq)
plot!(x_julia)
plot!(x_ie)

plot(iteration_error(errors_julia); yscale=:log10)
plot!(iteration_error(errors_diffeq); yscale=:log10)
plot!(iteration_error(errors_ie); yscale=:log10)
hline!([tol])

plot(iteration_number(errors_julia))
plot!(iteration_number(errors_diffeq))
plot!(iteration_number(errors_ie))



using ModelingToolkit
using DAE2AE

@parameters t
D = Differential(t)
vars = @variables x(t)=0

eqs =[
    D(x) ~ (F - k*x^1.5)/d
]

@named odesys = ODESystem(eqs, t, vars, [])

aesys = DAE2AE.dae_to_ae(odesys, Δt)
sys = structural_simplify(aesys)
prob = ODEProblem(sys,[], (0, Inf), [])
integrator = init(prob, ImplicitEuler(); adaptive=false, dt=Δt, abstol=tol, initializealg=NoInit())
step!(integrator)