using ForwardDiff
using Plots

d=1      # damping coefficient [N/(m/s)]
k=1000   # spring stiffness [N/m]
Δt=1e-3  # time step [s]
F = 100  # input force [N]

function f(xᵢ, xᵢ₋₁)

    ẋᵢ = (xᵢ - xᵢ₋₁)/Δt     # finite difference derivative
    lhs = d*ẋᵢ + k*xᵢ^1.5   # lhs --> left hand side
    rhs = F                 # rhs --> right hand side

    return lhs - rhs     # equation --> lhs = rhs, residual --> 0 = lhs - rhs
end

# Newton's Method
# first time step (i=2)
xᵢ₋₁ = 0.0
xᵢ = xᵢ₋₁ #<-- guess
g(xᵢ) = f(xᵢ, xᵢ₋₁)  # g(xᵢ) turns f(xᵢ, xᵢ₋₁) into a function of only xᵢ
# Run Newton Iterations
xᵢ -= g(xᵢ)/ForwardDiff.derivative(g, xᵢ) # iteration 1
xᵢ -= g(xᵢ)/ForwardDiff.derivative(g, xᵢ) # iteration 2
xᵢ -= g(xᵢ)/ForwardDiff.derivative(g, xᵢ) # iteration 3


tol = 1e-3
x = zeros(10)
for i=2:10
    g(xᵢ) = f(xᵢ, x[i-1])
    Δx = Inf
    while abs(Δx) > tol
        Δx = g(x[i])/ForwardDiff.derivative(g, x[i])
        x[i] -= Δx
    end
end

plot(x; ylabel="x [m]", xlabel="time step")


using DifferentialEquations

p  = xᵢ₋₁ = 0.0 # initial condition if i=2, x[1]=0
u0 = xᵢ = xᵢ₋₁  # guess value for x[i]

prob = NonlinearProblem(f, u0, p)
sol=solve(prob, NewtonRaphson(); abstol=tol)

x = zeros(10)
for i=2:10
    prob′ = remake(prob; u0=x[i], p=x[i-1])
    sol = solve(prob′, NewtonRaphson(); abstol=tol)
    x[i] = sol[]
end
plot(x; ylabel="x [m]", xlabel="time step")



function du_dt(u,p,t)
    F, k, d = p
    x = u
    return (F - k*x^1.5)/d
end
u0 = 0.0            # initial value for x
p = [F, k, d]       # parameters
tspan = (0.0, 0.01) # solution time span
prob = ODEProblem(du_dt, u0, tspan, p)
sol = solve(prob)
plot(sol; xlabel="time [s]", ylabel="x [m]")


function du_dt(u,p,t)
    F, k, d = p
    x, ẋ = u

    eqs = [
        ẋ                       # D(x) = ẋ
        (d*ẋ + k*x^1.5) - (F)   #    0 = ( lhs ) - ( rhs )
    ]

    return eqs
end

fmm = ODEFunction(du_dt; mass_matrix=[1 0; 0 0])
u0 = [0.0, F/d] # initial value for x,ẋ
prob = ODEProblem(fmm, u0, tspan, p)
sol = solve(prob)
plot(sol; idxs=1, xlabel="time [s]", ylabel="x [m]")


function du_dt(u,p,t)
    F, k, d = p
    x, ẋ, ẍ = u

    eqs = [
        ẋ                       # D(x) = ẋ
        ẍ                       # D(ẋ) = ẍ
        (d*ẋ + k*(x^1.5)) - (F)   #    0 = ( lhs ) - ( rhs )
    ]

    return eqs
end

fmm = ODEFunction(du_dt; mass_matrix=[1 0 0;0 1 0;0 0 0])
u0 = [0.0, F/d, 0.0] # initial value for x, ẋ, ẍ
prob = ODEProblem(fmm, u0, tspan, p)
sol = solve(prob);
sol.retcode