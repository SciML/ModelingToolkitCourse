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
    g(xᵢ) = f(xᵢ, x[i-1])
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

# ---------------------------------

x = zeros(10)
for i=2:10
    prob′ = remake(prob; u0=x[i], p=x[i-1])
    sol = solve(prob′, NewtonRaphson(); abstol=tol)
    x[i] = sol[]
end
plot(x; ylabel="x [m]", xlabel="time step")

# ---------------------------------

function du_dt(u,p,t)
    F, k, d = p
    x = u
    return (F - k*x^1.5)/d
end

prob = ODEProblem(du_dt, 0.0, (0.0, 0.01), [F, k, d])
sol = solve(prob)
plot(sol; xlabel="time [s]", ylabel="x [m]")

# ---------------------------------
tol = 1e-3
d=1
k=1000
F = 100

function du_dt1(u,p,t)
    F, k, d = p
    x, ẋ = u
    
    eqs = [
        ẋ                       # D(x) = ẋ
        (d*ẋ + k*x^1.5) - (F)   #    0 = ( lhs ) - ( rhs )
    ]

    return eqs
end

fmm = ODEFunction(du_dt1; mass_matrix=[1 0;0 0])
prob = ODEProblem(fmm, [0.0, F/d], (0.0, 0.01), [F, k, d])
sol = solve(prob; abstol=tol)
plot(sol; layout=2)


function du_dt1(du,u,p,t)
    F, k, d = p
    x, ẋ, ẍ = u
    
    du[1] = ẋ
    du[2] = ẍ
    du[3] = (d*ẋ + k*(x^1.5)) - (F)

end

fmm = ODEFunction(du_dt1; mass_matrix=[1 0 0;0 1 0;0 0 0])
prob = ODEProblem(fmm, [0.0, F/d, 0.0], (0.0, 0.01), [F, k, d])

sol = solve(prob, ImplicitEuler(autodiff=false); abstol=tol)
sol = solve(prob, ImplicitEuler(); abstol=tol, dt=0.001, adaptive=false)
sol = solve(prob, Rosenbrock23(); abstol=tol)

plot(sol; layout=3)

using ModelingToolkitComponents: SimpleImplicitEuler
sol = solve(prob, SimpleImplicitEuler(); abstol=tol)
plot(sol; layout=3)

function du_dt2(u,p,t)
    F, k, d = p
    x, ẋ, ẍ = u
    
    eqs = [
        ẋ                       # D(x) = ẋ
        ẍ                       # D(ẋ) = ẍ
        (d*ẍ + 1.5*k*(x^0.5)*ẋ) - (0)   #    0 = ( lhs ) - ( rhs )
    ]

    return eqs
end

fmm = ODEFunction(du_dt2; mass_matrix=[1 0 0;0 1 0;0 0 0])
prob = ODEProblem(fmm, [0.0, F/d, 0.0], (0.0, 0.01), [F, k, d])
sol = solve(prob, ImplicitEuler(); abstol=tol)



sol′ = solve(prob, RadauIIA5(); abstol=tol, reltol=10.0)

plot(sol′ ; layout=3)

plot(sol; idxs=1, xlabel="time [s]", ylabel="x [m]")
plot!(sol′; idxs=1, xlabel="time [s]", ylabel="x [m]")

plot(sol; idxs=2, xlabel="time [s]", ylabel="ẋ [m/s]")
plot!(sol′; idxs=2, xlabel="time [s]", ylabel="ẋ [m/s]")

plot(0:1e-5:0.01, x->ForwardDiff.derivative(t->sol(t), x)[2]; xlabel="time [s]", ylabel="ẍ [m/s^2]")
plot!(sol′; idxs=3, xlabel="time [s]", ylabel="ẍ [m/s^2]")
# ---------------------------------

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

vars = @variables x(t)=0.0 ẋ(t)=F/d ẍ(t)=0.0
eqs = [
    D(x) ~ ẋ
    D(ẋ) ~ ẍ
    d*ẋ + k*x^1.5 ~ F
]
@mtkbuild odesys = ODESystem(eqs, t, vars, [])
sys = structural_simplify(odesys)
prob = ODEProblem(sys, [], (0.0, 0.01))
sol = solve(prob; abstol=tol)
plot(sol; idxs=ẍ, xlabel="time [s]", ylabel="ẍ [m/s^2]")

ẍ_sol = sol[ẍ]
t_sol = sol.t


# ---------------------------------
@connector MechanicalPort begin
    v(t)
    f(t), [connect = Flow]
end

@mtkmodel Mass begin
    @parameters begin
        m = 10
    end
    @variables begin
        v(t) = 0
        f(t) = 0
    end
    @components begin
        port = MechanicalPort(;v=v, f=f)
    end
    @equations begin
        # connectors
        port.v ~ v
        port.f ~ f
        
        # physics
        f ~ m*D(v)
    end
end

@mtkmodel Damper begin
    @parameters begin
        d = 1
    end
    @variables begin
        v(t) = 0
        f(t) = d*v
    end
    @components begin
        port = MechanicalPort(;v=v, f=f)
    end
    @equations begin
        # connectors
        port.v ~ v
        port.f ~ f
        
        # physics
        f ~ d*v
    end
end

@mtkmodel System begin
    @parameters begin
        v
        m
        d
    end
    @components begin
        mass = Mass(;v,m)
        damper = Damper(;v,d)
    end
    @equations begin
        connect(mass.port, damper.port)
    end
end

@mtkbuild sys = System(;v=100, m=5, d=3)
full_equations(sys)

defs = ModelingToolkit.defaults(sys)
defs[sys.damper.d]





prob1 = ODEProblem(sys, [], (0,10))
sol1 = solve(prob1)
prob2 = remake(prob1; p=[sys.m=>3, sys.d=>4])
sol2 = solve(prob2)
plot(sol1; idxs=sys.mass.v)
plot!(sol2; idxs=sys.mass.v)





















using DAE2AE
aesys = dae_to_ae(odesys, Δt, 2)
sys = no_simplify(aesys)
prob = ODEProblem(sys, [], (0.0, 0.01))
sol = solve(prob, ImplicitEuler(nlsolve=NLNewton(check_div=false)); abstol=tol, dt=Δt, adaptive=false, initializealg=NoInit())
plot(sol; idxs=ẍ, xlabel="time [s]", ylabel="ẍ [m/s^2]")

ẍ_sol = sol[ẍ]
t_sol = sol.t




# -------------------
# -------------------
# -------------------
# ---   PROJECT -----
# -------------------
# -------------------
# -------------------

Δt = 1e-3


function f(Xᵢ, P)

    Xᵢ₋₁, Xᵢ₋₂ = P

    xᵢ, ẋᵢ, ẍᵢ = Xᵢ
    xᵢ₋₁, ẋᵢ₋₁, ẍᵢ₋₁ = Xᵢ₋₁
    xᵢ₋₂, ẋᵢ₋₂, ẍᵢ₋₂ = Xᵢ₋₂

    eqs = [
        ( ẋᵢ ) - ( (3*xᵢ - 4*xᵢ₋₁ + xᵢ₋₂)/(2*Δt) )
        ( ẍᵢ ) - ( (3*ẋᵢ - 4*ẋᵢ₋₁ + ẋᵢ₋₂)/(2*Δt) )
        ( d*ẋᵢ + k*xᵢ^1.5 ) - ( F )
    ]

    return eqs
end

T = (0:Δt:Δt*99)
X = zeros(100,3)
X[1,2] = F/d

prob = NonlinearProblem(f, X[1,:], (X[1,:],X[1,:]))
sol = solve(prob, NewtonRaphson())

for i=2:100
    prob′ = remake(prob; u0=X[i-1,:], p=(X[i-1,:],X[i>2 ? i-2 : i-1, :]))
    sol = solve(prob′, NewtonRaphson(); abstol=tol)
    X[i,:] = sol[:]
end
plot(X[:,1]; ylabel="x [m]", xlabel="time step")
plot(X[:,2]; ylabel="ẋ [m/s]", xlabel="time step")

plot(T, X[:,3]; ylabel="ẍ [m/s²]", xlabel="time [s]")
plot!(t_sol, ẍ_sol)























xo = zeros(11)
ẋo = zeros(11)
ẋo[1] = F/d
function du_dt3(u,p,t)
    F, k, d, Δt = p
    x, ẋ, ẍ = u
    
    i = round(Int, t/Δt)+1
    xo[i] = ForwardDiff.value(x)
    ẋo[i] = ForwardDiff.value(ẋ)

    eqs = [
        (ẋ) - (x-xo[i>1 ? i-1 : i])/Δt                       # D(x) = ẋ
        (ẍ) - (ẋ-ẋo[i>1 ? i-1 : i])/Δt                       # D(ẋ) = ẍ
        (d*ẋ + k*x^1.5) - (F)   #    0 = ( lhs ) - ( rhs )
    ]


    return eqs
end

fmm = ODEFunction(du_dt3; mass_matrix=[0 0 0;0 0 0;0 0 0])
prob = ODEProblem(fmm, [0.0, F/d, 0.0], (0.0, 0.01), [F, k, d, Δt])
sol′′ = solve(prob, ImplicitEuler(nlsolve=NLNewton(always_new=true, check_div=false, relax=4//10)); abstol=tol, adaptive=false, dt=Δt, initializealg=NoInit()) 
#
plot(sol′′ ; idxs=3, xlabel="time [s]", ylabel="ẍ [m/s^2]")



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
using ModelingToolkit: t_nounits as t, D_nounits as D

vars = @variables x(t)=0

eqs =[
    D(x) ~ (F - k*x^1.5)/d
]

@mtkbuild odesys = ODESystem(eqs, t, vars, [])

aesys = DAE2AE.dae_to_ae(odesys, Δt)
sys = structural_simplify(aesys)
prob = ODEProblem(sys,[], (0, Inf), [])
integrator = init(prob, ImplicitEuler(); adaptive=false, dt=Δt, abstol=tol, initializealg=NoInit())
step!(integrator)