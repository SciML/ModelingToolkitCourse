# Debugging difficult stiff ODE/DAE models

## Strategies
### 1. Is your model correct?
In the world of programming, debugging a model has got to be the most challenging because all equations must be solved together.  If any equation is wrong then not only will the model not solve, but there is very little that can be done to identify which equation is problematic.  Therefore the best that we can do is implement best practices to ensure the model is correct from the begining.  So 

**Use acausal modeling (i.e. ModelingToolkit.jl)**
As has been shown ModelingToolkit.jl will help in many ways with model definition.  One of the first programming practices that it enables is the DRY (Don't Repeat Yourself) principle.  By defining components once and reusing them, this helps reduce the chance of human error.  For example, when discovering a component level bug, it will be fixed at one source of truth and the fix will automatically propogate throughout.  

**Start small and verify components**
In using acausal modeling, the main focus for ensuring well definied models lies mainly at the component level.  Make sure to implement the rules of thumb discussed previously for number of equations and sign conventions.   Each component should have a well definied unit test.  When building your model start with the smallest subsystem possible and build from there.  Attempting to build a full system model before checking the pieces is doomed to fail, leaving little to no insight into what went wrong.  When a model fails to run, the error message will rarely give enough information to pinpoint the problem.  The best tool for debugging is taking small incremental steps which allows one to identify which change caused the problem.  

**Make sure equations match states**
It is not always the case, but for most models, the unsimplified system should give a match of equations and states.  Let's take the pendulum problem for example

```@example l6
using ModelingToolkit, DifferentialEquations, Plots
@parameters t
D = Differential(t)

pars = @parameters m = 1 g = 1 L = 1 Φ=0

vars = @variables begin
    x(t)=+L*cos(Φ)
    y(t)=-L*sin(Φ)
    dx(t)=0
    dy(t)=0
    λ(t) = 0
end

eqs = [
    D(x) ~ dx   
    D(y) ~ dy

    m*D(dx) ~ -λ*(x/L)
    m*D(dy) ~ -λ*(y/L) - m*g

    x^2 + y^2 ~ L^2 # algebraic constraint
]

@named pendulum = ODESystem(eqs, t, vars, pars)
nothing # hide
```

When we view the `ODESystem` we can see it has matching equations and states

```@repl
pendulum
```

Note: when using `@mtkbuild` then `structural_simplify` is automatically called and we therefore cannot see the unsimplify system.  Replace `@mtkbuild` with `@named` to generate an `ODESystem` without applying `structural_simplify`.


**Add compliance**
The pendulum problem as described above is derived assuming the following:

- a massless perfectly stiff and rigid string/rod connected to the mass
- a point mass
- a frictionless mechanism

If we attempt to solve this system we can see that it only solves up to the point that `x` crosses 0.

```@example l6
sys = structural_simplify(pendulum)
prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 10))
sol = solve(prob)# gives retcode: DtLessThanMin
plot(sol; idxs=[x,y])
```

The problem is rooted in the algebraic constraint which has `x^2`.  Having exponents (squares or square roots) can often cause issues with numerical solutions.  In this case the issue is that a unique solution cannot be found, `x` could be positive or negative.  There are different solutions to this problem, however lets consider the concept of adding compliance.  In reality is it really possible to have a massless, perfectly stiff and rigid string?  No.  Therefore let's consider adjusting the problem so the string has stiffness, which means we add `L` now as a variable.

```@example l6
pars = @parameters m = 1 g = 1 L_0 = 1 Φ=0 k=1e6

vars = @variables begin
    L(t)=L_0
    x(t)=+L*cos(Φ)
    y(t)=-L*sin(Φ)
    dx(t)=0
    dy(t)=0
    λ(t) = 0
end

eqs = [
    D(x) ~ dx   
    D(y) ~ dy

    m*D(dx) ~ -λ*(x/L)
    m*D(dy) ~ -λ*(y/L) - m*g

    x^2 + y^2 ~ L^2 # algebraic constraint

    λ ~ k*(L - L_0)
]

@named stiffness_pendulum = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(stiffness_pendulum)
prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 10))
sol = solve(prob)# Success
plot(sol; idxs=[x,y])
```

As can be seen, now we get a solution.  But is it correct?  


**Try `dae_index_lowering`**
In some cases we can apply `dae_index_lowering()` to further simplify the problem.  In this case ModelingToolkit.jl finds a better form of the equations which can be solved wihtout issue.

```@example l6
sys = structural_simplify(dae_index_lowering(pendulum))
prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 10))
ref = solve(prob)
plot(ref; idxs=x, label="dae_index_lowering")
plot!(sol; idxs=x, label="compliance")
```


**Design components with variable complexity**
In general this can be achieved with parameters.  For example, a *mass-spring-damper* system can easily become a *mass-damper* system by setting the spring stiffness to zero.  But in other cases we might want to *structurally* variable the complexity.  For example, the `ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible.Tube` component has 2 structural parameters:

- `N` for discretization
- `add_inertia` for including the wave equation

Based on the inputs of these structural parameters, the number of generated equations will be different.  Therefore, to start simple, one can set `N=0` and `add_inertia=false` to generate the simplest form of the problem.  Solving this case first and ensuring the model physical behavior is correct is a good best practice before attempting to increase the fidelity of the model.  


**Check the values of parameters**
Another possible cause of problems in your model can come not from the equations, but from the parameters that are supplied to the equations.  It's always a good idea to ensure your parameters match real life values to some degree, and to ensure human error is not factoring in, it can be a good idea to use units (note ModelingToolkit v9 will be enforcing units using Uniful.jl).  If you know all of your parameters are correct but still having issues, another debugging strategy is to reduce the energy input of your system.  Rather than starting at 100%, start at 10%.  This gives the model a better chance to solve and with a model solution this gives some insight to what the root cause problem might be.  For example, if working with a hydrualic system, turn the input pressure down to 10%.  If the model solves and you can see one of the pressure vessels is going into the negative, now you have a clue to a root cause problem.  

**Check acausal boundary conditions**
As discussed in Lecture 1, acausal connections always have a minumum of 2 variables.  Therefore, acausal input (or boundary condition) components will need to pay attention to what should be done to both variables.  As an example, refer to the hydraulic cylinder problem from Lecture 2 and consider the case where the position $x$ is supplied as the input boundary condition and the mass flow input $\dot{m}$ is set to an `Open()` boundary condition, thereby solving for $\dot{m}$ to give input $x$.

![example](../img/Example.svg)

We can assemble the problem as

```@example l6
import ModelingToolkitStandardLibrary.Mechanical.Translational as T
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B

include("volume.jl") # <-- missing Volume component from MTKSL (will be released in new version)

function MassVolume(solves_force = true; name)

    pars = @parameters begin
        A = 0.01 #m²
        x₀ = 1.0 #m
        M = 10_000 #kg
        g = 9.807 #m/s²
        amp = 5e-2 #m
        f = 15 #Hz
        p_int=M*g/A
        dx=0
        drho=0
        dm=0
    end
    vars = []
    systems = @named begin
        fluid = IC.HydraulicFluid(; density = 876, bulk_modulus = 1.2e9)
        mass = T.Mass(;v=dx,m=M,g=-g)
        vol = Volume(;area=A, x=x₀, p=p_int, dx, drho, dm) # <-- missing Volume component from MTKSL (will be released in new version)
        mass_flow = IC.Open(;p_int)
        position = T.Position(solves_force)
        position_input = B.TimeVaryingFunction(;f = t -> amp*sin(2π*t*f) + x₀)
    end

    eqs = [
        connect(mass.flange, vol.flange, position.flange)
        connect(vol.port, mass_flow.port)
        connect(position.s, position_input.output)
    ]

    return ODESystem(eqs, t, vars, pars; systems, name)
end

@named odesys = MassVolume()
nothing # hide
```

If we check the number of equations and states we see a mismatch!

```@repl l6
odesys
```

The reason for the mismatch is that the input boundary condition `Position()` needs to know what to do about the connection variable for force `f`.  In this problem, do we need a force introduced to the system to make the mass move as set by `Position()`?  The answer is no, the force causing the mass to move is already given by the hydraulic pressure and gravity.  If we look at the documentation for `Position()` we can see that it has a structural parameter `solves_force` which is defaulted to `true`.  Therefore, to assemble the proper system we set this to `false` and now have a properly defined system

```@repl l6
@named odesys = MassVolume(false)
```


### 2. My model is correct!  Now what?
Next step is to force a model solution.  It's still possible that something with the model is wrong, but the best way to know that is to see what the equations are outputting.  For example if the model is simulating negative pressure, but negative pressure is impossible, then this is a good clue of what is wrong with the model!  The strategies for forcing a model solve will come from a simple hydraulic system that is attempting to start a hydraulic cylinder at a high pressure differential.  See [ModelingToolkit Industrial Example](https://github.com/bradcarman/ModelingToolkitWebinar) for more information about the model.

```@example l6
@mtkmodel System begin
    @parameters begin
        res₁₊Cₒ = 2.7
        res₁₊Aₒ = 0.00094
        res₁₊ρ₀ = 1000
        res₁₊p′ = 3.0e7
        res₂₊Cₒ = 2.7
        res₂₊Aₒ = 0.00094
        res₂₊ρ₀ = 1000
        res₂₊p′ = 0
        act₊p₁′ = 3.0e7
        act₊p₂′ = 0
        act₊vol₁₊A = 0.1
        act₊vol₁₊ρ₀ = 1000
        act₊vol₁₊β = 2.0e9
        act₊vol₁₊direction = -1
        act₊vol₁₊p′ = act₊p₁′
        act₊vol₁₊x′ = 0.5
        act₊vol₂₊A = 0.1
        act₊vol₂₊ρ₀ = 1000
        act₊vol₂₊β = 2.0e9
        act₊vol₂₊direction = 1
        act₊vol₂₊p′ = act₊p₂′
        act₊vol₂₊x′ = 0.5
        act₊mass₊m = 100
        act₊mass₊f′ = 0.1(-act₊p₁′ + act₊p₂′)
        src₊p′ = 3.0e7
        snk₊p′ = 0
        dmp₊c = 1000
    end
    @variables begin
        act₊vol₁₊x(t) = act₊vol₁₊x′
        act₊vol₁₊r(t) = act₊vol₁₊ρ₀ * (1 + act₊vol₁₊p′ / act₊vol₁₊β)
        act₊vol₂₊x(t) = act₊vol₂₊x′
        act₊vol₂₊r(t) = act₊vol₂₊ρ₀ * (1 + act₊vol₂₊p′ / act₊vol₂₊β)
        act₊mass₊x(t) = 0
        act₊mass₊ẋ(t) = 0
        res₁₊ṁ(t) = 0
        res₂₊ṁ(t) = 0
        act₊vol₁₊ṙ(t) = 0
        act₊vol₂₊ṙ(t) = 0
    end
    @equations begin
        D(act₊vol₁₊x) ~ act₊vol₁₊direction * act₊mass₊ẋ
        D(act₊vol₁₊r) ~ act₊vol₁₊ṙ
        D(act₊vol₂₊x) ~ act₊vol₂₊direction * act₊mass₊ẋ
        D(act₊vol₂₊r) ~ act₊vol₂₊ṙ
        D(act₊mass₊x) ~ act₊mass₊ẋ
        D(act₊mass₊ẋ) ~ ((-act₊vol₂₊A * act₊vol₂₊β * (act₊vol₂₊ρ₀ - act₊vol₂₊r)) / (act₊vol₂₊direction * act₊vol₂₊ρ₀) + (-act₊vol₁₊A * act₊vol₁₊β * (act₊vol₁₊ρ₀ - act₊vol₁₊r)) / (act₊vol₁₊direction * act₊vol₁₊ρ₀) - dmp₊c * act₊mass₊ẋ) / act₊mass₊m
        0 ~ -src₊p′ + (-act₊vol₁₊β * (act₊vol₁₊ρ₀ - act₊vol₁₊r)) / act₊vol₁₊ρ₀ + 0.5res₁₊Cₒ * res₁₊ρ₀ * ((res₁₊ṁ / (res₁₊Aₒ * res₁₊ρ₀))^2)
        0 ~ snk₊p′ + (act₊vol₂₊β * (act₊vol₂₊ρ₀ - act₊vol₂₊r)) / act₊vol₂₊ρ₀ + 0.5res₂₊Cₒ * res₂₊ρ₀ * ((res₂₊ṁ / (res₂₊Aₒ * res₂₊ρ₀))^2)
        0 ~ -res₁₊ṁ + act₊vol₁₊A * act₊vol₁₊x * act₊vol₁₊ṙ + act₊vol₁₊A * act₊vol₁₊direction * act₊mass₊ẋ * act₊vol₁₊r
        0 ~ res₂₊ṁ + act₊vol₂₊A * act₊vol₂₊ṙ * act₊vol₂₊x + act₊vol₂₊A * act₊vol₂₊direction * act₊mass₊ẋ * act₊vol₂₊r
    end
end

@mtkbuild sys = System()
prob = ODEProblem(sys, [], (0, 1))
sol = solve(prob)
```

As can be seen, when attempting to solve we get an `Unstable` return code.  Let's explore strategies to fix the problem or find a forced numerical solution for debugging purposes.

**Initial Conditions**
First, let's check the initial conditions to see if at time 0 we are starting with zero residual for our algebraic equations.

```@example l6
eqs = equations(sys)
defs = ModelingToolkit.defaults(sys)
residuals = Float64[]
for eq in eqs
    if !ModelingToolkit.isdifferential(eq.lhs)
        push!(residuals, ModelingToolkit.fixpoint_sub(eq.rhs, defs))
    end
end
residuals
```

As can be seen, we have a problem with our first algebraic equation, the residual is not zero!  To solve this problem, ModelingToolkit v9 will be releasing a new feature to properly generate a non-linear system to calculate initial conditions.  We also can apply the [Initialization Schemes](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#Initialization-Schemes) provided from DifferentialEquations.jl.  The `BrownFullBasicInit` is the default algorithm used, and this did not work for our problem, so we will move to the `ShampineCollocationInit`.

```@example l6
dt = 1e-7
sol = solve(prob; initializealg=ShampineCollocationInit(dt))
```

The `ShampineCollocationInit` solves the initial conditions by essentially taking a small step forward in time and then updating the initial condtions with that solve.  If this doesn't work, we can instead do this manually. 

```@example l6
dt = 1e-7
prob = ODEProblem(sys, [], (0, dt))
sol = solve(prob, ImplicitEuler(nlsolve=NLNewton(check_div=false, always_new=true, relax=4/10, max_iter=100)); dt, adaptive=false)

# update u0 with the ImplicitEuler non-adaptive step
prob′ = ODEProblem(sys, sol[2], (0, 1))
sol = solve(prob′);
sol.retcode
```

As can be seen, now we have a successful solve.  We can see the change to the initial conditions is very minimal.  As can be seen, the solver needs the derivative terms to be offset by a small amount.

```@example l6
[println("$s $(round(x; digits=3)) -> $(round(y; digits=3))") for (s,x,y) in zip(states(sys), prob.u0, prob′.u0)];
```

Another strategy that can help is to offset any initial condtions from 0 by a small value.  

**Adjust tolerance**
Here we get a solve by increasing the `abstol` and `reltol` to very large values.  This is therefore understood to give us a very low resolution solution that is far from the true solution, but we can now at least see if the model is calculating generally correct values, at least with the correct sign.  Here we expect the `act₊mass₊ẋ` to be around -1 and that's exactly what we get.  

```@example l6
prob = ODEProblem(sys, [], (0, 1))
sol = solve(prob, ImplicitEuler(); abstol=10000, reltol=100.0, initializealg=NoInit())
plot(sol; idxs=sys.act₊mass₊ẋ)
```


**Turn off adaptivity**
Another strategy similar to adjusting tolerance is to turn off adaptivity.  This means we can no longer gaurantee tolerance and arrive at the most efficient numerical solution, but we can at least adjust the time step such that it's small enough to give a good solution, but large enough so that a reasonable amount of steps are taken.

```@example l6
prob = ODEProblem(sys, [], (0, 1))
sol = solve(prob, ImplicitEuler(nlsolve=NLNewton(check_div=false, always_new=true, relax=4/10, max_iter=100)); initializealg=NoInit(), adaptive=false, dt=1e-4)
plot(sol; idxs=sys.act₊mass₊ẋ)
```


- use a lower order solver
- check_div=false
- always_new=true
- relaxation (increase, decrease)
- tolerances (increase, decrease)
- timestep (increase, decrease)
- autodiff (true, false, analytical jacobian)

### 3. Experimental Strategies




