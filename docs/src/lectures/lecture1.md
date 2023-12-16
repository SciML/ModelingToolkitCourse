# Lecture 1: Introduction to acausal modeling for physical systems with ModelingToolkit.jl

## Background
### Julia
This course will use Julia as the fundamental tool to solve numerical problems.  ModelingToolkit.jl is a package writen in pure Julia and leverages the fundamental technologies of symbolic math from Symbolics.jl, numerical solvers from DifferentialEquations.jl, and automatic differentiation from ForwardDiff.jl.  To demonstrate an introduction to these technologeies, lets focus on one of the most fundamental engineering problems: the mass-spring-damper.  For now, let's leave the mass out of the system to avoid the 2nd derivative term

![](../img/spring_damper.svg)

This system can be represented by the ordinary differential equation (ODE):

```math
d*\dot{x} + k*x = F
```

To solve this in Julia we can apply finite differencing $\dot{x}_i = \frac{x_i - x_{i-1}}{\Delta t}$ and [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method).  Here we solve for the first time step...

```@example l1
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
```

Note that we can get the derivative for `f` from automatic differentiation using `ForwardDiff.derivative` (or using `ForwardDiff.jacobian` for a system of equations).  To solve for a series of time steps, we can simply update `x` and run again for each time step `Δt`.  

```@example l1
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
```


### DifferentialEquations.jl
For this simple problem it's easy enough to implement the Newton method and solve directly, however it's possible to instead use the solvers from DifferentialEquations.jl.  To do this, we simply need to defined a `NonlinearProblem` by supplying the function `f` of the form $f(u,p)$ where:

- $u$ is the variable (scalar or vector)
- $p$ is the parameters (scalar or vector)

In this case $u$ and $p$ corespond to `xᵢ` and `xᵢ₋₁`, respectfully.  This is refered to as the "out-of-place" form, where each call to `f` allocates, it is also possible to define $f(du,u,p)$ in "in-place" form that gives $du$ as a pre-allocated memory space to mutate.  

Then we can solve by specifying the method, in this case we specify `NewtonRaphson` to implement the same algorithm.

```@example l1
using DifferentialEquations

prob = NonlinearProblem(f, 0.0, 0.0)
sol=solve(prob, NewtonRaphson(); abstol=tol)
```

Note:  we get exactly the same result for the first time step.

TODO: now show how to use `remake` to solve for a series of time steps.

### Symbolics.jl
- Use Symbolics to build the function
- introduce @variables

### ModelingToolkit.jl

#### Solving a system of equations (`NonlinearSystem`)
- Use ModelingToolkit to solve the problem
- introduce @parameters
- show the f function
- show the jacobian
- introduce `equations`, `defaults`, and how to setup initial conditions and parameters

#### Using `structural_simplify()`
- demonstrate the problem simplificaiton

#### [Practice Problem]
- introduce remake
- plot a solution set (do something that errors once the initial conditions are not defaulted correctly)


#### Introducing Time: ODE's, DAE's (`ODESystem`)
- explain Differential
- explain Mass Matrix





## Acausal - Component Based Modeling
 



### Connections
- demonstrate the theory of connections
- TODO: reference to where nodal network modeling originated?
- thru variables sum
- accross variables are equal
- Question: does the across variable have to be velocity?  Could it be any other derivative?

### Components
- demonstrate how components are defined
- best practices

### [Practice Problem]
- solve the simple 2 component systems shown in connections help page

### Systems and Sub-Systems
- how to expose ports and create hierarchy

