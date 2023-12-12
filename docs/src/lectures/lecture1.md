# Introduction to acausal modeling for physical systems with ModelingToolkit.jl

## Background
### Julia
- solve a simple nonlinear problem with pure Julia
- demonstrate ForwardAD.jl and non-allocating code

### DifferentialEquations.jl
- Use NewtonRaphson to solve the problem

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

