# Debugging difficult stiff ODE/DAE models

## Strategies
### 1. Is your model correct?
In the world of programming, debugging a model has got to be the most challenging because all the equations must be solved together.  If any equation is wrong then not only will the model not solve, but there is very little that can be done to identify which equation is problematic.  So what can be done?

- use acausal modeling
- start small and verify components
- design components so that complexity can be adjusted (for example, turn off fluid inertia).  Note the solve difference between `sol_ṁ1` and `sol_ṁ2`

- be very careful about boundary conditions: explain `Position(solves_force = true; name)` when to use `solves_force = false`

### 2. My model is correct!  Now what?
Next step is to force a model solution.  It's still possible that something with the model is wrong, but the best way to know that is to see what the equations are outputing.  For example if the model is simulating negative pressure, but negative pressure is impossible, then this is a good clue of what is wrong with the model!

- initial conditions
    - check that the equations give 0 residuals at time zero
    - use a small low order time step to get an improved initial condition
    - implement a small offset

- use non-adaptive solver 
- use a lower order solver
- check_div=false
- always_new=true
- relaxation (increase, decrease)
- tolerances (increase, decrease)
- timestep (increase, decrease)
- autodiff (true, false, analytical jacobian)

### 3. Experimental Strategies
- explore equation formulation: `no_simplify`, `structural_simplify(dae_index_reduction())`, `structural_simplify(; allow_parameter=false))`. Discuss the pendulum problem: https://github.com/SciML/ModelingToolkit.jl/issues/2417
- Brad's singular Jacobian avoidance trick: `z .= (J_ .+ tol) \ (du_ .+ tol)`



