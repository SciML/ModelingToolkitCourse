# ModelingToolkitCourse

## 18.S191 Special Subject in Mathematics: Composable System Modeling and Its Compilation

Traditionally, modeling physical systems often requires a deep understanding of the physics and equations of motions or states, simplifying the differential equations using conservation laws and constraints, and finally implementing simplified equations in a scientific computing language to numerically solve them. However, this workflow is tedious and not expressive. A simple change in the underlying physical system often requires a complete re-derivation of the simplified equations. A composable modeling system frees domain experts from the time-consuming derivation, simplification, and implementation by allowing them to model each physical component separately and hierarchically, thereby enabling them to build more accurate and complex models without compromising the simulation performance. In this course, we will dive into the practice of implementing composable physical models and the compilation process of the model system using the ModelingToolkit.jl acausal modeling framework in Julia. Students will learn the mathematics and numerical methods behind solving industry-scale models, covering topics such as differential-algebraic equations (DAEs), modern techniques in implicit integrators (backwards differentiation formulae (BDFs)), symbolic manipulation of equations via techniques like Pantelides algorithm and tearing of nonlinear systems, and more. Applications for solving real-world problems like modeling battery systems of electric vehicles, efficient control of hydraulic and HVAC systems, and more will be used to demonstrate how these techniques are used in industrial settings.

## Syllabus

Prerequisites: While this course will be mixing ideas from symbolic computing, numerical differential equations, and topics from mechanical engineering, no one in the course is expected to have covered all of these topics before. Understanding of calculus, linear algebra, and programming is essential. The course is considered self-contained starting from the basic building blocks of undergraduate differential equations. Problem sets will involve use of Julia, a Matlab-like environment (little or no prior experience required; you will learn as you go), for doing acausal modeling via the ModelingToolkit.jl system.

Textbook & Other Reading: There is no textbook for this course. For a textbook that covers the practical parts of doing modeling and simulation in an acausal way, Michael Tiller's "Modelica by Example" is a good reference (see https://mbe.modelica.university/). For a textbook that covers the algorithms of acausal modeling compilers, there is no recommended textbook and lecture notes will be supplied as a primary source.

Grading: The final project proposal (due January 15th) is 25%, and 75% for the final project (due February 2nd). Final projects will be submitted electronically via email.

## Final Project

The final project can take two forms: 

1. Developing an acausal model of some real-world system. 
2. Implementation and analysis of a new acausal modeling compiler feature. 
3. Implementation and analysis of numerical methods for acausal models.

A final project proposal is due January 19th and the final project is due on the last day of the course. The last day will be final project presentations where the work is demonstrated to the class.

The final project's deliverable can take two different forms:

1. A final project writeup: a 5-10 page paper using the style template from the [_SIAM Journal on Numerical Analysis_](http://www.siam.org/journals/auth-info.php) (or similar) explaining what was done, along with a Github repository package with the components of the model and docs/tests which demonstrate the successful composed model.
2. A pull request to one of the libraries (ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq, NonlinearSolve, etc.). For this version of the project, it is sufficient to supply a pull request to MTK/MSL with a description of the feature being implemented, tests of the transformation, and documentation showcasing its correct action on test models.

We expect the work to be roughly the same for the two routes, where the 1st would entail more theory and mathematical writeup while the latter is more focused on code and documentation. Note that any project considering doing a new acausal modeling feature should heavily consider doing the pull request route as writing a toy acausal modeling compiler within the timeframe of the course is likely to be unsuccessful.

Note that many of these projects are starter projects towards publications. If you're interested in continuing this work after the IAP towards a publication, please discuss during the project selection page so the project can be appropriately scoped.

### Project Type 1 Ideas: Developing an Acausal Model of A Real-World System

The following sources can be used as inspiration:

* Modelica "other" libraries (https://modelica.org/libraries/) 
* Modelica Standard Library (https://github.com/modelica/ModelicaStandardLibrary)


### Project Type 2 Ideas: Implementation and analysis of a new acausal modeling compiler feature

* Automated Laplace and Fourier transforms
* Automated function transformation of observables (i.e. log-transform states to enforce positivity)
* Symbolic generation of sensitivity analysis equations (https://github.com/SciML/ModelingToolkit.jl/issues/39)
* Lamperti transformation of stochastic differential equations (https://github.com/SciML/ModelingToolkit.jl/issues/140)
* Automated conversion of distributed delay equations into ODEs (https://github.com/SciML/ModelingToolkit.jl/issues/45)
* Specialized nonlinear solvers based on strongly connected components
* Inline integration (https://people.inf.ethz.ch/fcellier/Pubs/OO/esm_95.pdf)
* Automated detection of events from discontinuities in the ODE/DAE definition
* Polynomial chaos expansions for fast uncertainty quantification
* DCP on OptimizationSystem to automatically transform nonlinear optimization problems to convex optimization problems (http://cvxr.com/cvx/doc/dcp.html)
* Common subexpression elimination in Symbolic code generation
* Extendable C code generation maps from Symbolics
* Direct-quadrature-zero transformation for multibody systems and robotics (https://en.wikipedia.org/wiki/Direct-quadrature-zero_transformation)
* Pryce's algorithm for DAE index reduction (https://link.springer.com/article/10.1023/A:1021998624799, https://inria.hal.science/hal-03104030v2/document)

### Project Type 3 Ideas: Implementation and analysis of numerical methods for acausal models

* Adaptive order Radau methods (https://www.sciencedirect.com/science/article/pii/S037704279900134X)
* Parallel Rosenbrock and FIRK methods
* Handling the difficulties of BDFs in DAE systems (i.e. handling known deficiencies in the DFBDF algorithm)
* New time stepping schema for Rosenbrock methods for DAE interpolation performance
* Investigation of nonlinear solver globalization schemes for difficult DAE initialization problems

## Tentative Schedule

* January 10th: Introduction to the course, Guest lecture: Brad Carman, introduction to acausal modeling for physical systems with ModelingToolkit
* January 12th: Guest lecture: Brad Carman, developing high-fidelity models of hydraulic systems 
* January 15th: Martin Luther King Day!
* January 17th: Real numerical methods for implicit equations and stiff ordinary differential equations (ODEs), i.e., Jacobian-free Newton-Krylov, adaptive time stepping, dense output, sparse automatic differentiation, event handling.
* January 18th (Make up day for MLK day): Continuing discussion of stiff ODEs and onto numerical methods for differential-algebraic equations (DAEs). Rosenbrock methods, Backwards-Differentiation Formulae (BDF), fully-implicit Runge-Kutta methods.
* January 19th: Finishing the discussion on stiff ODEs and DAEs. If time allows, discussion of handling inverse problems (parameter estimation), adjoint methods, uncertainty quantification, and the connections to reverse-mode AD.
* January 22nd: Discussion and interactive workshop on debugging difficult stiff ODE/DAE models (featuring Brad Carman and Yingbo Ma).
* January 24th: Guest Lecture: Yingbo Ma. How acausal model compilers work: index reduction. Pantelides algorithm, dummy derivatives, and demonstrations.
* January 26th: Guest Lecture: Yingbo Ma. How acausal model compilers work: Tearing of nonlinear systems and alias elimination.
* January 29th: Guest Lecture: Yingbo Ma. How acausal model compilers work: Loop rerolling, specialized optimizations for multibody systems, and other generated code robustness and performance optimizations.
* January 31st: TBD based on what is not sufficiently covered earlier in the course.
* February 2nd: Final project presentations!
