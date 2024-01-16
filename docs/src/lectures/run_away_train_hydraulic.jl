using ModelingToolkit
using DifferentialEquations
using ModelingToolkitStandardLibrary.Mechanical.Translational: MechanicalPort, Mass
using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible: Valve, DynamicVolume, HydraulicFluid, FixedPressure
using ModelingToolkitStandardLibrary.Blocks: Constant
using Plots

@parameters t
D = Differential(t)

# NOTE: How Decouple works to provide discontinuous behavior
@mtkmodel Decouple begin #Decouple(;k, d, v_a, v_b, x_a, x_b, f=0) | port_a, port_b
    @parameters begin
        k
        d
    end
    @variables begin
        v_a(t)
        v_b(t)
        x_a(t)
        x_b(t)
        f(t)=0
    end
    @components begin
        port_a = MechanicalPort(;v=v_a,f=+f)
        port_b = MechanicalPort(;v=v_b,f=-f)
    end
    @equations begin
        # differentials
        D(x_a) ~ v_a
        D(x_b) ~ v_b

        # connectors
        port_a.v ~ v_a
        port_b.v ~ v_b
        port_a.f ~ +f
        port_b.f ~ -f

        # physics
        f ~ ifelse(x_a >= x_b, (v_a - v_b)*d + k*(x_a - x_b), 0)
    end
end

@mtkmodel Mass begin #Mass(;m, v, f, a=f/m) | port
    @parameters begin
        m = 10
    end
    @variables begin
        v(t)
        f(t)
        a(t) = f/m
    end
    @components begin
        port = MechanicalPort(; v, f)
    end
    @equations begin
        # derivatives
        D(v) ~ a

        # connectors
        port.v ~ v
        port.f ~ f

        # physics
        f ~ m*D(v)
    end
end

@mtkmodel Damper begin #Damper(;d=1,v,f) | port_a, port_b
    @parameters begin
        d = 1
    end
    @variables begin
        v(t)
        f(t)
    end
    @components begin
        port_a = MechanicalPort()
        port_b = MechanicalPort()
    end
    @equations begin
        # connectors
        (port_a.v - port_b.v) ~ v
        port_a.f ~ +f
        port_b.f ~ -f

        # physics
        f ~ d*v
    end
end

@mtkmodel Spring begin #Spring(;k=100,x,v,f) | port_a, port_b
    @parameters begin
        k = 100
    end
    @variables begin
        x(t)
        v(t)
        f(t)
    end
    @components begin
        port_a = MechanicalPort()
        port_b = MechanicalPort()
    end
    @equations begin
        # derivatives
        D(x) ~ v

        # connectors
        (port_a.v - port_b.v) ~ v
        port_a.f ~ +f
        port_b.f ~ -f

        # physics
        f ~ k*x
    end
end

@mtkmodel Reference begin #Reference(;f) | port
    @parameters begin

    end
    @variables begin
        f(t)
    end
    @components begin
        port = MechanicalPort()
    end
    @equations begin
        # connectors
        port.v ~ 0
        port.f ~ -f
    end
end

# NOTE: How TrainCar works to track absolute position
@mtkmodel TrainCar begin #TrainCar(;m, v, x_a, x_b) | port_a, port_b
    @parameters begin
        m
        v
    end
    @variables begin
        x_a(t)
        x_b(t)
    end
    @components begin
        port_a = MechanicalPort(;v)
        port_b = MechanicalPort(;v)
        mass = Mass(;m,v)
    end
    @equations begin
        D(x_a) ~ port_a.v
        D(x_b) ~ port_b.v
        
        connect(mass.port, port_a, port_b)
    end
end

@mtkmodel Coupler begin #Coupler(;k,d,v) | port_a, port_b
    @parameters begin
        k
        d
        v
    end
    @variables begin
        
    end
    @components begin
        port_a = MechanicalPort(;v)
        port_b = MechanicalPort(;v)
        spring = Spring(;k,x=0,v=0,f=0)
        damper = Damper(;d,v=0,f=0)
    end
    @equations begin
        connect(port_a, spring.port_a, damper.port_a)
        connect(spring.port_b, damper.port_b, port_b)
    end
end

@mtkmodel Stopper begin #Stopper(;k,d) | port
    @parameters begin
        k
        d
    end
    @variables begin
        
    end
    @components begin
        port = MechanicalPort()
        spring = Spring(;k,x=0,v=0,f=0)
        damper = Damper(;d,v=0,f=0)
        ref = Reference()
    end
    @equations begin
        connect(port, spring.port_a, damper.port_a)
        connect(spring.port_b, damper.port_b, ref.port)
    end
end

# NOTE: How HydraulicStopper works
@component function HydraulicStopper(;name, area) #HydraulicStopper(;area)  | port
    pars = @parameters begin
        area = area
    end
    vars = @variables begin
        
    end
    systems = @named begin
        port = MechanicalPort()
        volume = DynamicVolume(5, true, false; p_int=0, area=0.1, x_int = 5, x_max = 5, x_min = 1, x_damp = 2, direction = -1)
        valve = Valve(;p_a_int=0, p_b_int=0, area_int=area, Cd=0.7)
        constarea = Constant(;k=area)
        open = FixedPressure(;p=0)
        fluid = HydraulicFluid()
    end
    eqs = [
        connect(port, volume.flange)
        connect(volume.port, valve.port_a)
        connect(valve.port_b, open.port)
        connect(constarea.output, valve.area)

        connect(fluid, volume.port)
    ]

    return ODESystem(eqs, t, vars, pars; systems, name)
end


@mtkmodel System begin
    @parameters begin
        v=10
    end
    @variables begin
        
    end
    @components begin
        car1 = TrainCar(;m=1000, v, x_a=0, x_b=0.9)
        cx1 = Coupler(;k=1e5,d=1e5,v)

        car2 = TrainCar(;m=1000, v, x_a=1.1, x_b=1.9)
        cx2 = Coupler(;k=1e5,d=1e5,v)

        car3 = TrainCar(;m=1000, v, x_a=2.1, x_b=2.9)
        cx3 = Coupler(;k=1e5,d=1e5,v)

        engine = TrainCar(;m=2000, v=10, x_a=3.1, x_b=3.9)
        decouple = Decouple(;k=1e6, d=1e6, v_a=v, v_b=0, x_a=3.9, x_b=10, f=0)
        stopper = HydraulicStopper(; area=1e-2)
    end
    @equations begin
        connect(car1.port_b, cx1.port_a)
        connect(car2.port_a, cx1.port_b)
        connect(car2.port_b, cx2.port_a)
        connect(car3.port_a, cx2.port_b)
        connect(car3.port_b, cx3.port_a)
        connect(engine.port_a, cx3.port_b)
        connect(engine.port_b, decouple.port_a)
        connect(decouple.port_b, stopper.port)
    end
end

@mtkbuild sys = System()

# NOTE: missing_variable_defaults()
prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 5))

# NOTE: strategies for challenging solve
sol = solve(prob, ImplicitEuler(nlsolve=NLNewton(check_div=false, always_new=true, relax=4/10, max_iter=100)); adaptive=false, dt=1e-5)

# NOTE: pressure wave and no negative pressure
plot(sol; 
        idxs=[sys.stopper.volume.v1.port.p/1e5, 
                sys.stopper.volume.v2.port.p/1e5, 
                sys.stopper.volume.v3.port.p/1e5, 
                sys.stopper.volume.v4.port.p/1e5, 
                sys.stopper.volume.v5.port.p/1e5],
        xlims=(0.72, 0.74),
        ylabel="pressure [bar]", 
        xlabel="time [s]")

# NOTE: pressure bump at the end when the valve closes (1m x_min)
plot(sol; 
        idxs=[sys.stopper.volume.v1.port.p/1e5, 
                sys.stopper.volume.v2.port.p/1e5, 
                sys.stopper.volume.v3.port.p/1e5, 
                sys.stopper.volume.v4.port.p/1e5, 
                sys.stopper.volume.v5.port.p/1e5],
        ylims=(0, 20),
        ylabel="pressure [bar]", 
        xlabel="time [s]")

plot(sol;
        idxs=[sys.engine.mass.v])

# NOTE: How to design safe acceleration limit?
g = 9.807
plot(sol;
        idxs=[
            sys.car1.mass.a/g,
            sys.car2.mass.a/g,
            sys.car3.mass.a/g,
            sys.engine.mass.a/g
        ],
        ylims=(-25,5),
        xlims=(0.7,0.85),
        ylabel="acceleration [g]", 
        xlabel="time [s]")
