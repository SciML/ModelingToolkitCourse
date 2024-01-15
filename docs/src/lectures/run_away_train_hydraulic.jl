using ModelingToolkit
using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible: Valve, DynamicVolume, HydraulicFluid, FixedPressure
using ModelingToolkitStandardLibrary.Blocks: Constant

@mtkmodel Decouple begin #Decouple(;k, d, v_a, v_b, x_a, x_b, f=0)
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

@mtkmodel Mass begin #Mass(;m, v, f)
    @parameters begin
        m = 10
    end
    @variables begin
        v(t)
        f(t)
    end
    @components begin
        port = MechanicalPort()
    end
    @equations begin
        # connectors
        port.v ~ v
        port.f ~ f

        # physics
        f ~ m*D(v)
    end
end

@mtkmodel Damper begin #Damper(;d,v,f)
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

@mtkmodel Spring begin #Spring(;k,x,v,f)
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

@mtkmodel Reference begin
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

@mtkmodel TrainCar begin #TrainCar(;m, v, x_a, x_b)
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

@mtkmodel Stopper begin #Stopper(;k,d)
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

@mtkmodel Coupler begin #Coupler(;k,d,v)
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

@component function HydraulicStopper(;name, area)
    pars = @parameters begin
        area = area
    end
    vars = @variables begin
        
    end
    systems = @named begin
        flange = MechanicalPort()
        volume = DynamicVolume(5, true, false; p_int=0, area=0.1, x_int = 5, x_max = 5, x_min = 1, x_damp = 2, direction = -1)
        valve = Valve(;p_a_int=0, p_b_int=0, area_int=area, Cd=0.7)
        constarea = Constant(;k=area)
        open = FixedPressure(;p=0)
        fluid = HydraulicFluid()
    end
    eqs = [
        connect(flange, volume.flange)
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
        connect(decouple.port_b, stopper.flange)
    end
end

@mtkbuild sys = System()
prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 5))
sol = solve(prob, ImplicitEuler(nlsolve=NLNewton(check_div=false, always_new=true, relax=4/10, max_iter=100)); adaptive=false, dt=1e-5)
plot(sol; idxs=[sys.decouple.x_a, sys.decouple.x_b, (5-sys.stopper.volume.x)+10])
plot!(sol; idxs=(5-sys.stopper.volume.x)+10+sys.stopper.volume.v2.x)


plot(sol; idxs=[sys.stopper.valve.port_a.p/1e5], ylims=(0, 20), xlims=(0.5,1))
plot!(sol; idxs=[sys.stopper.volume.v1.port.p/1e5], xlims=(0.5,1))
plot!(sol; idxs=[sys.stopper.volume.v2.port.p/1e5], xlims=(0.5,1))
plot!(sol; idxs=[sys.stopper.volume.v3.port.p/1e5], xlims=(0.5,1))
plot!(sol; idxs=[sys.stopper.volume.v4.port.p/1e5], xlims=(0.5,1))
plot!(sol; idxs=[sys.stopper.volume.v5.port.p/1e5], xlims=(0.5,1))





using CairoMakie
using Makie.GeometryBasics

function plot_train(sol, sys, tx=0.0)
    tm = Observable(tx)
    idx = Dict(reverse.(enumerate(states(sys))))

    fig = Figure()
    a = Axis(fig[1,1], aspect=DataAspect(), )
    hidedecorations!(a)
    slt(t,idxs) = sol(t; idxs) # =[sys.car1.x_a, sys.car1.x_b, sys.car2.x_a, sys.car2.x_b, sys.engine.x_a, sys.engine.x_b, sys.stopper.spring.x+sys.engine.barrier])

    c1xa(t) = slt(t, sys.car1.x_a)
    c1xb(t) = slt(t, sys.car1.x_b)

    flb1(t) = (c1xa(t), 0)
    frb1(t) = (c1xb(t), 0)
    frt1(t) = (c1xb(t), 0.25)
    delta1(t) = c1xb(t)-c1xa(t)
    flt1a(t) = (c1xb(t)-delta1(t)*0.75, 0.25)
    flt1b(t) = (c1xa(t), 0.1)

    ply1(t) = Polygon(Point2f[flb1(t), frb1(t), frt1(t), flt1a(t), flt1b(t)])
    poly!(a, lift(ply1, tm), color=:blue)

    c2xa(t) = slt(t, sys.car2.x_a)
    c2xb(t) = slt(t, sys.car2.x_b)

    flb2(t) = (c2xa(t), 0)
    frb2(t) = (c2xb(t), 0)
    frt2(t) = (c2xb(t), 0.25)
    flt2(t) = (c2xa(t), 0.25)

    ply2(t) = Polygon(Point2f[flb2(t), frb2(t), frt2(t), flt2(t)])
    poly!(a, lift(ply2, tm), color=:blue)


    c3xa(t) = slt(t, sys.car3.x_a)
    c3xb(t) = slt(t, sys.car3.x_b)

    flb3(t) = (c3xa(t), 0)
    frb3(t) = (c3xb(t), 0)
    frt3(t) = (c3xb(t), 0.25)
    flt3(t) = (c3xa(t), 0.25)

    ply3(t) = Polygon(Point2f[flb3(t), frb3(t), frt3(t), flt3(t)])
    poly!(a, lift(ply3, tm), color=:blue)



    exa(t) = slt(t, sys.engine.x_a)
    exb(t) = slt(t, sys.engine.x_b)

    flbe(t) = (exa(t), 0)
    frbe(t) = (exb(t), 0)
    deltae(t) = exb(t)-exa(t)
    frtea(t) = (exb(t), 0.1)
    frteb(t) = (exa(t)+deltae(t)*0.75, 0.25)
    flte(t) = (exa(t), 0.25)

    plye(t) = Polygon(Point2f[flbe(t), frbe(t), frtea(t), frteb(t), flte(t)])
    poly!(a, lift(plye, tm), color=:blue)


    br(t) = slt(t, (5-sys.stopper.volume.x)+10)
    brp(t) = Polygon(Point2f[
                            (br(t),0)
                            (br(t)+0.1,0)
                            (br(t)+0.1,0.1)
                            (br(t)+0.5,0.1)
                            (br(t)+0.5,0.2)
                            (br(t)+0.1,0.2)
                            (br(t)+0.1,0.3)
                            (br(t),0.3)
                            ])
    
    poly!(a, lift(brp, tm), color=(:red, 0.5))


    v1(t) = slt(t, 10-sys.stopper.volume.v1.x)
    v1c(t) = (:red, 0.1+slt(t, sys.stopper.volume.v1.port.p/1e5)/20 )
    v1p(t) = Polygon(Point2f[
                            (v1(t),0)
                            (11,0)
                            (11,0.3)
                            (v1(t),0.3)
                            ])
    poly!(a, lift(v1p, tm); color=lift(v1c, tm)) #, colormap = :jet, colorrange = (1, 20))


    v2(t) = slt(t, 10+1-sys.stopper.volume.v2.x)
    v2c(t) = (:red, 0.1+slt(t, sys.stopper.volume.v2.port.p/1e5)/20 )
    v2p(t) = Polygon(Point2f[
                            (v2(t),0)
                            (12,0)
                            (12,0.3)
                            (v2(t),0.3)
                            ])
    poly!(a, lift(v2p, tm); color=lift(v2c, tm)) #, colormap = :jet, colorrange = (1, 20))


    v3(t) = slt(t, 10+2-sys.stopper.volume.v3.x)
    v3c(t) = (:red, 0.1+slt(t, sys.stopper.volume.v3.port.p/1e5)/20 )
    v3p(t) = Polygon(Point2f[
                            (v3(t),0)
                            (13,0)
                            (13,0.3)
                            (v3(t),0.3)
                            ])
    poly!(a, lift(v3p, tm); color=lift(v3c, tm)) #, colormap = :jet, colorrange = (1, 20))

    v4(t) = slt(t, 10+3-sys.stopper.volume.v4.x)
    v4c(t) = (:red, 0.1+slt(t, sys.stopper.volume.v4.port.p/1e5)/20 )
    v4p(t) = Polygon(Point2f[
                            (v4(t),0)
                            (14,0)
                            (14,0.3)
                            (v4(t),0.3)
                            ])
    poly!(a, lift(v4p, tm); color=lift(v4c, tm)) #, colormap = :jet, colorrange = (1, 20))

    v5(t) = slt(t, 10+4-sys.stopper.volume.v5.x)
    v5c(t) = (:red, 0.1+slt(t, sys.stopper.volume.v5.port.p/1e5)/20 )
    v5p(t) = Polygon(Point2f[
                            (v5(t),0)
                            (15,0)
                            (15,0.3)
                            (v5(t),0.3)
                            ])
    poly!(a, lift(v5p, tm); color=lift(v5c, tm)) #, colormap = :jet, colorrange = (1, 20))

    # track
    lines!(a, [0,15], [-0.05,-0.05]; linewidth=1, color=:black)


    CairoMakie.ylims!(a, -5, 5)
    CairoMakie.xlims!(a, 0, 15)
    fig, tm
end

fig,tm = plot_train(sol, sys); fig
fig, = plot_train(sol, sys, 0.76); fig

function record_train(fig, tm; start, stop, step, framerate, filename)

    timestamps = range(start, stop; step)

    record(fig, filename, timestamps;
            framerate) do t
        tm[] = t
    end

    nothing
end

fig, tm = plot_train(sol, sys)
record_train(fig, tm; start=0.7125, stop=0.76, step=2e-4, framerate=30, filename="slow_motion.mp4")
record_train(fig, tm; start=0.0, stop=5, step=1/30, framerate=30, filename="run_away_train.mp4")

plot(sol; idxs=sys.stopper.volume.moving_volume.port.p/1e5, xlims=(0.7, 0.9))
plot!(sol; idxs=sys.stopper.volume.v1.port.p/1e5, xlims=(0.7, 0.9))
plot!(sol; idxs=sys.stopper.volume.v5.port.p/1e5, xlims=(0.7, 0.9))

plot(sol; idxs=sys.stopper.volume.moving_volume.dx, xlims=(0.7, 0.9))
plot!(sol; idxs=sys.stopper.volume.v1.dx, xlims=(0.7, 0.9))
plot!(sol; idxs=sys.stopper.volume.v5.dx, xlims=(0.7, 0.9))