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

    scaling(t,p) = (10+log((slt(t,p)/1e5+0.01)/100))/10

    v5(t) = slt(t, 10-sys.stopper.volume.v5.x)
    v5c(t) = (:red, scaling(t, sys.stopper.volume.v5.port.p) )
    v5p(t) = Polygon(Point2f[
                            (v5(t),0)
                            (11,0)
                            (11,0.3)
                            (v5(t),0.3)
                            ])
    poly!(a, lift(v5p, tm); color=lift(v5c, tm)) #, colormap = :jet, colorrange = (1, 20))


    v4(t) = slt(t, 10+1-sys.stopper.volume.v4.x)
    v4c(t) = (:red, scaling(t, sys.stopper.volume.v4.port.p) )
    v4p(t) = Polygon(Point2f[
                            (v4(t),0)
                            (12,0)
                            (12,0.3)
                            (v4(t),0.3)
                            ])
    poly!(a, lift(v4p, tm); color=lift(v4c, tm)) #, colormap = :jet, colorrange = (1, 20))


    v3(t) = slt(t, 10+2-sys.stopper.volume.v3.x)
    v3c(t) = (:red, scaling(t, sys.stopper.volume.v3.port.p) )
    v3p(t) = Polygon(Point2f[
                            (v3(t),0)
                            (13,0)
                            (13,0.3)
                            (v3(t),0.3)
                            ])
    poly!(a, lift(v3p, tm); color=lift(v3c, tm)) #, colormap = :jet, colorrange = (1, 20))

    v2(t) = slt(t, 10+3-sys.stopper.volume.v2.x)
    v2c(t) = (:red, scaling(t, sys.stopper.volume.v2.port.p) )
    v2p(t) = Polygon(Point2f[
                            (v2(t),0)
                            (14,0)
                            (14,0.3)
                            (v2(t),0.3)
                            ])
    poly!(a, lift(v2p, tm); color=lift(v2c, tm)) #, colormap = :jet, colorrange = (1, 20))

    v1(t) = slt(t, 10+4-sys.stopper.volume.v1.x)
    v1c(t) = (:red, scaling(t, sys.stopper.volume.v1.port.p) )
    v1p(t) = Polygon(Point2f[
                            (v1(t),0)
                            (15,0)
                            (15,0.3)
                            (v1(t),0.3)
                            ])
    poly!(a, lift(v1p, tm); color=lift(v1c, tm)) #, colormap = :jet, colorrange = (1, 20))

    # track
    lines!(a, [0,15], [-0.05,-0.05]; linewidth=1, color=:black)


    CairoMakie.ylims!(a, -5, 5)
    CairoMakie.xlims!(a, 0, 15)
    fig, tm
end

# fig,tm = plot_train(sol, sys); fig
# fig, = plot_train(sol, sys, 0.76); fig

function record_train(fig, tm; start, stop, step, framerate, filename)

    timestamps = range(start, stop; step)

    record(fig, filename, timestamps;
            framerate) do t
        tm[] = t
    end

    nothing
end

fig, tm = plot_train(sol, sys)
record_train(fig, tm; start=0.72, stop=0.729, step=1e-4, framerate=20, filename="run_away_train_slow.mp4")
record_train(fig, tm; start=0.0, stop=5, step=1/30, framerate=30, filename="run_away_train.mp4")

# plot(sol; idxs=sys.stopper.volume.moving_volume.port.p/1e5, xlims=(0.7, 0.9))
# plot!(sol; idxs=sys.stopper.volume.v1.port.p/1e5, xlims=(0.7, 0.9))
# plot!(sol; idxs=sys.stopper.volume.v5.port.p/1e5, xlims=(0.7, 0.9))

# plot(sol; idxs=sys.stopper.volume.moving_volume.dx, xlims=(0.7, 0.9))
# plot!(sol; idxs=sys.stopper.volume.v1.dx, xlims=(0.7, 0.9))
# plot!(sol; idxs=sys.stopper.volume.v5.dx, xlims=(0.7, 0.9))