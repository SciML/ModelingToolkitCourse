using ModelingToolkit
using ModelingToolkitStandardLibrary.Mechanical.Translational
using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible

using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible: liquid_density



@component function Volume(;
    #initial conditions
    x,
    dx=0,
    p,
    drho=0,
    dm=0,

    #parameters
    area,
    direction=+1, name)
    pars = @parameters begin
        area = area
    end

    vars = @variables begin
        x(t) = x
        dx(t) = dx
        p(t) = p
        f(t) = p * area
        rho(t)
        drho(t) = drho
        dm(t) = dm
    end

    systems = @named begin
        port = HydraulicPort(; p_int=p)
        flange = MechanicalPort(; f, v=dx)
    end

    eqs = [
        # connectors
        port.p ~ p
        port.dm ~ dm
        flange.v * direction ~ dx
        flange.f * direction ~ -f

        # differentials
        D(x) ~ dx
        D(rho) ~ drho

        # physics
        rho ~ liquid_density(port, p)
        f ~ p * area
        dm ~ drho * x * area + rho * dx * area]

    ODESystem(eqs, t, vars, pars; name, systems, defaults=[rho => liquid_density(port)])
end