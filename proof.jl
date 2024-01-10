using ModelingToolkit

@parameters t
D = Differential(t)

@connector MechanicalPort begin
    v(t)
    f(t), [connect = Flow]
end

@mtkmodel ConstantForce begin
    @parameters begin
        f
    end
    @components begin
        port = MechanicalPort()
    end
    @equations begin
        # connectors
        port.f ~ f  
    end
end

@mtkmodel Mass begin
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
        port.f ~ -f
        
        # physics
        f ~ m*D(v)
    end
end

@mtkmodel System begin
    @components begin
        mass = Mass()
        force = ConstantForce()
    end
    @equations begin
        connect(mass.port, force.port)
    end
end
@mtkbuild sys = System()

full_equations(sys)