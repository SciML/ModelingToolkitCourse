# Lecture 2: Developing high-fidelity models of hydraulic systems
Why focus on hydraulics?  The answer is essentially hydraulic modelling is really hard (in numerical computing terms, hydraulic models are often refered to as "stiff" ODE's, which require more rigorous solvers from standard ODE's).  Solving the challenges of modeling hydraulics is applicable to the numerical modeling challenges of all other domains.  Let's first start with the concept of *compressibility*.  Often we think of a liquid as incompressible, imagine attempting to "squeeze" water, it can be done but takes some very high forces.  Therefore, if the model in question won't be solving a problem with high forces, it can be assumed incompressible.  However, most hydrulic industrial models will involve high forces, this is precisly the area where most hydraulic machines are used.  

## Compressibility

### Density

Density is simply mass over volume

```math
\rho = m/V
```

Given a volume and mass of liquid, if the volume were to change from ``V_0`` to ``V``, we know that the pressure would increase, and since the mass in this case was constant, the density will increase as well.

![volume change](../img/VolumeChange.svg)

The change in pressure for an isothermal compressible process is typically given as

```math
\Delta p = -\beta \frac{\Delta V}{V_0}
```

### Calculating Density as a Function of Pressure

Substituting ``\Delta p`` and ``\Delta V``

```math
p - p_0 = -\beta \frac{V - V_0}{V_0}
```

substituting ``V = m / \rho ``

```math
p - p_0 = -\beta (1 - \rho/\rho_0)  
``` 

Solving for ``\rho``

```math
\rho = \rho_0 (1 + (p - p_0)/\beta)
```

Taking a known ``\rho_0`` when ``p_0`` is 0 (at gage pressure), simplifies to

```math
\rho = \rho_0 (1 + p/\beta) 
```

### Change in Mass

Conservation of mass gives us

```math
m_{in} - m_{out} = m_s 
```

The stored mass of oil is simply

```math
m_s = \rho V 
```

Taking the derivative gives us the rate of mass change

```math
\dot{m}_{in} - \dot{m}_{out} = \frac{\delta (\rho V)}{\delta t} 
```

Here is where the standard hydraulic modeling gets the physics wrong...

Correct Derivation (1):  

```math
\frac{\delta (\rho V)}{\delta t} = \dot{\rho} V + \rho \dot{V} 
```

Standard Practice[^1] (2):  

```math
\color{red} \frac{\delta (\rho V)}{\delta t} = \dot{\rho} V + \rho_0 \dot{V}   
```

Given ``\dot{\rho} = \rho_0 (\dot{p} / \beta)``, and ``q = \dot{m}/\rho_0`` the above is often written as

```math
\color{red} q_{in} - q_{out} = (\dot{p} / \beta) V + \dot{V} 
```

[^1]: See [simscape hydraulic chamber](https://www.mathworks.com/help/simscape/ref/variablehydraulicchamber.html).  Note the deprication warning moving to isothermal liquid library which uses the correct derivation.

### Example
Problem Definition - Given

- ``ṁ_{in} = f(t)``
- ``m = 3,000 kg``
- ``A = 900 cm^2``
- ``\rho_0 = 876 kg/m^3``
- ``\beta = 1.2e9 Pa/m^3``
- ``g = 9.807 m/s^2``

![example](../img/Example.svg)

