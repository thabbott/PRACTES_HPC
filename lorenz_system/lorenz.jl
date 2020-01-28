using Random

@doc raw"""
    struct Position

A position (`x`, `y`, `z`) in 3-space
"""
mutable struct Position
    x
    y
    z
end

@doc raw"""
    function advance_position!(r, σ, β, ρ, dt)

Advance a position `r` with time step `dt` according to

```math
\frac{dr_1}{dt} = \sigma(r_2 - r_1)
\frac{dr_2]{dt} = r_1 (\rho - r_3) - r_2
\frac{dr_3}{dt} = r_1 r_2 - \beta r_3
```

"""
function advance_position!(r, σ, β, ρ, dt)
    dr1 = σ*(r.y - r.x)
    dr2 = r.x*(ρ - r.z) - r.y
    dr3 = r.x*r.y - β*r.z
    r.x += dr1*dt
    r.y += dr2*dt
    r.z += dr3*dt
    return nothing
end

@doc raw"""
    struct LorenzEnsemble

Positions `r` and parameters `σ`, `β`, and `ρ` for an ensemble of
Lorenz systems

"""
struct LorenzEnsemble
    N
    r
    σ
    β
    ρ
end

@doc raw"""
    function LorenzSystem(N; x = 1.0, y = 1.0, z = 1.0, σx = 1.0,
                             σ = 10.0, β = 8.0/3.0, ρ = 28.0)

Construct a Lorenz system with `N` ensemble members. Default parameters 
and initial conditions can be modified with keyword arguments
"""
function LorenzEnsemble(N; x = 1.0, y = 1.0, z = 1.0, σx = 1.0,
                           σ = 10.0, β = 8.0/3.0, ρ = 28.0)
    
    # Set initial conditions
    rng = MersenneTwister(747)
    r = Array{Position,1}(undef, N)
    for i = 1:N
        r[i] = Position(x + σx * rand(rng), y, z)
    end

    # Return ensemble
    return LorenzEnsemble(N, r, σ, β, ρ)

end

@doc raw"""
    function advance_ensemble!(ensemble, dt)

Advance each position in an ensemble of Lorenz systems
"""
function advance_ensemble!(ensemble, dt)
    for pos in ensemble.r
        advance_position!(pos, ensemble.σ, ensemble.β, ensemble.ρ, dt)
    end
    return nothing
end

@doc raw"""
    function fill!(x, z, ensemble)

Fill pre-allocated arrays with `x` and `z` positions
"""
function fill!(x, z, ensemble)
    for i in 1:ensemble.N
        x[i] = ensemble.r[i].x
        z[i] = ensemble.r[i].z
    end
    return nothing
end
