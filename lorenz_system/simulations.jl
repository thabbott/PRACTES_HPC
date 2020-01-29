using Random
using PyPlot
using IJulia

include("lorenz.jl")

@doc raw"""
    function animate_ensemble(; N = 1000)

Animate the evolution of an ensemble of `N` points with similar
(but not identical) initial conditions
"""
function animate_ensemble(; N = 1000)

    # Create the ensemble
    ensemble = LorenzEnsemble(N)

    # Pre-allocate arrays for plotting
    x = zeros(N)
    z = zeros(N)

    # Set up plot
    fig = figure()
    fill!(x, z, ensemble)
    scat = plot(x, z,
                color = "black", linestyle = "None",
                marker = "o", markersize = 0.5)[1]
    scat.set_zorder(20)
    ax = gca()
    xlim([-30,30])
    ylim([-10,60])
    ax.axis("off")
    fig.canvas.draw()
    sleep(0.01)
    
    # Integrate while plotting
    dt = 0.01
    try
        while true
            advance_ensemble!(ensemble, dt)
            fill!(x, z, ensemble)
            scat.set_xdata(x)
            scat.set_ydata(z)
            ax.draw_artist(scat)
            fig.canvas.blit(ax.bbox)
            IJulia.clear_output(true)
            display(fig)
            sleep(0.00001)
        end
    catch exception
    end

    # Plot and return final state
    return ensemble

end

@doc raw"""
    function benchmark(ensemble)

Run benchmarking targets. Note that macros to evaluate performance
(e.g. @time) should be applied to a call to this function and are
not included inside of it.
"""
function benchmark(ensemble)
    dt = 0.01
    for i = 1:1000
        advance_ensemble!(ensemble, dt)
    end
    return nothing
end
