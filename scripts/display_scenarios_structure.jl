"""
    display_scenarios_structure.jl

Visualizes the structure and distribution of load scenarios generated from 
the NORTA model, including Day-Ahead Market (DAM) scenarios and Intra-Day 
Market (IDM) scenarios at different hours.
"""

using CairoMakie

# Read data
dam_load = deserialize(joinpath(pwd(), "results", "load_average_weather_3d.jls"))
dam_load = dam_load[1, :, :]

idm_18_load = deserialize(joinpath(pwd(), "results", "load_IDM_average_weather_3d_intraday_hour_18.jls"))
idm_18_load = idm_18_load[1, :, :]


function plot_lines_for_scenarios(
    x::Matrix{Float64};
    title::String,
    xlabel::String="Operation hours",
    ylabel::String="Load (1000 x MW)",
    savepath=nothing
)

    nrows, T = size(x)

    fig = Figure(size=(900, 450))
    ax = Axis(fig[1, 1],
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )

    for i in 1:nrows
        lines!(ax, 1:T, x[i, :] ./ 1000, color=:gray)
    end

    # Spacing between x ticks
    ax.xticks = 0:5:T

    # Save if requested
    if savepath !== nothing
        mkpath(dirname(savepath))
        save(savepath, fig)
    end

    return fig
end


plot_lines_for_scenarios(
    dam_load;  title="DAM Load Scenarios",
    savepath=joinpath(pwd(), "results", "compare_dam_idm_scenario_tree_structure", "dam_load_scenarios.png")
)

plot_lines_for_scenarios(
    idm_18_load; title="IDM Load Scenarios at Hour 18",
    savepath=joinpath(pwd(), "results", "compare_dam_idm_scenario_tree_structure", "idm_load_scenarios.png")
)


