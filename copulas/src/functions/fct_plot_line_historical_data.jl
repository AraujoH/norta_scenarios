function plot_historical_data()
    """
    plot_all_historical_data(): create line plot for actual values of 
        load, solar, and wind. Stack plots into one column.

    # Arguments: path where to save image. 
    - 
    # Returns: .png file save to given path.
    - 
    """
    p1 = plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation=:left,
        legend=false,
        legendcolumns=3,
        labelfontsize=7,
        background_color_legend=nothing)
    plot!(load_actuals[:, :time_index], load_actuals[:, :values] / 1000, lc=:seagreen, dpi=300)
    ylabel!("Load profile [GW]")

    p2 = plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation=:left,
        legend=false,
        legendcolumns=3,
        labelfontsize=7,
        background_color_legend=nothing)
    plot!(solar_actuals[:, :time_index], solar_actuals[:, :values] / 1000, lc=:darkorange2, dpi=300)
    ylabel!("Solar profile [GW]")

    p3 = plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation=:left,
        legend=false,
        legendcolumns=3,
        labelfontsize=7,
        background_color_legend=nothing)
    plot!(wind_actuals[:, :time_index], wind_actuals[:, :values] / 1000, lc=:brown2, dpi=300)
    ylabel!("Wind profile [GW]")

    l = @layout [a; b; c]
    plot(p1, p2, p3, layout=l)
    savefig(joinpath(path,"all_historical_data.png"))
end