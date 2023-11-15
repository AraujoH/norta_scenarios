function plot_scenarios_and_actual(scenarios, actuals, type, scenario_year, scenario_month, scenario_day, scenario_hour)
    """
    """
    # Obtain the historical data
    # scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour)
    # scenario_timestamp_end = scenario_timestamp_begin + Hour(23)
    scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour);
    scenario_timestamp_end = scenario_timestamp_begin + Hour(47);

    #= dt returns 48 rows. Every odd-numbered row is a forecast from 24h
    to 48h ahead while every even-numbered row is a forecast up until 23h
    ahead. =#
    #dt = filter(:forecast_time => x -> scenario_timestamp_begin <= x <= scenario_timestamp_end, actuals);
    #dt_historical = dt[:, :avg_actual];

    dt = filter(:forecast_time => x -> scenario_timestamp_begin <= x <= scenario_timestamp_end, actuals);
    dt_onedayahead = dt[dt.ahead_factor .== "one",:];
    dt_historical = dt_onedayahead[:, :avg_actual];

    # Parameters for the plot
    scenarios = Matrix(scenarios);
    title = monthname(scenario_month) * " " * string(scenario_day) * ", " * string(scenario_year) * " | " * string(type);
    upper_lim_y = maximum([maximum(scenarios), maximum(dt_historical)]);
    
    # Initialize plot
    plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation =:left,
        legend=:outerbottom,
        legendcolumns=3,
        background_color_legend=nothing,
        size=(1000,700)) #Removes border from legend


    # Add scenarios
    plot!(scenarios[1,:], lc=:lightgrey, label = "Scenarios")
    for i in 2:size(scenarios,1) 
        plot!(scenarios[i,:], lc=:lightgrey, label = "")
    end

    # Add historical data
    plot!(dt_historical, lw = 2, lc=:seagreen, label = "Historical")

    # Add scenarios average
    scen_means = transpose(mean(scenarios, dims=1))
    plot!(scen_means, lw = 2, lc=:violetred, label = "Scenarios average")

    xlabel!("Hours ahead")
    ylabel!(string(type) * " [MW]")
    ylims!(0, 1.1*upper_lim_y)
    title!(title)

    # Save plot
    plot_name = string(scenario_year)*"_"*string(scenario_month)*"_"*string(scenario_day)*"_"*string(type)*"_"*".png"
    savefig(here("plots", plot_name))
end