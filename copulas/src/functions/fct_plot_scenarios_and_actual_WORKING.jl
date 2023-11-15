function plot_scenarios_and_actual_w(hist, scenarios, actuals, type, scenario_year, scenario_month, scenario_day, scenario_hour)
    """
    """
    # Set DateTime to fetch data.
    # Historical data are shown from midnight onwards. This means that 
    # if the scenarios start at any given hour H on a day D, the historical
    # line will start at midnight on day D regardless of H. 
    # scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour)
    # scenario_timestamp_end = scenario_timestamp_begin + Hour(23)
    historical_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day);
    scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour);
    timestamp_end = scenario_timestamp_begin + Hour(47);

    scenario_range = collect(scenario_timestamp_begin:Hour(1):timestamp_end);
    historical_range = collect(historical_timestamp_begin:Hour(1):timestamp_end);

    #= dt returns 48 rows. Every odd-numbered row is a forecast from 24h
    to 48h ahead while every even-numbered row is a forecast up until 23h
    ahead. =#
    # Fetching forecasts
    dt = filter(:forecast_time => x -> historical_timestamp_begin <= x <= timestamp_end, actuals);
    dt_onedayahead = dt[dt.ahead_factor .== "one",:];
    dt_historical = dt_onedayahead[:, :avg_actual];
    
    hist2 = filter(:time_index => x -> historical_timestamp_begin <= x <= timestamp_end, hist);
    hist2 = hist2[minute.(hist2.time_index) .== 0,:];
    # Parameters for the plot
    scenarios = Matrix(scenarios);
    
    title_hour = Dates.format(scenario_timestamp_begin, "Ip");
    title = monthname(scenario_month) * " " * string(scenario_day) * ", " * string(scenario_year) * " @ " *title_hour* " | " * string(type);
    upper_lim_y = maximum([maximum(scenarios), maximum(dt_historical)]);
    
    # X-axis text
    t = historical_timestamp_begin:Hour(6):timestamp_end;
    time_ticks = Dates.format.(t,"Ip");

    # Initialize plot
    plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation =:left,
        legend=:outerbottom,
        legendcolumns=3,
        background_color_legend=nothing,
        size=(1000,700))

    # Add scenarios
    plot!(scenario_range, scenarios[1,:], lc=:lightgrey, label = "Scenarios")
    for i in 2:size(scenarios,1) 
        plot!(scenario_range,scenarios[i,:], lc=:lightgrey, label = "")
    end


    # Add historical data
    plot!(historical_range, dt_historical,  lw = 2, lc=:seagreen, label = "Historical hourly average")

    # Add historical data
    scatter!(historical_range, hist2[:, :values], mc=:seagreen, ms=1.75, markerstrokewidth = 0, label = "Historical")

    # Add scenarios average
    scen_means = transpose(mean(scenarios, dims=1))
    plot!(scenario_range, scen_means, lw = 2, lc=:violetred, label = "Scenarios average")

    plot!(xticks = (t, time_ticks));
    xlabel!("Hours ahead");
    ylabel!(string(type) * " [MW]");
    ylims!(0, 1.1*upper_lim_y);
    title!(title);
    #vline!([scenario_timestamp_begin], linecolor = "grey", linestyle=:dash);
    

    # Save plot
    plot_name = string(scenario_year)*"_"*string(scenario_month)*"_"*string(scenario_day)*"_"*string(type)*"_"*".png"
    savefig(here("plots", plot_name))
end