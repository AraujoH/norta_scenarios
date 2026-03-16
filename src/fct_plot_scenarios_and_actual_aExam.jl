function plot_scenarios_and_actual_aExam(hist, scenarios, actuals, type, scenario_year, scenario_month, scenario_day, scenario_hour)
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
        framestyle=:grid,
        titlelocation =:left,
        legend=:outerbottom,
        legendcolumns=3,
        background_color_legend=nothing,
        foreground_color_legend = nothing,
        titlefontsize=10,
        guidefontsize=8,
        tickfontsize=8,
        legendfontsize=6,
        dpi = 300,
        bg = :transparent,)

    # Add scenarios
    plot!(scenario_range, scenarios[1,:], lc=:navajowhite2, label = "Scenarios")
    for i in 2:size(scenarios,1) 
        plot!(scenario_range,scenarios[i,:], lc=:navajowhite2, label = "")
    end


    # Add historical data
    plot!(historical_range, dt_historical,  lw = 2.5, lc=:blue3, label = "Historical hourly average")

    # Add historical data
    scatter!(historical_range, hist2[:, :values], mc=:blue3, ms=1.75, markerstrokewidth = 0, label = "Historical")

    # Add scenarios average
    scen_means = transpose(mean(scenarios, dims=1))
    plot!(scenario_range, scen_means, lw = 2.5, lc=:orangered3, label = "Scenarios average")

    plot!(xticks = (t, time_ticks));
    xlabel!("Hours ahead");
    ylabel!(string(type) * " [MW]");
    ylims!(0, 1.1*upper_lim_y);
    title!(title);
    #vline!([scenario_timestamp_begin], linecolor = "grey", linestyle=:dash);
    

    # Save plot
    # Save plot
    filename = string(scenario_year)*"_"*string(scenario_month)*"_"*string(scenario_day)*"_"*string(type)*".png"
    savefig(joinpath("C:\\Users\\hsd36\\OneDrive - Cornell University\\_PhD\\important_exams\\A\\plots", filename))
end