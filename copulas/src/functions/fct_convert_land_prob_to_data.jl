function convert_land_prob_to_data(data, prob_scenarios, scenario_year, scenario_month, scenario_day, scenario_hour)
    """
    convert_land_prob_to_data:

    # Arguments
    - data::DataFrame
    - prob_scenarios::Matrix{Float64}
    - scenario_year::Int64
    - scenario_month::Int64
    - scenario_day::Int64
    - scenario_hour::Int64
    - 
    # Returns
    - 
    """
    # Create timestamp vectors
    #scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, 0)
    #scenario_timestamp_end = DateTime(scenario_year, scenario_month, scenario_day, 23)
    scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour);
    scenario_timestamp_end = scenario_timestamp_begin + Hour(23);
    
    #= dt returns 48 rows. Every odd-numbered row is a forecast from 24h
    to 48h ahead while every even-numbered row is a forecast up until 23h
    ahead. =#
    dt = filter(:forecast_time => x -> scenario_timestamp_begin <= x <= scenario_timestamp_end, data)

    up_to_24h = dt[2:2:48, :]
    after_24h = dt[1:2:48, :]

    indexes = startswith.(names(dt), "p_")
    scen_data = Matrix{Float64}(undef, number_of_scenarios, scenario_length)

    for i in 1:1:scenario_length÷2
        scen_data[:, i] = quantile(up_to_24h[i, indexes], prob_scenarios[:, i])
        scen_data[:, i+24] = quantile(after_24h[i, indexes], prob_scenarios[:, i])
    end

    return(scen_data)
end