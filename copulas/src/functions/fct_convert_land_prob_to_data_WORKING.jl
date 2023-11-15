function convert_land_prob_to_data_w(data, prob_scenarios, scenario_year, scenario_month, scenario_day, scenario_hour)
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
    scenario_timestamp_end = scenario_timestamp_begin + Hour(47);
    
    #= dt returns 48 rows. Every odd-numbered row is a forecast from 24h
    to 48h ahead while every even-numbered row is a forecast up until 23h
    ahead. =#
    #= WORKING
    FOR A GIVEN DATE RANGE I WANT SCENARIOS FOR, I WILL SELECT THE 
    MARGINAL DISTRIBUTIONS THAT ARE CLOSER IN TIME TO THE FORECAST TIME
    I WANT. MY THINKING IS THAT IF I HAVE TWO OPTIONS: (A) FORECASTS
    ISSUED AT TIME X FOR A FORECAST TIME F, AND (B) FORECASTS ISSUED 
    AT A TIME Y FOR A FORECAST TIME F, I WILL PICK THE MARGINAL DISTRIBUTION
    ASSOCIATED WITH THE FORECAST ISSUED AT THE CLOSEST TIME TO MY F. 
    IF F IS MAR. 29TH 6AM, X MAR.28 6PM, AND Y IS MAR. 29 1AM, I WILL
    PICK THE MARGINAL DISTRIBUTIONS ASSOCIATED WITH Y
    =#
    dt = filter(:forecast_time => x -> scenario_timestamp_begin <= x <= scenario_timestamp_end, data);
    dt = dt[dt.ahead_factor .== "one",:];

    marg_distributions = select(dt, r"p_");
    scen_data = Matrix{}(undef, number_of_scenarios, scenario_length);
    for scen in axes(prob_scenarios, 1)
        for j in axes(prob_scenarios, 2)
            scen_data[scen, j] = quantile(marg_distributions[j,:], prob_scenarios[scen,j]);           
        end
    end

    return(scen_data)
    #= Old code 
    up_to_24h = dt[2:2:48, :]
    after_24h = dt[1:2:48, :]

    indexes = startswith.(names(dt), "p_")
    scen_data = Matrix{Float64}(undef, number_of_scenarios, scenario_length)

    for i in 1:1:scenario_length÷2
        scen_data[:, i] = quantile(up_to_24h[i, indexes], prob_scenarios[:, i])
        scen_data[:, i+24] = quantile(after_24h[i, indexes], prob_scenarios[:, i])
    end

    return(scen_data)
    --- Old code =#
end