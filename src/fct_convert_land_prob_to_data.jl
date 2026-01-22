# function convert_land_prob_to_data_w(data, prob_scenarios, scenario_year, scenario_month, scenario_day, scenario_hour)
#     """
#     convert_land_prob_to_data:

#     # Arguments
#     - data::DataFrame
#     - prob_scenarios::Matrix{Float64}
#     - scenario_year::Int64
#     - scenario_month::Int64
#     - scenario_day::Int64
#     - scenario_hour::Int64
#     - 
#     # Returns
#     - 
#     """
#     # Create timestamp vectors
#     #scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, 0)
#     #scenario_timestamp_end = DateTime(scenario_year, scenario_month, scenario_day, 23)
#     scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour);
#     scenario_timestamp_end = scenario_timestamp_begin + Hour(47);

#     #= dt returns 48 rows. Every odd-numbered row is a forecast from 24h
#     to 48h ahead while every even-numbered row is a forecast up until 23h
#     ahead. =#
#     #= WORKING
#     FOR A GIVEN DATE RANGE I WANT SCENARIOS FOR, I WILL SELECT THE 
#     MARGINAL DISTRIBUTIONS THAT ARE CLOSER IN TIME TO THE FORECAST TIME
#     I WANT. MY THINKING IS THAT IF I HAVE TWO OPTIONS: (A) FORECASTS
#     ISSUED AT TIME X FOR A FORECAST TIME F, AND (B) FORECASTS ISSUED 
#     AT A TIME Y FOR A FORECAST TIME F, I WILL PICK THE MARGINAL DISTRIBUTION
#     ASSOCIATED WITH THE FORECAST ISSUED AT THE CLOSEST TIME TO MY F. 
#     IF F IS MAR. 29TH 6AM, X MAR.28 6PM, AND Y IS MAR. 29 1AM, I WILL
#     PICK THE MARGINAL DISTRIBUTIONS ASSOCIATED WITH Y
#     =#
#     dt = filter(:forecast_time => x -> scenario_timestamp_begin <= x <= scenario_timestamp_end, data);
#     dt = dt[dt.ahead_factor .== "one",:];

#     marg_distributions = select(dt, r"p_");
#     scen_data = Matrix{}(undef, number_of_scenarios, scenario_length);
#     for scen in axes(prob_scenarios, 1)
#         for j in axes(prob_scenarios, 2)
#             scen_data[scen, j] = quantile(marg_distributions[j,:], prob_scenarios[scen,j]);           
#         end
#     end

#     return(scen_data)
#     #= Old code 
#     up_to_24h = dt[2:2:48, :]
#     after_24h = dt[1:2:48, :]

#     indexes = startswith.(names(dt), "p_")
#     scen_data = Matrix{Float64}(undef, number_of_scenarios, scenario_length)

#     for i in 1:1:scenario_length÷2
#         scen_data[:, i] = quantile(up_to_24h[i, indexes], prob_scenarios[:, i])
#         scen_data[:, i+24] = quantile(after_24h[i, indexes], prob_scenarios[:, i])
#     end

#     return(scen_data)
#     --- Old code =#
# end

function convert_land_prob_to_data(
    data::DataFrame,
    prob_scenarios::AbstractMatrix{<:Real},
    scenario_year::Int,
    scenario_month::Int,
    scenario_day::Int,
    scenario_hour::Int
)
    # Build the time window for the horizon
    scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour)
    scenario_timestamp_end = scenario_timestamp_begin + Hour(size(prob_scenarios, 2) - 1)

    # Select the 48 relevant forecast times, then filter to the ahead group you want
    dt = filter(:forecast_time => t -> scenario_timestamp_begin <= t <= scenario_timestamp_end, data)
    dt = dt[dt.ahead_factor.=="one", :]

    # Keep only the marginal distribution columns (p_1, ..., p_T).
    # Each row corresponds to one forecast lead time, and each column
    # contains the empirical marginal samples for that lead time.
    marg = select(dt, r"p_")

    # Number of scenarios (rows) and time periods (columns)
    # implied by the probability scenarios matrix
    ns = size(prob_scenarios, 1)
    T = size(prob_scenarios, 2)

    # Consistency check:
    # For each time period j = 1,...,T, we need exactly one row in `marg`
    # defining the empirical marginal distribution used to invert
    # prob_scenarios[:, j]. If this is not true, the time alignment between
    # probability scenarios and marginal distributions is broken.
    if size(marg, 1) != T
        error(
            "Mismatch between marginal distributions and probability scenarios: " *
            "marg has $(size(marg, 1)) rows but prob_scenarios has $T columns."
        )
    end

    # Allocate output using prob_scenarios size
    scen_data = Matrix{Float64}(undef, ns, T)

    # For each time j, take the empirical distribution in row j of marg
    # and invert the probabilities for all scenarios
    for j in 1:T
        rowvals = collect(marg[j, :])  # Vector of p values for that time step
        for scen in 1:ns
            scen_data[scen, j] = quantile(rowvals, prob_scenarios[scen, j])
        end
    end

    return scen_data
end


function convert_land_prob_cube_to_data(
    data::DataFrame,
    scenarios_4d::AbstractArray{<:Real,4},   # scenarios_4d[i, s, scen, t] in (0,1)
    scenario_year::Int,
    scenario_month::Int,
    scenario_day::Int,
    scenario_hour::Int;
    save_path_weather_4d::AbstractString
)
    # Horizon window
    scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour)
    T = size(scenarios_4d, 4)
    scenario_timestamp_end = scenario_timestamp_begin + Hour(T - 1)

    # Select relevant times and ahead group
    dt = filter(:forecast_time => t -> scenario_timestamp_begin <= t <= scenario_timestamp_end, data)
    dt = dt[dt.ahead_factor .== "one", :]

    # Marginal samples per lead time (row j corresponds to time j)
    marg = select(dt, r"p_")

    if size(marg, 1) != T
        error(
            "Mismatch between marginal distributions and scenarios_4d: " *
            "marg has $(size(marg, 1)) rows but scenarios_4d has $T time steps."
        )
    end

    Nit     = size(scenarios_4d, 1)
    nsheets = size(scenarios_4d, 2)
    Nscen   = size(scenarios_4d, 3)

    # Output: weather_4d[i, s, scen, t]
    weather_4d = Array{Float64}(undef, Nit, nsheets, Nscen, T)

    # Convert probabilities to data via empirical quantiles, time by time
    for j in 1:T
        rowvals = collect(marg[j, :])  # empirical samples for time j

        for i in 1:Nit
            for s in 1:nsheets
                @inbounds for scen in 1:Nscen
                    p = scenarios_4d[i, s, scen, j]
                    weather_4d[i, s, scen, j] = quantile(rowvals, p)
                end
            end
        end
    end

    serialize(save_path_weather_4d, weather_4d)
    return weather_4d
end