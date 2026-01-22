"""    convert_land_prob_to_data(
        data::DataFrame,
        prob_scenarios::AbstractMatrix{<:Real},
        scenario_year::Int,
        scenario_month::Int,
        scenario_day::Int,
        scenario_hour::Int
    )

Convert probability scenarios to physical values using empirical marginal
distributions from historical forecast data.

Transforms probability scenarios (values in [0,1]) into actual physical quantities
(e.g., MW for load/wind/solar) by applying empirical quantile functions from
historical data for each time step in the forecast horizon.

# Arguments
- `data`: DataFrame containing historical forecast data with columns `forecast_time`,
  `ahead_factor`, and marginal distribution samples (`p_*` columns)
- `prob_scenarios`: Matrix `(number_of_scenarios, scenario_length)` with probability
  values in [0,1]
- `scenario_year`, `scenario_month`, `scenario_day`, `scenario_hour`: Starting timestamp
  for the forecast horizon

# Returns
- `scen_data`: Matrix `(number_of_scenarios, scenario_length)` containing physical
  values obtained by inverting probabilities through empirical quantiles

# Behavior
- Filters data to match the forecast horizon from the starting timestamp
- Selects one-ahead forecasts only (`ahead_factor == "one"`)
- For each time step `j`, uses the empirical distribution from `data` to map
  `prob_scenarios[:, j]` to physical values via quantile inversion
- Throws an error if time alignment between scenarios and marginals is inconsistent

# Example
```julia
scen_data = convert_land_prob_to_data(
    historical_df, prob_scenarios,
    2018, 7, 15, 12  # July 15, 2018 at noon
)
```
"""
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


"""    convert_land_prob_cube_to_data(
        data::DataFrame,
        scenarios_4d::AbstractArray{<:Real,4},
        scenario_year::Int,
        scenario_month::Int,
        scenario_day::Int,
        scenario_hour::Int;
        save_path_weather_4d::AbstractString
    )

Convert 4D probability scenarios to physical values and compute iteration-wise
averages across sheets.

Extends `convert_land_prob_to_data` to handle 4D scenario ensembles, transforming
probability values to physical quantities via empirical quantiles and computing
sheet-averaged scenarios for each iteration.

# Arguments
- `data`: DataFrame with historical forecast data including `forecast_time`,
  `ahead_factor`, and marginal samples (`p_*` columns)
- `scenarios_4d`: 4D array `[i, s, scen, t]` with probability values in [0,1]
  where `i` = iteration, `s` = sheet, `scen` = scenario, `t` = time
- `scenario_year`, `scenario_month`, `scenario_day`, `scenario_hour`: Forecast
  horizon starting timestamp
- `save_path_weather_4d`: Output path for serializing weather data arrays

# Returns
- `weather_average_3d`: Array `[i, scen, t]` with sheet-averaged physical values
  for each iteration
- `weather_4d`: Array `[i, s, scen, t]` with full physical value scenarios

# Behavior
- Converts all probability scenarios to physical values using empirical quantiles
- Computes mean across sheets for each iteration: `mean(weather_4d[i, :, :, :])`
- Serializes both `weather_4d` and `weather_average_3d` to disk
- Throws an error if time step count mismatches between scenarios and marginals

# Example
```julia
avg_weather, full_weather = convert_land_prob_cube_to_data(
    historical_df, scenarios_4d,
    2018, 7, 15, 12;
    save_path_weather_4d = "weather_data_4d.jls"
)
```
"""
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
    dt = dt[dt.ahead_factor.=="one", :]

    # Marginal samples per lead time (row j corresponds to time j)
    marg = select(dt, r"p_")

    if size(marg, 1) != T
        error(
            "Mismatch between marginal distributions and scenarios_4d: " *
            "marg has $(size(marg, 1)) rows but scenarios_4d has $T time steps."
        )
    end

    Nit = size(scenarios_4d, 1)
    nsheets = size(scenarios_4d, 2)
    Nscen = size(scenarios_4d, 3)

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

    # Average all scenarios for each iteration
    weather_average_3d = Array{Float64}(undef, Nit, Nscen, T)

    for i in 1:Nit
        weather_average_3d[i, :, :] .= 0.0
        for s in 1:nsheets
            weather_average_3d[i, :, :] .+= weather_4d[i, s, :, :]
        end
        weather_average_3d[i, :, :] ./= nsheets
    end


    serialize(save_path_weather_4d, weather_4d)
    serialize(replace(save_path_weather_4d, "weather_4d" => "average_weather_3d"), weather_average_3d)
    return weather_average_3d, weather_4d
end