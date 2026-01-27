#= Scenario Generation with Copulas 
.......................................................
Hugo S. de Araujo
Nov. 14th, 2022 | Mays Group | Cornell University
.......................................................
=#

"""
Analyzes the optimal number of sheets (nsheets) for copula-based scenario 
generation by evaluating stability metrics across different nsheets values.
"""

#=======================================================================
PROJECT SETUP
=======================================================================#
using Pkg
Pkg.activate("norta_scenarios")

# Import all required packages. 
begin
    using CairoMakie
    using CSV
    using DataFrames
    using Dates
    using DelimitedFiles
    using Distributions
    using HDF5
    using LinearAlgebra
    using LinearSolve
    using Printf
    using Random
    using RCall
    using Serialization
    using Statistics
    using StatsBase
    using Tables
    using TSFrames
    using TimeZones
end

# Include functions 
include(joinpath(pwd(), "src", "fct_bind_historical_forecast.jl"));
include(joinpath(pwd(), "src", "fct_check_stability_of_scenarios.jl"));
include(joinpath(pwd(), "src", "fct_compute_hourly_average_actuals.jl"));
include(joinpath(pwd(), "src", "fct_compute_landing_probability.jl"));
include(joinpath(pwd(), "src", "fct_convert_hours_2018.jl"));
include(joinpath(pwd(), "src", "fct_convert_ISO_standard.jl"));
include(joinpath(pwd(), "src", "fct_convert_land_prob_to_data.jl"));
include(joinpath(pwd(), "src", "fct_generate_probability_scenarios.jl"));
include(joinpath(pwd(), "src", "fct_read_h5_file.jl"));
include(joinpath(pwd(), "src", "fct_read_input_file.jl"));
include(joinpath(pwd(), "src", "fct_transform_landing_probability.jl"));
include(joinpath(pwd(), "src", "fct_write_percentiles.jl"));
include(joinpath(pwd(), "src", "fct_write_scenarios.jl"));
include(joinpath(pwd(), "src", "fct_generate_IDM_scenarios.jl"));
#=======================================================================
AUXILIARY FUNCTIONS
=======================================================================#

function projectdir(x::String)
    return (joinpath(pwd(), x))
end

function datadir(x::String)
    return (joinpath(pwd(), "data", x))
end

function plotsdir(x::String)
    return (joinpath(pwd(), "plots", x))
end

#=======================================================================
CREATE OUTPUT DIRECTORY FOR THE TASK
=======================================================================#
output_dir = mkpath(joinpath(pwd(), "results", "statibility_of_scenarios"))


#=======================================================================
READ INPUT FILE
=======================================================================#
input_file_path = joinpath(output_dir, "copulas.txt")

data_type,
scenario_length,
number_of_scenarios,
number_of_sheets,
number_of_iterations,
scenario_hour,
scenario_day,
scenario_month,
scenario_year,
intraday_hours,
read_locally,
historical_load,
forecast_load,
historical_solar,
forecast_da_solar,
forecast_2da_solar,
historical_wind,
forecastd_da_wind,
forecast_2da_wind,
write_percentile = read_input_file(input_file_path);

#=======================================================================
READ INPUT DATA: ARPA-E PERFORM PROJECT H5 FILES
=======================================================================#
# Function that reads the .h5 file and binds the time index and the actuals/fore-
# cast values into a single dataframe.

# Load data
load_actuals = read_h5_file(datadir("ercot_BA_load_actuals_2018.h5"), "load");
load_forecast = read_h5_file(datadir("ercot_BA_load_forecast_day_ahead_2018.h5"), "load", false);

# Solar data
solar_actuals = read_h5_file(datadir("ercot_BA_solar_actuals_Existing_2018.h5"), "solar");
solar_forecast_dayahead = read_h5_file(datadir("ercot_BA_solar_forecast_day_ahead_existing_2018.h5"), "solar", false);
solar_forecast_2dayahead = read_h5_file(datadir("ercot_BA_solar_forecast_2_day_ahead_existing_2018.h5"), "solar", false);

# Wind data
wind_actuals = read_h5_file(datadir("ercot_BA_wind_actuals_Existing_2018.h5"), "wind");
wind_forecast_dayahead = read_h5_file(datadir("ercot_BA_wind_forecast_day_ahead_existing_2018.h5"), "wind", false);
wind_forecast_2dayahead = read_h5_file(datadir("ercot_BA_wind_forecast_2_day_ahead_existing_2018.h5"), "wind", false);

#=======================================================================
Compute the hourly average for the actuals data
=======================================================================#
# Load
aux = compute_hourly_average_actuals(load_actuals);
load_actual_avg = DataFrame();
time_index = aux[:, :Index];
avg_actual = aux[:, :values_mean];
load_actual_avg[!, :time_index] = time_index;
load_actual_avg[!, :avg_actual] = avg_actual;

# Solar
aux = compute_hourly_average_actuals(solar_actuals);
time_index = aux[:, :Index];
avg_actual = aux[:, :values_mean];
solar_actual_avg = DataFrame();
solar_actual_avg[!, :time_index] = time_index;
solar_actual_avg[!, :avg_actual] = avg_actual;

# Wind
aux = compute_hourly_average_actuals(wind_actuals);
time_index = aux[:, :Index];
avg_actual = aux[:, :values_mean];
wind_actual_avg = DataFrame();
wind_actual_avg[!, :time_index] = time_index;
wind_actual_avg[!, :avg_actual] = avg_actual;

#=======================================================================
ADJUST THE TIME 
=======================================================================#
#= For the year of 2018, adjust the time to Texas' UTC (UTC-6 or UTC-5)
depending on daylight saving time =#

# Load data
load_actuals = convert_hours_2018(load_actuals);
load_actual_avg = convert_hours_2018(load_actual_avg);
load_forecast = convert_hours_2018(load_forecast, false);

# Solar data
solar_actuals = convert_hours_2018(solar_actuals);
solar_actual_avg = convert_hours_2018(solar_actual_avg);
solar_forecast_dayahead = convert_hours_2018(solar_forecast_dayahead, false);
solar_forecast_2dayahead = convert_hours_2018(solar_forecast_2dayahead, false);

# Wind data
wind_actuals = convert_hours_2018(wind_actuals);
wind_actual_avg = convert_hours_2018(wind_actual_avg);
wind_forecast_dayahead = convert_hours_2018(wind_forecast_dayahead, false);
wind_forecast_2dayahead = convert_hours_2018(wind_forecast_2dayahead, false);

#=======================================================================
BIND HOURLY HISTORICAL DATA WITH FORECAST DATA
========================================================================#
#= The binding is made by ("forecast_time" = "time_index"). This causes the 
average actual value to be duplicated, which is desired, given the # of rows
in the load_forecast is double that of load_actual. To distinguish a 
one-day-ahead forecast from a two-day-ahead forecast, the column "ahead_factor"
is introduced. Bind the day-ahead and two-day-ahead forecasts for wind and solar
to get all the forecast data into one object as it is for load forecast =#
load_data = bind_historical_forecast(true,
    load_actual_avg,
    load_forecast);

solar_data = bind_historical_forecast(false,
    solar_actual_avg,
    solar_forecast_dayahead,
    solar_forecast_2dayahead);

wind_data = bind_historical_forecast(false,
    wind_actual_avg,
    wind_forecast_dayahead,
    wind_forecast_2dayahead);


#=======================================================================
Write forecast percentile to files 
=======================================================================#
#write_percentile(load_data, "load", scenario_year, scenario_month, scenario_day, scenario_hour);
write_percentile = false
if write_percentile
    write_percentiles(load_data, "load", scenario_year, scenario_month, scenario_day, scenario_hour)
    write_percentiles(solar_data, "solar", scenario_year, scenario_month, scenario_day, scenario_hour)
    write_percentiles(wind_data, "wind", scenario_year, scenario_month, scenario_day, scenario_hour)
end

#=======================================================================
Landing probability
=======================================================================#
#= This section holds the calculation of the probability that the actual
value was equaled or superior than the forecast percentiles for a given
day. This is made possible by the estimation of an approximate CDF
computed on the forecast percentiles. Once estimated, this function is
used to find the "landing probability"; the prob. that the actual value
is equal or greater than a % percentage of the forecast percentile.
=#
#include(here("src", "functions", "fct_compute_landing_probability.jl"))
landing_probability_load = compute_landing_probability(load_data);
landing_probability_solar = compute_landing_probability(solar_data);
landing_probability_wind = compute_landing_probability(wind_data);

#=======================================================================
ADJUST LANDING PROBABILITY DATAFRAME
=======================================================================#
#= Analysis to address point J.Mays raised on Slack on Dec. 29,2022.
Sort the landing_probability dataframe by issue time. Then group the 
dataset by issue_time and count how many observations exist per 
issue_time. We're only interested in keeping the forecasts that share
the same issue_time 48 times since 48 is the length for the generation=#
lp_load = transform_landing_probability(landing_probability_load);
lp_solar = transform_landing_probability(landing_probability_solar);
lp_wind = transform_landing_probability(landing_probability_wind);

#=======================================================================
CORRELATION HEATMAP FOR THE LANDING PROBABILITIES
=======================================================================#
plot_correlogram = false;

if plot_correlogram
    plot_correlogram_landing_probability(lp_load, "Load")
    plot_correlogram_landing_probability(lp_solar, "Solar")
    plot_correlogram_landing_probability(lp_wind, "Wind")
end

#=======================================================================
SIMULATE INPUT THROUGH NORTA-LIKE APPROACH
=======================================================================#
load_scenarios_4d, load_w_4d = generate_probability_scenarios_cube!(
    lp_load,
    scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
    save_path_scenarios_4d=joinpath(pwd(), "results", "load_scenarios_4d.jls"),
    save_path_W_4d=joinpath(output_dir, "load_w_4d.jls"),
);

solar_scenarios_4d, solar_w_4d = generate_probability_scenarios_cube!(
    lp_solar,
    scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
    save_path_scenarios_4d=joinpath(pwd(), "results", "solar_scenarios_4d.jls"),
    save_path_W_4d=joinpath(output_dir, "solar_w_4d.jls"),
);

wind_scenarios_4d, wind_w_4d = generate_probability_scenarios_cube!(
    lp_wind,
    scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
    save_path_scenarios_4d=joinpath(pwd(), "results", "wind_scenarios_4d.jls"),
    save_path_W_4d=joinpath(output_dir, "wind_w_4d.jls"),
);

#=======================================================================
CONVERT PROBABILITY SCENARIOS INTO DATA SCENARIOS
=======================================================================#
load_weather_avg_scenarios, load_weather_scenarios =
    convert_land_prob_cube_to_data(
        load_data, load_scenarios_4d,
        scenario_year, scenario_month, scenario_day, scenario_hour;
        save_path_weather_4d=joinpath(output_dir, "load_weather_4d.jls"),
    );

solar_weather_avg_scenarios, solar_weather_scenarios =
    convert_land_prob_cube_to_data(
        solar_data, solar_scenarios_4d,
        scenario_year, scenario_month, scenario_day, scenario_hour;
        save_path_weather_4d=joinpath(output_dir, "solar_weather_4d.jls"),
    );

wind_weather_avg_scenarios, wind_weather_scenarios =
    convert_land_prob_cube_to_data(
        wind_data, wind_scenarios_4d,
        scenario_year, scenario_month, scenario_day, scenario_hour;
        save_path_weather_4d=joinpath(output_dir, "wind_weather_4d.jls"),
    );

#=======================================================================
SIMULATE SCENARIOS FOR PENALTY PRICES COMPUTATION
=======================================================================#
idm_seed = 29031990
if !isempty(intraday_hours)
    for hour in intraday_hours
        # Generate the probability scenarios ...................................
        idm_load_scenarios, _ = generate_probability_IDM_scenarios_cube!(
            hour, lp_load,
            joinpath(output_dir, "load_w_4d.jls"),
            scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
            seed=idm_seed,
            save_path_scenarios_4d=joinpath(output_dir, "load_scenarios_4d_idm_hour_$(hour).jls"),
            save_path_W_4d=joinpath(output_dir, "load_w_4d_idm_hour_$(hour).jls"),
        )

        idm_solar_scenarios, _ = generate_probability_IDM_scenarios_cube!(
            hour, lp_solar,
            joinpath(output_dir, "solar_w_4d.jls"),
            scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
            seed=idm_seed,
            save_path_scenarios_4d=joinpath(output_dir, "solar_scenarios_4d_idm_hour_$(hour).jls"),
            save_path_W_4d=joinpath(output_dir, "solar_w_4d_idm_hour_$(hour).jls"),
        )

        idm_wind_scenarios, _ = generate_probability_IDM_scenarios_cube!(
            hour, lp_wind,
            joinpath(output_dir, "wind_w_4d.jls"),
            scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
            seed=idm_seed,
            save_path_scenarios_4d=joinpath(output_dir, "load_scenarios_4d_idm_hour_$(hour).jls"),
            save_path_W_4d=joinpath(output_dir, "load_w_4d_idm_hour_$(hour).jls"),
        )

        # Transform the probability scenarios into data scenarios ..............
        load_weather_IDM_avg_scenarios, load_weather_scenarios =
            convert_land_prob_cube_to_data(
                load_data, idm_load_scenarios,
                scenario_year, scenario_month, scenario_day, scenario_hour;
                save_path_weather_4d=joinpath(output_dir, "load_IDM_weather_4d_intraday_hour_$(hour).jls"),
            )

        solar_weather_avg_scenarios, solar_weather_scenarios =
            convert_land_prob_cube_to_data(
                solar_data, idm_solar_scenarios,
                scenario_year, scenario_month, scenario_day, scenario_hour;
                save_path_weather_4d=joinpath(output_dir, "solar_IDM_weather_4d_intraday_hour_$(hour).jls"),
            )

        wind_weather_avg_scenarios, wind_weather_scenarios =
            convert_land_prob_cube_to_data(
                wind_data, idm_wind_scenarios,
                scenario_year, scenario_month, scenario_day, scenario_hour;
                save_path_weather_4d=joinpath(output_dir, "wind_IDM_weather_4d_intraday_hour_$(hour).jls"),
            )
    end
end

# ==============================================================================
# CHECK IDEAL NUMBER OF SHEETS BASED ON STABILITY OF SCENARIOS
# ==============================================================================
tol_vector = collect(0.001:0.005:0.1)
k_values = [5, 10, 15, 20]

load_weather_scenarios = deserialize(joinpath(output_dir, "load_weather_4d.jls"))

for k_consecutive in k_values
    stable_sheet_number = fill(NaN, length(tol_vector))

    best_x   = typemax(Int)
    best_tol = NaN
    best_y   = Float64[]   # will be replaced by y from the best tol
    best_idx = 0

    for (idx, tol) in enumerate(tol_vector)
        x, y = check_stability_of_scenarios_mean(
            load_weather_scenarios,
            1;
            tol_rel = tol,
            k_consecutive = k_consecutive
        )

        # x can be nothing if it never stabilizes
        stable_sheet_number[idx] = x === nothing ? NaN : x
        
        if x !== nothing && x < best_x
            best_x   = x
            best_tol = tol
            best_y   = y
            best_idx = idx
        end
    end

    # Plot ---------------------------------------------------------------------
    f = Figure(resolution = (800, 600))
    ax = Axis(
        f[1, 1],
        xlabel = "Number of scenario sheets",
        ylabel = "Tolerance",
        title  = "Number of Scenario Sheets vs. Tolerance (k = $(k_consecutive))"
    )

    barplot!(
        ax,
        tol_vector,
        stable_sheet_number;
        color = :turquoise4,
        width = 0.004,
        direction = :x
    )

    # x axis ticks every 5
    max_sheets = maximum(stable_sheet_number[.!isnan.(stable_sheet_number)]; init = 0.0)
    ax.xticks = 0:5:ceil(Int, max_sheets)

    # y axis: show all tolerance values
    ax.yticks = (tol_vector, [@sprintf("%.3f", v) for v in tol_vector])

    save(joinpath(output_dir, "sheets_vs_tolerance__k$(k_consecutive).png"), f)
end
