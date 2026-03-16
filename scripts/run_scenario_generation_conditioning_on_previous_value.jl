#= Scenario Generation with Copulas 
.......................................................
Hugo S. de Araujo
Nov. 14th, 2022 | Mays Group | Cornell University
.......................................................

This script runs the scenario generation procedure to create S sequences or 
time series with length T. Using the approach in this script, a sequence s in 
S will have dependence on the previous s values in the same sequence s. 

Because of this, once the |S| = T, no new sequences will be generated, therefore,
requiring a new seed. For example, 

Given S = 5 and T = 5, the sequences will be generated as follows:

s1 = 7, 6, 12, 4, 9
s2 = 7, 5, 8, 2, 11
s3 = 7, 5, 6, 1, 5
s4 = 7, 5, 6, 9, 15
s5 = 7, 5, 6, 9, 2 

That is one batch of sequences.


############################################################################# =#

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
include(joinpath(pwd(), "src", "fct_compute_hourly_average_actuals.jl"));
include(joinpath(pwd(), "src", "fct_compute_landing_probability.jl"));
include(joinpath(pwd(), "src", "fct_convert_hours_2018.jl"));
include(joinpath(pwd(), "src", "fct_convert_ISO_standard.jl"));
include(joinpath(pwd(), "src", "fct_convert_land_prob_to_data.jl"));
include(joinpath(pwd(), "src", "fct_generate_probability_scenarios.jl"));
include(joinpath(pwd(), "src", "fct_generate_IDM_scenarios.jl"));
include(joinpath(pwd(), "src", "fct_getplots.jl"));
include(joinpath(pwd(), "src", "fct_plot_historical_landing.jl"));
include(joinpath(pwd(), "src", "fct_plot_historical_synthetic_autocorrelation.jl"));
include(joinpath(pwd(), "src", "fct_plot_correlogram_landing_probability.jl"));
include(joinpath(pwd(), "src", "fct_plot_scenarios_and_actual.jl"));
include(joinpath(pwd(), "src", "fct_read_h5_file.jl"));
include(joinpath(pwd(), "src", "fct_read_input_file.jl"));
include(joinpath(pwd(), "src", "fct_transform_landing_probability.jl"));
include(joinpath(pwd(), "src", "fct_write_percentiles.jl"));
include(joinpath(pwd(), "src", "fct_write_scenarios.jl"));
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
READ INPUT FILE
=======================================================================#
input_file_path = projectdir("copulas.txt")

data_type,
scenario_length,
number_of_scenarios,
number_of_sheets,
number_of_iterations,
number_of_iterations_IDM,
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
Variance of percentiles 
=======================================================================#
#= Brief graphical analysis for the variance of the percentile forecasts.
The idea for this came from a meeting with Jacob Mays on March. 30th,
2023.

THIS FUNCTION HAS TO BE REWRITTEN USING PLOTS.JL
=#
# print_graphical_analysis = true
# if print_graphical_analysis
#     getplots(load_data, "Load data", "load")
#     getplots(solar_data, "Solar data", "solar")
#     getplots(wind_data, "Wind data", "wind")
# end

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
# load_prob_scen = generate_probability_scenarios(
#     lp_load, scenario_length, number_of_scenarios;
#     intraday_market_scenarios=true, num_batches=10);
# solar_prob_scen = generate_probability_scenarios(
#     lp_solar, scenario_length, number_of_scenarios;
#     intraday_market_scenarios=true, num_batches=10);
# wind_prob_scen = generate_probability_scenarios(
#     lp_wind, scenario_length, number_of_scenarios;
#     intraday_market_scenarios=true, num_batches=10);

load_scenarios_4d, load_w_4d = generate_probability_scenarios_cube!(
    lp_load,
    scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
    save_path_scenarios_4d=joinpath(pwd(), "results", "load_scenarios_4d.jls"),
    save_path_W_4d=joinpath(pwd(), "results", "load_w_4d.jls"),
);

solar_scenarios_4d, solar_w_4d = generate_probability_scenarios_cube!(
    lp_solar,
    scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
    save_path_scenarios_4d=joinpath(pwd(), "results", "solar_scenarios_4d.jls"),
    save_path_W_4d=joinpath(pwd(), "results", "solar_w_4d.jls"),
);

wind_scenarios_4d, wind_w_4d = generate_probability_scenarios_cube!(
    lp_wind,
    scenario_length, number_of_scenarios, number_of_iterations, number_of_sheets;
    save_path_scenarios_4d=joinpath(pwd(), "results", "wind_scenarios_4d.jls"),
    save_path_W_4d=joinpath(pwd(), "results", "wind_w_4d.jls"),
);

#=======================================================================
CONVERT PROBABILITY SCENARIOS INTO DATA SCENARIOS
=======================================================================#
# load_scen = convert_land_prob_to_data(
#     load_data, load_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);
# solar_scen = convert_land_prob_to_data(
#     solar_data, solar_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);
# wind_scen = convert_land_prob_to_data(
#     wind_data, wind_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);

load_weather_avg_scenarios, load_weather_scenarios =
    convert_land_prob_cube_to_data(
        load_data, load_scenarios_4d,
        scenario_year, scenario_month, scenario_day, scenario_hour;
        save_path_weather_4d=joinpath(pwd(), "results", "load_weather_4d.jls"),
    );

solar_weather_avg_scenarios, solar_weather_scenarios =
    convert_land_prob_cube_to_data(
        solar_data, solar_scenarios_4d,
        scenario_year, scenario_month, scenario_day, scenario_hour;
        save_path_weather_4d=joinpath(pwd(), "results", "solar_weather_4d.jls"),
    );

wind_weather_avg_scenarios, wind_weather_scenarios =
    convert_land_prob_cube_to_data(
        wind_data, wind_scenarios_4d,
        scenario_year, scenario_month, scenario_day, scenario_hour;
        save_path_weather_4d=joinpath(pwd(), "results", "wind_weather_4d.jls"),
    );

#=======================================================================
WRITE SCENARIOS TO FILE
=======================================================================#
write_scenarios_to_file(
    load_weather_avg_scenarios,
    scenario_day, scenario_month, scenario_year,
    "load"
)

write_scenarios_to_file(
    solar_weather_avg_scenarios,
    scenario_day, scenario_month, scenario_year,
    "solar"
)

write_scenarios_to_file(
    wind_weather_avg_scenarios,
    scenario_day, scenario_month, scenario_year,
    "wind"
)

#=======================================================================
SIMULATE SCENARIOS FOR PENALTY PRICES COMPUTATION
=======================================================================#
idm_seed = 29031990
if !isempty(intraday_hours)
    # Pre-allocate 5D arrays for weather data (iterations, IDM_iterations, sheets, scenarios, timesteps)
    # and 4D arrays for averages (iterations, IDM_iterations, 48, 48)
    load_weather_5d = Dict()
    solar_weather_5d = Dict()
    wind_weather_5d = Dict()

    load_avg_4d = Dict()
    solar_avg_4d = Dict()
    wind_avg_4d = Dict()

    for hour in intraday_hours
        # Allocate 5D: (iterations, IDM_iterations, sheets, scenarios, timesteps)
        load_weather_5d[hour] = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, number_of_sheets, number_of_scenarios, scenario_length)
        solar_weather_5d[hour] = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, number_of_sheets, number_of_scenarios, scenario_length)
        wind_weather_5d[hour] = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, number_of_sheets, number_of_scenarios, scenario_length)

        # Allocate 4D: (iterations, IDM_iterations, timesteps, timesteps)
        load_avg_4d[hour] = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, scenario_length, scenario_length)
        solar_avg_4d[hour] = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, scenario_length, scenario_length)
        wind_avg_4d[hour] = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, scenario_length, scenario_length)
    end

    # Loop over each iteration and hour, filling the arrays
    for iter in 1:number_of_iterations
        println("Generating IDM scenarios for Iteration $(iter)")

        for hour in intraday_hours
            println("  Processing intraday hour $(hour)")

            # Generate the probability scenarios ...................................
            # Don't save intermediate probability files - we'll save consolidated weather arrays at the end
            idm_load_scenarios, _ = generate_probability_IDM_scenarios_cube!(
                hour, lp_load,
                joinpath(pwd(), "results", "load_w_4d.jls"),
                scenario_length, number_of_scenarios, number_of_iterations_IDM, number_of_sheets;
                iteration_index=iter,
                seed=idm_seed
            )

            idm_solar_scenarios, _ = generate_probability_IDM_scenarios_cube!(
                hour, lp_solar,
                joinpath(pwd(), "results", "solar_w_4d.jls"),
                scenario_length, number_of_scenarios, number_of_iterations_IDM, number_of_sheets;
                iteration_index=iter,
                seed=idm_seed
            )

            idm_wind_scenarios, _ = generate_probability_IDM_scenarios_cube!(
                hour, lp_wind,
                joinpath(pwd(), "results", "wind_w_4d.jls"),
                scenario_length, number_of_scenarios, number_of_iterations_IDM, number_of_sheets;
                iteration_index=iter,
                seed=idm_seed
            )

            # Transform the probability scenarios into data scenarios ..............
            # Don't save individual files - we'll save consolidated arrays at the end
            load_weather_IDM_avg_scenarios, load_weather_scenarios =
                convert_land_prob_cube_to_data(
                    load_data, idm_load_scenarios,
                    scenario_year, scenario_month, scenario_day, scenario_hour
                )

            solar_weather_IDM_avg_scenarios, solar_weather_scenarios =
                convert_land_prob_cube_to_data(
                    solar_data, idm_solar_scenarios,
                    scenario_year, scenario_month, scenario_day, scenario_hour
                )

            wind_weather_IDM_avg_scenarios, wind_weather_scenarios =
                convert_land_prob_cube_to_data(
                    wind_data, idm_wind_scenarios,
                    scenario_year, scenario_month, scenario_day, scenario_hour
                )

            # Store in the 5D and 4D arrays
            load_weather_5d[hour][iter, :, :, :, :] = load_weather_scenarios
            solar_weather_5d[hour][iter, :, :, :, :] = solar_weather_scenarios
            wind_weather_5d[hour][iter, :, :, :, :] = wind_weather_scenarios

            load_avg_4d[hour][iter, :, :, :] = load_weather_IDM_avg_scenarios
            solar_avg_4d[hour][iter, :, :, :] = solar_weather_IDM_avg_scenarios
            wind_avg_4d[hour][iter, :, :, :] = wind_weather_IDM_avg_scenarios
        end
    end

    # After all loops, save consolidated 5D and 4D arrays (one file per intraday hour)
    for hour in intraday_hours
        println("Saving consolidated arrays for intraday hour $(hour)")

        # Save 5D weather arrays
        serialize(joinpath(pwd(), "results", "load_IDM_weather_5d_hour_$(hour).jls"), load_weather_5d[hour])
        serialize(joinpath(pwd(), "results", "solar_IDM_weather_5d_hour_$(hour).jls"), solar_weather_5d[hour])
        serialize(joinpath(pwd(), "results", "wind_IDM_weather_5d_hour_$(hour).jls"), wind_weather_5d[hour])

        # Save 4D average arrays
        serialize(joinpath(pwd(), "results", "load_IDM_avg_4d_hour_$(hour).jls"), load_avg_4d[hour])
        serialize(joinpath(pwd(), "results", "solar_IDM_avg_4d_hour_$(hour).jls"), solar_avg_4d[hour])
        serialize(joinpath(pwd(), "results", "wind_IDM_avg_4d_hour_$(hour).jls"), wind_avg_4d[hour])
    end

    # Zip 5D files by type and delete originals
    GC.gc()  # release file handles from serialize() calls before zipping
    let results_dir = joinpath(pwd(), "results")
        for type in ["load", "solar", "wind"]
            files = filter(f -> startswith(f, "$(type)_IDM_weather_5d") && endswith(f, ".jls"),
                           readdir(results_dir))
            files_full = joinpath.(results_dir, files)
            if !isempty(files_full)
                zip_path = joinpath(results_dir, "$(type)_IDM_weather_5d.zip")
                println("Zipping $(length(files_full)) files into $(zip_path)")
                file_list = join(["\"$f\"" for f in files_full], ",")
                ok = success(`powershell -Command "Compress-Archive -Path $file_list -DestinationPath \"$zip_path\" -Force"`)
                if ok && isfile(zip_path) && filesize(zip_path) > 0
                    for f in files_full
                        rm(f)
                    end
                    println("Deleted original 5D files for type: $(type)")
                else
                    println("WARNING: zip failed for type $(type) — original files NOT deleted")
                end
            end
        end
    end
end
