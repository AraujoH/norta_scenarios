#= IDM-only driver: regenerates the IDM penalty files from an existing base cube.
.......................................................
Hugo S. de Araujo
Oct. 6th, 2026 | Mays Group | Cornell University
.......................................................

What this script does (and does not do), relative to main.jl:
  - skips the base (day-ahead NORTA) step entirely
  - reads load_w_4d.jls, solar_w_4d.jls, wind_w_4d.jls from <base_folder>
  - runs the IDM loop of main.jl unchanged (seed = idm_seed + 1000*hour + iter
    when IDM SEED ITER HOUR = 1 in copulas.txt), writing the 4D averages
    {load,solar,wind}_IDM_avg_4d_hour_<h>.jls into <output_folder>
  - skips the 5D save (the full scenario cubes are never written)
All other parameters come from the committed copulas.txt. The scenario date
is overridden to Apr. 3 or Jul. 18 when the output folder basename starts with
"april" or "july"; otherwise the copulas.txt date is used.

usage: julia +1.12.1 --project=norta_scenarios scripts/run_idm_only.jl <base_folder> <output_folder>
(the root Project.toml environment fails to precompile CairoMakie; use norta_scenarios/)

Run with Julia 1.12.1 on 2026-10-06 with the following arguments, producing
results/<month>_10outerIter_1idmIter_50idmPath_seedIterHour:
  july : results/10outerIter_40idmIter_48scenarioPath_48idmScenarioPath
         results/july_10outerIter_1idmIter_50idmPath_seedIterHour
  april: results/april_10outerIter_50idmIter_48scenarioPath_48idmScenarioPath
         results/april_10outerIter_1idmIter_50idmPath_seedIterHour
############################################################################# =#
t_wall0 = time()
cd(dirname(@__DIR__))   # project root (parent of scripts/)
length(ARGS) == 2 || error("usage: julia +1.12.1 --project=norta_scenarios scripts/run_idm_only.jl <base_folder> <output_folder>")
base_dir = abspath(ARGS[1])
results_dir = mkpath(abspath(ARGS[2]))
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
function projectdir(x::String)
    return (joinpath(pwd(), x))
end

function datadir(x::String)
    return (joinpath(pwd(), "data", x))
end

function plotsdir(x::String)
    return (joinpath(pwd(), "plots", x))
end


#=== parameters: copulas.txt, with the scenario date overridden per month ===#
data_type, scenario_length, number_of_scenarios, number_of_scenarios_idm, number_of_sheets,
number_of_iterations, number_of_iterations_IDM, scenario_hour, scenario_day, scenario_month,
scenario_year, intraday_hours, read_locally, historical_load, forecast_load, historical_solar,
forecast_da_solar, forecast_2da_solar, historical_wind, forecastd_da_wind, forecast_2da_wind,
write_percentile, idm_seed_by_iter_hour = read_input_file(projectdir("copulas.txt"));

# month tag from the output folder name selects the scenario date (april: Apr. 3, july: Jul. 18)
month = lowercase(basename(results_dir))
if startswith(month, "july")
    scenario_day, scenario_month = 18, 7
elseif startswith(month, "april")
    scenario_day, scenario_month = 3, 4
end
@assert normpath(results_dir) != normpath(base_dir)
@assert (number_of_iterations, number_of_iterations_IDM, number_of_scenarios_idm, number_of_sheets, number_of_scenarios) == (10, 1, 50, 50, 48)
@assert intraday_hours == [0, 1, 6, 12, 18]
@assert idm_seed_by_iter_hour == 1
@assert scenario_hour == 0 && scenario_year == 2018
println("date=$(scenario_year)-$(scenario_month)-$(scenario_day) base=$base_dir out=$results_dir")
println("params: iters=$number_of_iterations idmIters=$number_of_iterations_IDM idmPaths=$number_of_scenarios_idm sheets=$number_of_sheets hours=$intraday_hours seedflag=$idm_seed_by_iter_hour")

t_data0 = time()
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

t_data = time() - t_data0
println("data+landing-probability setup: $(round(t_data, digits=1)) s")

#=== IDM step only (identical to the main.jl loop, 5D arrays skipped) ===#
t_idm0 = time()
idm_seed = 29031990
load_avg_4d = Dict(); solar_avg_4d = Dict(); wind_avg_4d = Dict()
for hour in intraday_hours
    load_avg_4d[hour]  = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, number_of_scenarios_idm, scenario_length)
    solar_avg_4d[hour] = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, number_of_scenarios_idm, scenario_length)
    wind_avg_4d[hour]  = Array{Float64}(undef, number_of_iterations, number_of_iterations_IDM, number_of_scenarios_idm, scenario_length)
end

for iter in 1:number_of_iterations
    println("Generating IDM scenarios for Iteration $(iter)")
    for hour in intraday_hours
        idm_seed_iter_hour = idm_seed_by_iter_hour == 1 ? idm_seed + 1000 * hour + iter : (hour == 0 ? idm_seed + iter : idm_seed)
        println("  hour $(hour) seed $(idm_seed_iter_hour)")
        idm_load_scenarios, _ = generate_probability_IDM_scenarios_cube!(
            hour, lp_load, joinpath(base_dir, "load_w_4d.jls"),
            scenario_length, number_of_scenarios_idm, number_of_iterations_IDM, number_of_sheets;
            iteration_index=iter, seed=idm_seed_iter_hour)
        idm_solar_scenarios, _ = generate_probability_IDM_scenarios_cube!(
            hour, lp_solar, joinpath(base_dir, "solar_w_4d.jls"),
            scenario_length, number_of_scenarios_idm, number_of_iterations_IDM, number_of_sheets;
            iteration_index=iter, seed=idm_seed_iter_hour)
        idm_wind_scenarios, _ = generate_probability_IDM_scenarios_cube!(
            hour, lp_wind, joinpath(base_dir, "wind_w_4d.jls"),
            scenario_length, number_of_scenarios_idm, number_of_iterations_IDM, number_of_sheets;
            iteration_index=iter, seed=idm_seed_iter_hour)

        load_avg, _  = convert_land_prob_cube_to_data(load_data,  idm_load_scenarios,  scenario_year, scenario_month, scenario_day, scenario_hour)
        solar_avg, _ = convert_land_prob_cube_to_data(solar_data, idm_solar_scenarios, scenario_year, scenario_month, scenario_day, scenario_hour)
        wind_avg, _  = convert_land_prob_cube_to_data(wind_data,  idm_wind_scenarios,  scenario_year, scenario_month, scenario_day, scenario_hour)

        load_avg_4d[hour][iter, :, :, :]  = load_avg
        solar_avg_4d[hour][iter, :, :, :] = solar_avg
        wind_avg_4d[hour][iter, :, :, :]  = wind_avg
    end
end

for hour in intraday_hours
    serialize(joinpath(results_dir, "load_IDM_avg_4d_hour_$(hour).jls"),  load_avg_4d[hour])
    serialize(joinpath(results_dir, "solar_IDM_avg_4d_hour_$(hour).jls"), solar_avg_4d[hour])
    serialize(joinpath(results_dir, "wind_IDM_avg_4d_hour_$(hour).jls"),  wind_avg_4d[hour])
end
t_idm = time() - t_idm0
println("IDM loop + save: $(round(t_idm, digits=1)) s")
println("WALL: total $(round(time() - t_wall0, digits=1)) s (package load $(round(t_data0 - t_wall0, digits=1)) s, data setup $(round(t_data, digits=1)) s, IDM $(round(t_idm, digits=1)) s)")
