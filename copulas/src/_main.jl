# Scenario Generation with Copulas 
# 
# Hugo S. de Araujo
# Nov. 14th, 2022 | Mays Group | Cornell University
################################################################################

#=======================================================================
PROJECT SETUP
=======================================================================#
using Pkg
Pkg.activate("copulas");
Pkg.instantiate();
# Import "here" function. Wrapper to allow easy path concatenation.
include(joinpath(@__DIR__, "functions", "fct_here.jl"))

# Import all required packages. 
begin
    using AWSS3
    using CSV
    using DataFrames
    using Dates
    using DelimitedFiles
    using Distributions
    using HDF5
    using JuliaFormatter
    using LaTeXStrings
    using LinearAlgebra
    using LinearSolve
    #using Measures
    using Random
    using RCall
    using Revise
    using Statistics
    using StatsBase
    #using StatsPlots
    using OhMyREPL
    using Plots
    #using PrettyTables
    using Tables
    using TSFrames
    using TimeZones
end

# Include functions 
#= functions_dirpath = joinpath(pwd(),"src", "functions");
function_paths = readdir(functions_dirpath, join=true);
function_index = occursin.(".jl", function_paths);
functions_only = function_paths[function_index];

for str in functions_only
    include(str)
end =#

include(here("src", "functions", "fct_bind_historical_forecast.jl"));
include(here("src", "functions", "fct_compute_hourly_average_actuals.jl"));
include(here("src", "functions", "fct_compute_landing_probability.jl"));
include(here("src", "functions", "fct_convert_hours_2018.jl"));
include(here("src", "functions", "fct_convert_ISO_standard.jl"));
include(here("src", "functions", "fct_convert_land_prob_to_data.jl"));
include(here("src", "functions", "fct_generate_probability_scenarios.jl"));
include(here("src", "functions", "fct_getplots.jl"));
include(here("src", "functions", "fct_plot_correlation_heatmap.jl"));
include(here("src", "functions", "fct_plot_historical_landing.jl"));
include(here("src", "functions", "fct_plot_historical_synthetic_autocorrelation.jl"));
include(here("src", "functions", "fct_plot_correlogram_landing_probability.jl"));
include(here("src", "functions", "fct_plot_scenarios_and_actual.jl"));
include(here("src", "functions", "fct_read_h5_file.jl"));
include(here("src", "functions", "fct_read_input_file.jl"));
include(here("src", "functions", "fct_transform_landing_probability.jl"));
include(here("src", "functions", "fct_write_percentiles.jl"));
#=======================================================================
READ INPUT FILE
=======================================================================#
input_file_path = here("src\\copulas.txt")

# XXX Needs to be updated to be a hardcoded instead of reading in a text file
data_type,
scenario_length,
number_of_scenarios,
scenario_hour,
scenario_day,
scenario_month,
scenario_year,
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
load_actuals_raw = read_h5_file(here("data", historical_load), "load");
load_forecast_raw = read_h5_file(here("data", "ercot_BA_load_forecast_day_ahead_2018.h5"), "load", false);

# Solar data
solar_actuals_raw = read_h5_file(here("data", "ercot_BA_solar_actuals_Existing_2018.h5"), "solar");
solar_forecast_dayahead_raw = read_h5_file(here("data", "ercot_BA_solar_forecast_day_ahead_existing_2018.h5"), "solar", false);
solar_forecast_2dayahead_raw = read_h5_file(here("data", "ercot_BA_solar_forecast_2_day_ahead_existing_2018.h5"), "solar", false);

# Wind data
wind_actuals_raw = read_h5_file(here("data", "ercot_BA_wind_actuals_Existing_2018.h5"), "wind");
wind_forecast_dayahead_raw = read_h5_file(here("data", "ercot_BA_wind_forecast_day_ahead_existing_2018.h5"), "wind", false);
wind_forecast_2dayahead_raw = read_h5_file(here("data", "ercot_BA_wind_forecast_2_day_ahead_existing_2018.h5"), "wind", false);

#=======================================================================
Compute the hourly average for the actuals data
=======================================================================#
# Load
aux = compute_hourly_average_actuals(load_actuals_raw);
load_actual_avg_raw = DataFrame();
time_index = aux[:, :Index];
avg_actual = aux[:, :values_mean];
load_actual_avg_raw[!, :time_index] = time_index;
load_actual_avg_raw[!, :avg_actual] = avg_actual;

# Solar
aux = compute_hourly_average_actuals(solar_actuals_raw);
time_index = aux[:, :Index];
avg_actual = aux[:, :values_mean];
solar_actual_avg_raw = DataFrame();
solar_actual_avg_raw[!, :time_index] = time_index;
solar_actual_avg_raw[!, :avg_actual] = avg_actual;

# Wind
aux = compute_hourly_average_actuals(wind_actuals_raw);
time_index = aux[:, :Index];
avg_actual = aux[:, :values_mean];
wind_actual_avg_raw = DataFrame();
wind_actual_avg_raw[!, :time_index] = time_index;
wind_actual_avg_raw[!, :avg_actual] = avg_actual;

#=======================================================================
ADJUST THE TIME 
=======================================================================#
#= For the year of 2018, adjust the time to Texas' UTC (UTC-6 or UTC-5)
depending on daylight saving time =#

# Load data
load_actuals = convert_hours_2018(load_actuals_raw);
load_actual_avg = convert_hours_2018(load_actual_avg_raw);
load_forecast = convert_hours_2018(load_forecast_raw, false);

# Solar data
solar_actuals = convert_hours_2018(solar_actuals_raw);
solar_actual_avg = convert_hours_2018(solar_actual_avg_raw);
solar_forecast_dayahead = convert_hours_2018(solar_forecast_dayahead_raw, false);
solar_forecast_2dayahead = convert_hours_2018(solar_forecast_2dayahead_raw, false);

# Wind data
wind_actuals = convert_hours_2018(wind_actuals_raw);
wind_actual_avg = convert_hours_2018(wind_actual_avg_raw);
wind_forecast_dayahead = convert_hours_2018(wind_forecast_dayahead_raw, false);
wind_forecast_2dayahead = convert_hours_2018(wind_forecast_2dayahead_raw, false);

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
DEFINE INDICES, DATETIMES, ISSUE SETS, ETC FOR USE AS INPUTS TO FUNCTIONS
=======================================================================#
# Set the date and time for the forecasts
start_date = DateTime(string(scenario_year) * "-" * string(scenario_month) * "-" * string(scenario_day) * "T" * string(scenario_hour));

wind_times = wind_data[!, [:forecast_time, :issue_time]]
unique_wind = unique(wind_data[!,:issue_time])

solar_times = solar_data[!, [:forecast_time, :issue_time]]
unique_solar = unique(solar_times[!,:issue_time])

load_times = load_data[!, [:forecast_time, :issue_time]]
unique_load = unique(load_times[!,:issue_time])

# XXX load data needs to be fixed, for now just use wind_times
unique_wind == unique_load

times = wind_times
hour_index = findall(x -> x == start_date, times[!,:forecast_time])[1]
hour_index = findall(x -> x == start_date, times)[1]
current_issue = times[hour_index, :issue_time]
# find index of current issue time
issue_index = findall(x -> x == current_issue, unique_issue_times)[1]
# find next issue time and compare to start date
next_issue = unique_issue_times[issue_index + 1]
# define active issues (XXX workshop name) as set of forecasts that are available for the hours to use for lookahead
active_issues = [current_issue]
# get the indices of the forecasts of the active issue times
all_indices = findall(x -> x in active_issues, times[!,:issue_time])


#=======================================================================
Write forecast percentile to files 
=======================================================================#
#write_percentile(load_data, "load", scenario_year, scenario_month, scenario_day, scenario_hour);
write_percentile = 0
if write_percentile == 1
    write_percentile(load_data, "load", scenario_year, scenario_month, scenario_day, scenario_hour)
    write_percentile(solar_data, "solar", scenario_year, scenario_month, scenario_day, scenario_hour)
    write_percentile(wind_data, "wind", scenario_year, scenario_month, scenario_day, scenario_hour)
end


#=======================================================================
Variance of percentiles 
=======================================================================#
#= Brief graphical analysis for the variance of the percentile forecasts.
The idea for this came from a meeting with Jacob Mays on March. 30th,
2023.
=#
print_graphical_analysis = 0
if print_graphical_analysis == 1
    getplots(load_data, "Load data", "load")
    getplots(solar_data, "Solar data", "solar")
    getplots(wind_data, "Wind data", "wind")
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

<<<<<<< Updated upstream:copulas/src/_main.jl
#=======================================================================
CORRELATION HEATMAP FOR THE LANDING PROBABILITIES
=======================================================================#
plot_correlogram = 0;
include(here("src", "functions", "fct_plot_correlogram_landing_probability.jl"));

if plot_correlogram == 1
    plot_correlogram_landing_probability(lp_load, "Load")
    plot_correlogram_landing_probability(lp_solar, "Solar")
    plot_correlogram_landing_probability(lp_wind, "Wind")
end
=======
# #=======================================================================
# CORRELATION HEATMAP FOR THE LANDING PROBABILITIES
# =======================================================================#
# plot_correlogram = true;

# if plot_correlogram
#     plot_correlogram_landing_probability(lp_load, "Load")
#     plot_correlogram_landing_probability(lp_solar, "Solar")
#     plot_correlogram_landing_probability(lp_wind, "Wind")
# end
>>>>>>> Stashed changes:scripts/_main.jl

#=======================================================================
SIMULATE INPUT THROUGH NORTA-LIKE APPROACH
=======================================================================#
load_prob_scen = generate_probability_scenarios(lp_load, scenario_length, number_of_scenarios, start_date);
solar_prob_scen = generate_probability_scenarios(lp_solar, scenario_length, number_of_scenarios, start_date);
wind_prob_scen = generate_probability_scenarios(lp_wind, scenario_length, number_of_scenarios, start_date);

#=======================================================================
CONVERT PROBABILITY SCENARIOS INTO DATA SCENARIOS
=======================================================================#
include(here("src", "functions", "fct_convert_land_prob_to_data_WORKING.jl"));
load_scen = convert_land_prob_to_data_w(load_data, load_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);
solar_scen = convert_land_prob_to_data_w(solar_data, solar_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);
wind_scen = convert_land_prob_to_data_w(wind_data, wind_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);

# #=======================================================================
# PLOT HISTORICAL LANDING
# =======================================================================#
# # Plot the forecasts for each hour and the historical data for same hour
# plot_historical_landing_data = true
# if plot_historical_landing_data
#     plot_historical_landing(load_data, "Load")
#     plot_historical_landing(load_data, "Load", false)
#     plot_historical_landing(solar_data, "Solar")
#     plot_historical_landing(solar_data, "Solar", false)
#     plot_historical_landing(wind_data, "Wind")
#     plot_historical_landing(wind_data, "Wind", false)
# end

<<<<<<< Updated upstream:copulas/src/_main.jl
#=======================================================================
PLOT SYNTHETIC AND HISTORICAL DATA
=======================================================================#
plot_scenarios_and_historical_data = false
include(here("src", "functions", "fct_plot_scenarios_and_actual.jl"))
if plot_scenarios_and_historical_data
    plot_scenarios_and_actual(load_scen, load_data, "Load", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual(solar_scen, solar_data, "Solar", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual(wind_scen, wind_data, "Wind", scenario_year, scenario_month, scenario_day, scenario_hour)
end

#=======================================================================
PLOT SYNTHETIC AND HISTORICAL DATA --- WORKING
=======================================================================#
include(here("src", "functions", "fct_plot_scenarios_and_actual_WORKING.jl"))
plot_scenarios_and_historical_data = true
if plot_scenarios_and_historical_data
    plot_scenarios_and_actual_w(load_actuals, load_scen, load_data, "Load", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual_w(solar_actuals, solar_scen, solar_data, "Solar", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual_w(wind_actuals, wind_scen, wind_data, "Wind", scenario_year, scenario_month, scenario_day, scenario_hour)
end

#=======================================================================
PLOT SYNTHETIC AND HISTORICAL AUTOCORRELATION
=======================================================================#
plot_autocorrelation = true
#scenarios, actuals, type, scenario_year, scenario_month, scenario_day, scenario_hour
include(here("src", "functions", "fct_plot_historical_synthetic_autocorrelation.jl"));
if plot_autocorrelation
    plot_hist_synth_autocor(load_scen, load_data, "Load",scenario_year, scenario_month, scenario_day, scenario_hour);
    plot_hist_synth_autocor(solar_scen, solar_data, "Solar",scenario_year, scenario_month, scenario_day, scenario_hour);
    plot_hist_synth_autocor(wind_scen, wind_data, "Wind",scenario_year, scenario_month, scenario_day, scenario_hour);
end

#=======================================================================
PLOT SYNTHETIC AND HISTORICAL PARTIAL AUTOCORRELATION
=======================================================================#
plot_partial_autocorrelation = false



#=======================================================================
CREATE DATA MATRIX FOR INTEGRATION WITH GEN X MODEL
=======================================================================#
myrange = collect(1:number_of_scenarios)

# Load -----------------------------------------------------------------
load_names = "Load_MW_z1_" .* string.(myrange);
load_header = [
    "Voll",
    "Demand_Segment",
    "Cost_of_Demand_Curtailment_per_MW",
    "Max_Demand_Curtailment",
    "Rep_Periods",
    "Timesteps_per_Rep_Period",
    "Sub_Weights",
    "Time_Index",
    "Load_MW_z1",
];
load_header = vcat(load_header, load_names);

load_genx = zeros(scenario_length, length(load_header));
load_genx[:, 10:length(load_header)] = transpose(load_scen);
load_genx = DataFrame(load_genx, :auto);
rename!(load_genx, load_header);
CSV.write("load_data_scenarios.csv", load_genx);

# Renewables -----------------------------------------------------------

# RENEWABLES NEED TO BE IMPROVED/CORRECTED. 

# renewable_df = zeros(scenario_length, 2 * number_of_scenarios);
# renewable_df[:, solar_index] = transpose(copy(solar_scen));
# renewable_df[:, wind_index] = transpose(copy(wind_scen));

# column_names = vec(["solar_pv_", "onshore_wind_"] .* string.(myrange'));
# renewable_header = [
=======
# #=======================================================================
# PLOT SYNTHETIC AND HISTORICAL DATA
# =======================================================================#
# plot_scenarios_and_historical_data = true
# if plot_scenarios_and_historical_data
#     plot_scenarios_and_actual(load_actuals, load_scen, load_data, "Load", scenario_year, scenario_month, scenario_day, scenario_hour)
#     plot_scenarios_and_actual(solar_actuals, solar_scen, solar_data, "Solar", scenario_year, scenario_month, scenario_day, scenario_hour)
#     plot_scenarios_and_actual(wind_actuals, wind_scen, wind_data, "Wind", scenario_year, scenario_month, scenario_day, scenario_hour)
# end

# #=======================================================================
# PLOT SYNTHETIC AND HISTORICAL AUTOCORRELATION
# =======================================================================#
# plot_autocorrelation = true
# if plot_autocorrelation
#     plot_hist_synth_autocor(load_scen, load_data, "Load",scenario_year, scenario_month, scenario_day, scenario_hour);
#     plot_hist_synth_autocor(solar_scen, solar_data, "Solar",scenario_year, scenario_month, scenario_day, scenario_hour);
#     plot_hist_synth_autocor(wind_scen, wind_data, "Wind",scenario_year, scenario_month, scenario_day, scenario_hour);
# end

# #=======================================================================
# PLOT SYNTHETIC AND HISTORICAL PARTIAL AUTOCORRELATION
# =======================================================================#
# # TBD


# #=======================================================================
# CREATE DATA MATRIX FOR INTEGRATION WITH GEN X MODEL
# =======================================================================#
# myrange = collect(1:number_of_scenarios);

# # Load -----------------------------------------------------------------
# load_names = "Load_MW_z1_" .* string.(myrange);
# load_header = [
#     "Voll",
#     "Demand_Segment",
#     "Cost_of_Demand_Curtailment_per_MW",
#     "Max_Demand_Curtailment",
#     "Rep_Periods",
#     "Timesteps_per_Rep_Period",
#     "Sub_Weights",
>>>>>>> Stashed changes:scripts/_main.jl
#     "Time_Index",
#     "Load_MW_z1",
# ];
# load_header = vcat(load_header, load_names);

# load_genx = zeros(scenario_length, length(load_header));
# load_genx[:, 10:length(load_header)] = transpose(load_scen);
# load_genx = DataFrame(load_genx, :auto);
# rename!(load_genx, load_header);
# filepath = mkpath(datadir("exp_pro","gen_x_integration"));
# CSV.write(joinpath(filepath,"load_data_scenarios.csv"), load_genx);

<<<<<<< Updated upstream:copulas/src/_main.jl
=======
# # Renewables -----------------------------------------------------------
# # RENEWABLES NEED TO BE IMPROVED/CORRECTED. 

# renewable_df = zeros(scenario_length, 2 * number_of_scenarios);
# # renewable_df[:, solar_index] = transpose(copy(solar_scen));
# # renewable_df[:, wind_index] = transpose(copy(wind_scen));

# # column_names = vec(["solar_pv_", "onshore_wind_"] .* string.(myrange'));
# # renewable_header = [
# #     "Time_Index",
# #     "natural_gas_combined_cycle",
# # ]

# # renewable_header = vcat(renewable_header, column_names, "battery")
# # renewable_genx = zeros(scenario_length, length(renewable_header));
# # renewable_genx[:, 3:length(renewable_header)-1] = renewable_df;
# # renewable_genx = DataFrame(renewable_genx, :auto);
# # rename!(renewable_genx, renewable_header);
# # CSV.write("generators_variability_scenarios.csv", renewable_genx)

# ### Generate single dataframe with concatenated scenarios of Solar, Wind, and Load, in that order
# master_scen_mat = zeros(scenario_length, 3 * number_of_scenarios)
# master_scen_df = DataFrame()

# # loop over the number of scenarios and concatenate scenario to the master scenario datafarme
# for scen_id in range(1,number_of_scenarios)
#     id_solar = copy(solar_scen[scen_id,:]) # using transpose does not work
#     id_wind = copy(wind_scen[scen_id,:])
#     id_load = copy(load_scen[scen_id,:])
#     insertcols!(master_scen_df, :S => id_solar, :W => id_wind, :L => id_load, makeunique=true )
# end
>>>>>>> Stashed changes:scripts/_main.jl
