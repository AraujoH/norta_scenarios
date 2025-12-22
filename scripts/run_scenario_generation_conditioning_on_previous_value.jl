# Scenario Generation with Copulas 
# 
# Hugo S. de Araujo
# Nov. 14th, 2022 | Mays Group | Cornell University
################################################################################

#=======================================================================
PROJECT SETUP
=======================================================================#
using Pkg
Pkg.activate("norta_scenarios")

# Import all required packages. 
begin
    using CSV
    using DataFrames
    using Dates
    using DelimitedFiles
    using Distributions
    using HDF5
    # using LaTeXStrings
    using LinearAlgebra
    using LinearSolve
    using Random
    using RCall
    using Statistics
    using StatsBase
    using Plots
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
include(joinpath(pwd(), "src", "fct_getplots.jl"));
# include(joinpath(pwd(), "src", "fct_plot_correlation_heatmap.jl"));
include(joinpath(pwd(), "src", "fct_plot_historical_landing.jl"));
include(joinpath(pwd(), "src", "fct_plot_historical_synthetic_autocorrelation.jl"));
include(joinpath(pwd(), "src", "fct_plot_correlogram_landing_probability.jl"));
include(joinpath(pwd(), "src", "fct_plot_scenarios_and_actual.jl"));
include(joinpath(pwd(), "src", "fct_read_h5_file.jl"));
include(joinpath(pwd(), "src", "fct_read_input_file.jl"));
include(joinpath(pwd(), "src", "fct_transform_landing_probability.jl"));
include(joinpath(pwd(), "src", "fct_write_percentiles.jl"));
include(joinpath(pwd(), "src", "fct_write_scenarios.jl"));
include(joinpath(pwd(), "src", "fct_plot_scenarios_and_actual_aExam.jl"));
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
load_prob_scen = generate_probability_scenarios(
    lp_load, scenario_length, number_of_scenarios; intraday_market_scenarios=true);
solar_prob_scen = generate_probability_scenarios(
    lp_solar, scenario_length, number_of_scenarios; intraday_market_scenarios=true);
wind_prob_scen = generate_probability_scenarios(
    lp_wind, scenario_length, number_of_scenarios; intraday_market_scenarios=true);

#=======================================================================
CONVERT PROBABILITY SCENARIOS INTO DATA SCENARIOS
=======================================================================#
load_scen = convert_land_prob_to_data(
    load_data, load_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);
solar_scen = convert_land_prob_to_data(
    solar_data, solar_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);
wind_scen = convert_land_prob_to_data(
    wind_data, wind_prob_scen, scenario_year, scenario_month, scenario_day, scenario_hour);

#=======================================================================
WRITE SCENARIOS TO FILE
=======================================================================#
write_scenarios(load_scen, "load_idm_50batches_scenarios")
write_scenarios(solar_scen, "solar_idm_50batches_scenarios")
write_scenarios(wind_scen, "wind_idm_50batches_scenarios")

#=======================================================================
PLOT HISTORICAL LANDING
=======================================================================#
# Plot the forecasts for each hour and the historical data for same hour
plot_historical_landing_data = true
if plot_historical_landing_data
    plot_historical_landing(load_data, "Load")
    plot_historical_landing(load_data, "Load", false)
    plot_historical_landing(solar_data, "Solar")
    plot_historical_landing(solar_data, "Solar", false)
    plot_historical_landing(wind_data, "Wind")
    plot_historical_landing(wind_data, "Wind", false)
end

#=======================================================================
PLOT SYNTHETIC AND HISTORICAL DATA
=======================================================================#
plot_scenarios_and_historical_data = true
if plot_scenarios_and_historical_data
    plot_scenarios_and_actual(load_actuals, load_scen, load_data, "Load", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual(solar_actuals, solar_scen, solar_data, "Solar", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual(wind_actuals, wind_scen, wind_data, "Wind", scenario_year, scenario_month, scenario_day, scenario_hour)
end


#=======================================================================
PLOT SYNTHETIC AND HISTORICAL DATA FOR A-EXAM PRESENTATION
=======================================================================#
plot_scenarios_and_historical_data = true
if plot_scenarios_and_historical_data
    plot_scenarios_and_actual_aExam(load_actuals, load_scen, load_data, "Load", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual_aExam(solar_actuals, solar_scen, solar_data, "Solar", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_scenarios_and_actual_aExam(wind_actuals, wind_scen, wind_data, "Wind", scenario_year, scenario_month, scenario_day, scenario_hour)
end


#=======================================================================
PLOT SYNTHETIC AND HISTORICAL AUTOCORRELATION
=======================================================================#
plot_autocorrelation = true
if plot_autocorrelation
    plot_hist_synth_autocor(load_scen, load_data, "Load", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_hist_synth_autocor(solar_scen, solar_data, "Solar", scenario_year, scenario_month, scenario_day, scenario_hour)
    plot_hist_synth_autocor(wind_scen, wind_data, "Wind", scenario_year, scenario_month, scenario_day, scenario_hour)
end

#=======================================================================
PLOT SYNTHETIC AND HISTORICAL PARTIAL AUTOCORRELATION
=======================================================================#
# TBD


#=======================================================================
CREATE DATA MATRIX FOR INTEGRATION WITH GEN X MODEL
=======================================================================#
myrange = collect(1:number_of_scenarios);

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
filepath = mkpath(datadir("exp_pro", "gen_x_integration"));
CSV.write(joinpath(filepath, "load_data_scenarios.csv"), load_genx);

# Renewables -----------------------------------------------------------
# RENEWABLES NEED TO BE IMPROVED/CORRECTED. 

renewable_df = zeros(scenario_length, 2 * number_of_scenarios);
# renewable_df[:, solar_index] = transpose(copy(solar_scen));
# renewable_df[:, wind_index] = transpose(copy(wind_scen));

# column_names = vec(["solar_pv_", "onshore_wind_"] .* string.(myrange'));
# renewable_header = [
#     "Time_Index",
#     "natural_gas_combined_cycle",
# ]

# renewable_header = vcat(renewable_header, column_names, "battery")
# renewable_genx = zeros(scenario_length, length(renewable_header));
# renewable_genx[:, 3:length(renewable_header)-1] = renewable_df;
# renewable_genx = DataFrame(renewable_genx, :auto);
# rename!(renewable_genx, renewable_header);
# CSV.write("generators_variability_scenarios.csv", renewable_genx)

### Generate single dataframe with concatenated scenarios of Solar, Wind, and Load, in that order
master_scen_mat = zeros(scenario_length, 3 * number_of_scenarios)
master_scen_df = DataFrame()

# loop over the number of scenarios and concatenate scenario to the master scenario datafarme
for scen_id in range(1, number_of_scenarios)
    id_solar = copy(solar_scen[scen_id, :]) # using transpose does not work
    id_wind = copy(wind_scen[scen_id, :])
    id_load = copy(load_scen[scen_id, :])
    insertcols!(master_scen_df, :S => id_solar, :W => id_wind, :L => id_load, makeunique=true)
end
