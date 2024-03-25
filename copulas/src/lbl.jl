# Scenario Generation with Copulas 
# 
# Hugo S. de Araujo
# Nov. 14th, 2022 | Mays Group | Cornell University
# Edited: Kaleb Smith
# Jan. 10th, 2024
################################################################################

#=======================================================================
PROJECT SETUP
=======================================================================#
using DrWatson
@quickactivate("norta_scenarios")

# Import all required packages. 
begin
    using CSV
    using DataFrames
    using Dates
    using DelimitedFiles
    using Distributions
    using HDF5
    using LaTeXStrings
    using LinearAlgebra
    using LinearSolve
    using Random
    using Statistics
    using StatsBase
    using Plots
    using Tables
    using TSFrames
    using TimeZones
end

# Include functions 
include(projectdir("src", "fct_bind_historical_forecast.jl"));
include(projectdir("src", "fct_compute_hourly_average_actuals.jl"));
include(projectdir("src", "fct_compute_landing_probability.jl"));
include(projectdir("src", "fct_convert_hours_2018.jl"));
include(projectdir("src", "fct_convert_ISO_standard.jl"));
include(projectdir("src", "fct_convert_land_prob_to_data.jl"));
include(projectdir("src", "fct_generate_probability_scenarios.jl"));
#include(projectdir("src", "fct_getplots.jl"));
#include(projectdir("src", "fct_plot_correlation_heatmap.jl"));
include(projectdir("src", "fct_plot_historical_landing.jl"));
include(projectdir("src", "fct_plot_historical_synthetic_autocorrelation.jl"));
include(projectdir("src", "fct_plot_correlogram_landing_probability.jl"));
include(projectdir("src", "fct_plot_scenarios_and_actual.jl"));
include(projectdir("src", "fct_read_h5_file.jl"));
include(projectdir("src", "fct_read_input_file.jl"));
include(projectdir("src", "fct_transform_landing_probability.jl"));
include(projectdir("src", "fct_write_percentiles.jl"));

#=======================================================================
READ INPUT FILE
=======================================================================#
input_file_path = projectdir("copulas.txt")

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
load_actuals = read_h5_file(datadir("exp_raw", historical_load), "load");
load_forecast = read_h5_file(datadir("exp_raw", "ercot_BA_load_forecast_day_ahead_2018.h5"), "load", false);

# Solar data
solar_actuals = read_h5_file(datadir("exp_raw", "ercot_BA_solar_actuals_Existing_2018.h5"), "solar");
solar_forecast_dayahead = read_h5_file(datadir("exp_raw", "ercot_BA_solar_forecast_day_ahead_existing_2018.h5"), "solar", false);
solar_forecast_2dayahead = read_h5_file(datadir("exp_raw", "ercot_BA_solar_forecast_2_day_ahead_existing_2018.h5"), "solar", false);

# Wind data
wind_actuals = read_h5_file(datadir("exp_raw","ercot_BA_wind_actuals_Existing_2018.h5"), "wind");
wind_forecast_dayahead = read_h5_file(datadir("exp_raw", "ercot_BA_wind_forecast_day_ahead_existing_2018.h5"), "wind", false);
wind_forecast_2dayahead = read_h5_file(datadir("exp_raw", "ercot_BA_wind_forecast_2_day_ahead_existing_2018.h5"), "wind", false);

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
# load_data = bind_historical_forecast(true,
#     load_actual_avg,
#     load_forecast);
is_load,
historical_data,
forecast_da_data,
forecast_2da_data=true, load_actual_avg, load_forecast, "none"

ahead_factor = repeat(["two", "one"], size(forecast_da_data, 1) ÷ 2)
forecast_da_data[!, :ahead_factor] = ahead_factor
load_da_data = filter(:ahead_factor => ==("one"), forecast_da_data)
load_2da_data = filter(:ahead_factor => ==("two"), forecast_da_data)
full_data = leftjoin(forecast_da_data, historical_data, on=[:forecast_time => :time_index])

one_data = leftjoin(load_da_data, historical_data, on=[:forecast_time => :time_index]) # XXX something weird happens here

# XXX weird becase the time index and forecast time arrays for load da and historical data are identical
load_da_data[:,:forecast_time] == historical_data[:, :time_index]
# XXX so the leftjoin should be 8760...

# XXX check to see what is included in one_data now that is not in load_da_data

one_data[one_data[:,:forecast_time] .∉ Ref(load_da_data[:,:forecast_time]),:forecast_time]
test = one_data[:,:forecast_time][(!in).(one_data[:,:forecast_time],Ref(load_da_data[:,:forecast_time]))]
setdiff(one_data[:,:forecast_time], load_da_data[:,:forecast_time])
setdiff(load_da_data[:,:forecast_time], one_data[:,:forecast_time])
# XXX none of these indicate any issue...
# is it possible that it is on the issue time and not the forecast time?
setdiff(load_da_data[:,:issue_time], one_data[:,:issue_time])
# XXX try unique... see if that gets rid of one of the 8764, could have duplicates....
unique(one_data[:,:forecast_time])
# XXX why are there only 8758 unique elements when there are 8764 in the dataframe
unique(load_da_data[:,:forecast_time])
# XXX but load_da_data has the same number of unique....
# I think the it is the daylight savings time is causing some duplicates
# XXX but would that affect the leftjoin? potentially yes becasue load might be done differently


# XXX try extract duplicates of forecast time




two_data = leftjoin(load_2da_data, historical_data, on=[:forecast_time => :time_index])

load_data = full_data
load_1d_data = one_data

solar_data = bind_historical_forecast(false,
    solar_actual_avg,
    solar_forecast_dayahead,
    solar_forecast_2dayahead);

wind_data = bind_historical_forecast(false,
    wind_actual_avg,
    wind_forecast_dayahead,
    wind_forecast_2dayahead);


data, type, 
scenario_year, 
scenario_month, 
scenario_day, 
scenario_hour = load_data, "load", scenario_year, scenario_month, scenario_day, scenario_hour

# Set the date and time for the forecasts
start_date = DateTime(string(scenario_year) * "-" * string(scenario_month) * "-" * string(scenario_day) * "T" * string(scenario_hour));


# calculate lookahead horizon length
correlated_times = solar_forecast_dayahead[!,[:forecast_time, :issue_time]]
# ensure identical to wind forecast wind_forecast_dayahead
correlated_times == wind_forecast_dayahead[!,[:forecast_time, :issue_time]]
# ensure identical to extracted load forecast wind_forecast_dayahead
correlated_times == load_1d_data[!,[:forecast_time, :issue_time]] #XXX this should be true
# get list of unique issue times from correlated times
unique_issue_times = unique(correlated_times[!,:issue_time])
# get list of forecast times too
forecast_times = load_1d_data[!,:forecast_time]
# find index of scenario time=forecast_time
hour_index = findall(x -> x == start_date, forecast_times)[1]
# find issue_time of current scenario_time=forecast_time
current_issue = correlated_times[hour_index, :issue_time]
# find index of current issue time
issue_index = findall(x -> x == current_issue, unique_issue_times)[1]
# find next issue time and compare to start date
next_issue = unique_issue_times[issue_index + 1]
# has the next issue time occurred in the real time (RT) rolling horizon index?
# if true, start_date is after next issue date and therefore next forecast can be included
start_date > next_issue
### calculate the lookahead duration
# define active issues (XXX workshop name) as set of forecasts that are available for the hours to use for lookahead
active_issues = [current_issue]

# if 
if start_date > next_issue
    # push!(active_issues, next_issue) # use this if 
    active_issues = [next_issue]
else
    active_issues = [current_issue]
end

### get the indices of the forecasts of the active issue times
all_indices = findall(x -> x in active_issues, correlated_times[!,:issue_time])
all_times = forecast_times[all_indices]

# filter for the times that are after the current time, (greater than the current hour)
active_times = all_times[all_times .>= start_date]

# calculate length as the lookahead length
lookahead_length = length(active_times)

end_date = start_date + Hour(lookahead_length)

#=======================================================================
Write forecast percentile to files 
=======================================================================#
#write_percentile(load_data, "load", scenario_year, scenario_month, scenario_day, scenario_hour);
write_percentile = true
if write_percentile
    write_percentiles(load_data, "load", scenario_year, scenario_month, scenario_day, scenario_hour, lookahead_length)
    write_percentiles(solar_data, "solar", scenario_year, scenario_month, scenario_day, scenario_hour, lookahead_length)
    write_percentiles(wind_data, "wind", scenario_year, scenario_month, scenario_day, scenario_hour, lookahead_length)
end

# landing_probability_load = compute_landing_probability(load_1d_data);
percentile_column_index = startswith.(names(load_1d_data), "p_");
landing_probability = Vector{Float64}(undef, size(load_1d_data, 1));

for i in range(1, size(load_1d_data, 1))
    quantiles = collect(load_1d_data[i, percentile_column_index]);
    empirical_cdf = ecdf(quantiles);
    landing_probability[i] = empirical_cdf(load_1d_data[i, :avg_actual]); # XXX why is this based on the average actual?
end


# Create a DataFrame with the forecast time, issue time, ahead factor
# and the landing probability.
landing_probability = DataFrame(landing_probability=landing_probability)
landing_probability[!, :forecast_time] = load_1d_data[!, :forecast_time]
landing_probability[!, :issue_time] = load_1d_data[!, :issue_time]
landing_probability[!, :ahead_factor] = load_1d_data[:, :ahead_factor]
select!(landing_probability, [:issue_time, :forecast_time, :landing_probability, :ahead_factor])

landing_probability_load = landing_probability



landing_probability_solar = compute_landing_probability(solar_data);
landing_probability_wind = compute_landing_probability(wind_data);

# lp_load = transform_landing_probability(landing_probability_load);
# lp_solar = transform_landing_probability(landing_probability_solar);
# lp_wind = transform_landing_probability(landing_probability_wind);

x_data = landing_probability_load
sort!(x_data, :issue_time)

df = combine(groupby(x_data, [:issue_time]), DataFrames.nrow => :count)

# df_filtered = filter(:count => ==(48), df); # this is no longer 48 
# because I separated the day ahead and two day ahead forecasts
# also, this would neglect the initial and final forecasts

# if it were 24, it would look like this and only keep the rows 2 through 365
df_filtered = filter(:count => ==(24), df)
# actually, I think the daylights savings and back is incorrect, because these have been removed too

issue_times_interest = df_filtered[!, :issue_time];
landing_probability_filtered = innerjoin(x_data, df_filtered, on=:issue_time);

# 
landing_probability_filtered_matrix = reshape(landing_probability_filtered[!, :landing_probability], (24, size(df_filtered, 1)));
landing_probability_filtered_matrix = transpose(landing_probability_filtered_matrix);

# end of transform landing probability

# yeah this is not going to work. The process has to be based on my active_issues and extracting, 
# joining, reshaping to the lookahead length, etc



### generate_probability_scenarios
# scenario_legth = lookahead_length
scenario_length = 24

 ### XXX What is this supposed to be doing???
x = copy(landing_probability_filtered_matrix)
allequal_set = Set(findall(allequal, eachcol(x)));
allequal_ind = sort(collect(allequal_set));
allindex_set = Set(collect(1:24));
alldifferent_ind = sort(collect(setdiff(allindex_set, allequal_set))); # Index for columns whose st. dev. isn't zero
x = x[:, alldifferent_ind];

if ishermitian(cor(x))
    Σ_Z = LinearAlgebra.cholesky(cor(x));
else
    Σ_Z = factorize(cor(x));
end
M = Σ_Z.L;
dim_M = size(M, 1);

Random.seed!(12345)
Y = Matrix{Float64}(undef, number_of_scenarios, scenario_length)

need_expansion = 0 # This is specially important for solar data with several columns whose st. dev. is zero
if dim_M != scenario_length # why would dim_M ever not be equal to scenario length
    original_scen_length = scenario_length
    scenario_length = dim_M
    Y = Matrix{Float64}(undef, number_of_scenarios, scenario_length)
    need_expansion = 1 # what is this need expansion mean?
end

for nscen in 1:number_of_scenarios
    # Set up normal distribution with mean 0 and sd equal to 1
    d = Normal(0,1);

    #Generate vector W = (W_1, ..., W_k) ~ iid standard normal
    W = rand(d, scenario_length);

    # Create vector Z such that Z <- MW
    Z = M * W;

    #Compute the CDF of Z
    #cdf_Z = sort(cdf.(d, Z));
    cdf_Z = cdf.(d, Z);
    
    for k in 1:scenario_length
        #Apply the inverse CDF for X_k
        # Y[nscen, k] = quantile(x[:, k], cdf_Z[k])
        Y[nscen, k] = quantile(x[:, k], cdf_Z[k]);
    end
end

#= If there is the need for expansion, expansion is done in the 
following block of code.
=#
if need_expansion == 1
    Y_aux = Matrix{}(undef, number_of_scenarios, original_scen_length)
    Y_aux[:, allequal_ind] .= 0
    Y_aux[:, alldifferent_ind] .= Y
else
    # do nothing
end