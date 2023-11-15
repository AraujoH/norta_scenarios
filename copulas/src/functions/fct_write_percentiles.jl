function write_percentiles(data, type, scenario_year, scenario_month, scenario_day, scenario_hour)
    """
    # Return 
    Input:
        -
        -
    Output:
        -
    """
    # Set the date and time for the forecasts
    start_date = DateTime(string(scenario_year) * "-" * string(scenario_month) * "-" * string(scenario_day) * "T" * string(scenario_hour));
    end_date = start_date + Hour(47)

    # Filter dataframe by dates
    x = data[ start_date .<= data.forecast_time .<= end_date ,:];

    # Save dataframe
    filename = type * "_" * string(scenario_year) * "_" * string(scenario_month) *"_"* string(scenario_day) * ".csv";
    CSV.write(joinpath(pwd(),"output",filename), x)
end