function compute_hourly_average_actuals(actual_data)
    """
    compute_hourly_average_actuals: take the average of the actual/
        historical data. Historical data has 5-min granularity. The 
        average sets the length (number of entries) of the historical
        data eqaul to the forecasts, existing only for exact hours.
        
    # Arguments:
    - actual_data::String - variable name of the dataframe with historical
        data
    
    # Returns:
    - Dataframe df with averaged historical values for on-the-hour values
    """
    df = actual_data
    df = TSFrame(df)
    return (apply(df, Hour(1), mean))
end