function bind_historical_forecast(
    is_load,
    historical_data,
    forecast_da_data,
    forecast_2da_data="none")
    """
    bind_historical_forecast: bind the historical averaged data to the 
        forecast data.

    # Arguments
    - historical_data::Dataframe - dataframe with the hourly averaged 
        historical data
    - forecast_da_data::Dataframe - dataframe with the forecast data. 
        If data is load, then only this one is necessary
    - forecast_2da_data::Dataframe - dataframe with the two-day ahead
        data. Necessary if data is either solar or wind

    # Returns
    - full_data::Dataframe - dataframe with bound historical average data
        and forecast data
    """
    if is_load
        ahead_factor = repeat(["two", "one"], size(forecast_da_data, 1) ÷ 2)
        forecast_da_data[!, :ahead_factor] = ahead_factor
        full_data = leftjoin(forecast_da_data, historical_data, on=[:forecast_time => :time_index])
    else
        forecast_da_data[!, :ahead_factor] = repeat(["one"], size(forecast_da_data, 1))
        forecast_2da_data[!, :ahead_factor] = repeat(["two"], size(forecast_2da_data, 1))
        forecast = [forecast_da_data; forecast_2da_data]
        full_data = leftjoin(forecast, historical_data, on=[:forecast_time => :time_index])
    end
    return full_data
end