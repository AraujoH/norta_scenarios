function read_h5_file(filepath::String,
    type::String=["load", "wind", "solar"],
    actuals::Bool=true)
    """
    read_h5_file: read .h5 files parsing through the type of data being
                  read

    # Arguments
    - `filepath::String`: path to the h5 file
    - `type::String`: name of the data type to be read: (load, wind, solar)
    - 'actuals::Bool': if true, data is historical, not forecast

    # Returns
    Dataframe with data and its time-related information. If actuals = 
    true, time index. If actuals = false, forecast time and issue time.
    """
    fid = h5open(filepath)
    x = read(fid)

    # Create vector to name columns of forecast-valued dataframes
    column_names = Vector{String}(undef, 101)
    for i in 1:101
        if i == 1
            column_names[i] = "forecast_time"
        elseif i == 2
            column_names[i] = "issue_time"
        else
            column_names[i] = string("p_", i - 2)
        end
    end


    # Load data --------------------------------------------------------
    if type == "load"
        if actuals == true
            values = x["actuals"]
            timestamp = x["time_index"]
            timestamp = convert_ISO_standard(timestamp)

            df = DataFrame(time_index=timestamp,
                values=values)
        else
            values = transpose(x["forecasts"])
            forecast_time = x["forecast_time"]
            forecast_time = convert_ISO_standard(forecast_time)
            issue_time = x["issue_time"]
            issue_time = convert_ISO_standard(issue_time)

            df = DataFrame(values, :auto) #:auto to automatically name columns
            insertcols!(df, 1, :forecast_time => forecast_time)
            insertcols!(df, 2, :issue_time => issue_time)
            rename!(df, column_names)
        end
    else
        # Solar and Wind data ------------------------------------------
        if actuals == true
            values = x["actuals"]
            values = transpose(values)[:, 1]

            timestamp = x["time_index"]
            timestamp = convert_ISO_standard(timestamp)

            df = DataFrame(time_index=timestamp,
                values=values)
        else
            values = x["forecasts"]
            forecast_time = x["forecast_time"]
            forecast_time = convert_ISO_standard(forecast_time)
            issue_time = x["issue_time"]
            issue_time = convert_ISO_standard(issue_time)

            df = DataFrame(values, :auto) #:auto to automatically name columns
            insertcols!(df, 1, :forecast_time => forecast_time)
            insertcols!(df, 2, :issue_time => issue_time)
            rename!(df, column_names)
        end
    end
    return (df)
end