function convert_hours_2018(data, is_actual = true)
    """
    """
    #= In 2018, daylight savings time in the US began on Mar. 11 at 02:00h, 
    and ended on Nov. 4 at 02:00h. In the UK, it began on Mar. 25 at 01:00h;
    it ended on Oct. 28 at 02:00h. 

    =====================  London ====== Houston |             
    Jan. 01 ---> Mar. 11 | UTC+0  ====== UTC-6   | 
    Mar. 11 ---> Mar. 25 | UTC+0  ====== UTC-5   | 
    Mar. 25 ---> Oct. 28 | UTC+1  ====== UTC-5   | 
    Oct. 28 ---> Nov. 04 | UTC+0  ====== UTC-5   | 
    Nov. 04 ---> Dec. 31 | UTC+0  ====== UTC-6   | 
    ----------------------------------------------
    =#

    # Find the time (Y-M-D-H) limits when dayligh savings time changes 
    # and break it in groups. Since forecasts have forecast time and 
    # issue time, forecasts have been classified into two groups. 
    # The issue time takes the forecast "y" after the group number.
    if is_actual
        x = copy(data[:, :time_index])
        group_1 = x .<= DateTime(2018, 03, 11, 2) #from Jan. till Mar. 11
        group_2 = DateTime(2018, 03, 11, 2) .< x .<= DateTime(2018, 03, 25, 1) #from Mar. 11 to Mar. 25
        group_3 = DateTime(2018, 03, 25, 1) .< x .<= DateTime(2018, 10, 28, 2) #from Mar. 25 to Oct. 28
        group_4 = DateTime(2018, 10, 28, 2) .< x .<= DateTime(2018, 11, 04, 2) #from Oct. 28 to Nov. 04
        group_5 = x .>= DateTime(2018, 11, 04, 2) #from Nov. 04 ahead
    else
        x = copy(data[:, :forecast_time])
        group_1 = x .<= DateTime(2018, 03, 11, 2) #from Jan. till Mar. 11
        group_2 = DateTime(2018, 03, 11, 2) .< x .<= DateTime(2018, 03, 25, 1) #from Mar. 11 to Mar. 25
        group_3 = DateTime(2018, 03, 25, 1) .< x .<= DateTime(2018, 10, 28, 2) #from Mar. 25 to Oct. 28
        group_4 = DateTime(2018, 10, 28, 2) .< x .<= DateTime(2018, 11, 04, 2) #from Oct. 28 to Nov. 04
        group_5 = x .>= DateTime(2018, 11, 04, 2) #from Nov. 04 ahead

        y = copy(data[:, :issue_time])
        group_1y = y .<= DateTime(2018, 03, 11, 2) #from Jan. till Mar. 11
        group_2y = DateTime(2018, 03, 11, 2) .< y .<= DateTime(2018, 03, 25, 1) #from Mar. 11 to Mar. 25
        group_3y = DateTime(2018, 03, 25, 1) .< y .<= DateTime(2018, 10, 28, 2) #from Mar. 25 to Oct. 28
        group_4y = DateTime(2018, 10, 28, 2) .< y .<= DateTime(2018, 11, 04, 2) #from Oct. 28 to Nov. 04
        group_5y = y .>= DateTime(2018, 11, 04, 2) #from Nov. 04 ahead
    end

    # Groups are set up. Carry out transformations
    if is_actual
        data[group_1, :time_index] = data[group_1, :time_index] - Hour(6)
        data[group_2, :time_index] = data[group_2, :time_index] - Hour(5)
        data[group_3, :time_index] = data[group_3, :time_index] - Hour(6)
        data[group_4, :time_index] = data[group_4, :time_index] - Hour(5)
        data[group_5, :time_index] = data[group_5, :time_index] - Hour(6)
    else
        # Forecast
        data[group_1, :forecast_time] = data[group_1, :forecast_time] - Hour(6)
        data[group_2, :forecast_time] = data[group_2, :forecast_time] - Hour(5)
        data[group_3, :forecast_time] = data[group_3, :forecast_time] - Hour(6)
        data[group_4, :forecast_time] = data[group_4, :forecast_time] - Hour(5)
        data[group_5, :forecast_time] = data[group_5, :forecast_time] - Hour(6)

        # Issue
        data[group_1y, :issue_time] = data[group_1y, :issue_time] - Hour(6)
        data[group_2y, :issue_time] = data[group_2y, :issue_time] - Hour(5)
        data[group_3y, :issue_time] = data[group_3y, :issue_time] - Hour(6)
        data[group_4y, :issue_time] = data[group_4y, :issue_time] - Hour(5)
        data[group_5y, :issue_time] = data[group_5y, :issue_time] - Hour(6)
    end
    return(data)
end