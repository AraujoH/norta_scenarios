function plot_historical_landing(data, type, dayahead = true)

    # Selecting data ---------------------------------------------------
    x = copy(data)

    if dayahead
        df = x[x.ahead_factor.=="one", :]
    else
        df = x[x.ahead_factor.=="two", :]
    end

    # Count the number of times the historical landed outside of 
    # the forecast range
    number_of_misses = 0
    for i in 1:size(df,1)
        if (df.avg_actual[i] .< (df[i,3])) + (df.avg_actual[i] .> (df[i,101])) !=0
            number_of_misses = number_of_misses + 1
        end
    end
    annotation = string(number_of_misses)*" misses"

    # Plot parameters --------------------------------------------------
    # Limit y-axis 
    y_ub = 1.1 * maximum(maximum.(eachcol(df[:, 3:101])))

    # Title 
    if dayahead
        title = "Historical and day-ahead forecasts | " * type
    else
        title = "Historical and 2 day-ahead forecasts | " * type
    end

    # Date format
    time_ticks = range(first(df.forecast_time), last(df.forecast_time), step = Month(1))
    ticks = Dates.format.(time_ticks,"u/dd")

    # Plotting -------------------------------------------------------------
    # Create plot background
    plot(fontfamily="Computer Modern",
        ylim=(0, y_ub),
        xticks=(time_ticks,ticks),
        grid=false,
        framestyle=:box,
        titlelocation=:left,
        titlefontsize=18,
        guidefontsize=16,
        tickfontsize=12,
        legendfontsize=12,
        legend=:outerbottom,
        legendcolumns=2,
        size=(1400, 800),
        foreground_color_legend=nothing)

    # Add forecasts 
    scatter!(
        df.forecast_time,
        df[:,3],
        mc=:whitesmoke,
        markersize=0.8,
        markerstrokewidth=0,
        label = "Forecasts")

    for i in 4:101
        scatter!(
            df.forecast_time,
            df[:, i],
            mc=:whitesmoke,
            markersize=0.8,
            markerstrokewidth=0,
            label = "")
    end

    # Add historical data
    scatter!(
        df.forecast_time,
        df.avg_actual,
        mc=:brown3,
        markersize=2,
        markerstrokewidth=0,
        label = "Historical")
    
    # Add annotation
    annotate!(DateTime(2018,11,30), 0.95*y_ub,annotation)
    # Add titles
    title!(title)
    xlabel!("Days of 2018")
    ylabel!(type * " [MW]")

    # Filename 
    if dayahead
        name = lowercase(type)*"_"*"day-ahead"*"_historical_landing.png"
    else
        name = lowercase(type)*"_"*"2day-ahead"*"_historical_landing.png"
    end
    savefig(here("plots",name))
end