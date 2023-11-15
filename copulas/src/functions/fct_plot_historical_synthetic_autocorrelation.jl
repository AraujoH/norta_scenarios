 function plot_hist_synth_autocor(scenarios, actuals, type, scenario_year, scenario_month, scenario_day, scenario_hour)
    """
    Plot autocorrelation for historical and synthetic data

    Input: 
        - scenarios::DataFrame - synthetic data
        - actuals::DataFrame - historical data
        - type::String - type of data: Load, Solar, or Wind
        - scenario_year::Int64 year for the synthetic data
        - scenario_month::Int64 month for the synthetic data
        - scenario_day::Int64 first day for the synthetic data
        - scenario_hour::Int64 first hour for the synthetic data
    Output:
        - plot (.png) with synthetic and historical autocorrelation
    """
    # Number of scenarios 
    nscen, aux = size(scenarios); 

    # Obtain the historical data
    scenario_timestamp_begin = DateTime(scenario_year, scenario_month, scenario_day, scenario_hour);
    scenario_timestamp_end = scenario_timestamp_begin + Hour(47);

    #= dt returns 48 rows. Every odd-numbered row is a forecast from 24h
    to 48h ahead while every even-numbered row is a forecast up until 23h
    ahead. =#
    dt = filter(:forecast_time => x -> scenario_timestamp_begin <= x <= scenario_timestamp_end, actuals);
    dt = dt[dt.ahead_factor .== "one",:];
    dt_historical = dt[:, :avg_actual];

    # Compute autocorrelation
    lags = 12;    
    cor_hist = Matrix{Float64}(undef,lags,1);
    cor_scen = Matrix{Float64}(undef,nscen,lags);

    cor_hist = autocor(Float64.(dt_historical), range(1,lags));
    for i in 1:nscen
        cor_scen[i,:] =  autocor(hcat(scenarios[i,:]), collect(range(1,lags)))
    end # The use of hcat() is to cast the vector in a type that autocor() takes.
    
    # Initialize plot
    plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation =:left,
        ylim = (-1,1),
        xlim = (0.9,12.1),
        xticks = collect(1:12);
        legend=:outerbottom,
        legendcolumns=2,
        background_color_legend=nothing,
        fg_legend = :transparent,
        dpi=300,
        size=(800,600)) #Removes border from legend

    # Plot scatter for scenario
    x = range(1, lags)
    scatter!(x, cor_scen[1,:], mc =:lightskyblue1, markerstrokewidth=0, msa = .5, markerstrokealpha = .5, label = "Synthetic")
    for i in 2:nscen
        scatter!(x, cor_scen[i,:], mc =:lightskyblue1, markerstrokewidth=0, msa = .5, markerstrokealpha = .5, label=false)
    end

    # Plot scatter for historical
    scatter!(x, cor_hist, mc =:navyblue, legend=:outerbottom, markersize = 6, label = "Historical", markerstrokewidth=0)

    # Add horizontal line at y = 0
    hline!([0], lc=:grey, linestyle=:dash, label = false)

    # Add labels
    xlabel!("Lags")
    ylabel!("Autocorrelation")

    title!(type *" ACF")

    # Save fig
    filename = joinpath(pwd(),"output","figs", type*"_acf.png");
    png(filename);
end