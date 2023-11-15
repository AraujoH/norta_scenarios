function compute_landing_probability(data)
    """
    compute_landing_probability:

    # Arguments
    - data::DataFrame: dataframe with the average historical data and the forecast
        data

    # Returns
    - landing_probability::DataFrame: dataframe with the landing 
        probability. DF has 4 columns: issue time, forecast_time, 
        landing_probability, ahead_factor    
    """
    percentile_column_index = startswith.(names(data), "p_");
    landing_probability = Vector{Float64}(undef, size(data, 1));

    # The BaseStats.ecdf() function applied to the forecasts returns a 
    # function. The returned function is then saved as empirical_cdf 
    # and used to calculate the empirical CDF of the historical 
    # averaged actual data.
    #mm = Matrix{}(undef, size(data,1), 4);
    for i in range(1, size(data, 1))
        quantiles = collect(data[i, percentile_column_index]);
        empirical_cdf = ecdf(quantiles);
        landing_probability[i] = empirical_cdf(data[i, :avg_actual]);
        # -----
        #= Storing min and max of adjust empirical_cdf function
        mm[i,3] = empirical_cdf(0); #Min. value emp. CDF can compute
        mm[i,4] = empirical_cdf(1); #Max. value emp. CDF can compute
        =#
    end

    # -----
    #= Adding forecast and issue times to mm
    println(pwd())
    mm[:,1] = data[!, :forecast_time];
    mm[:,2] = data[!, :issue_time];
    writedlm(joinpath(pwd(),"output","max_min_ecdf.txt"), mm, ";")=#

    # Create a DataFrame with the forecast time, issue time, ahead factor
    # and the landing probability.
    landing_probability = DataFrame(landing_probability=landing_probability)
    landing_probability[!, :forecast_time] = data[!, :forecast_time]
    landing_probability[!, :issue_time] = data[!, :issue_time]
    landing_probability[!, :ahead_factor] = data[:, :ahead_factor]
    select!(landing_probability, [:issue_time, :forecast_time, :landing_probability, :ahead_factor])
    return (landing_probability)
end