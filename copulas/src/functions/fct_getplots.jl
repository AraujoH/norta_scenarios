function getplots(data, title, type)
    """
    getplots: print two graphs - scatterplot and histogram - for the 
    variance of the percentile forecasts. 

    # Arguments:
    - data::DataFrame: load_data, or solar_data, or wind_data
    - title::String: string with the name for the plot title

    # Output
    Two graphs: scatterplot for the variance, histogram for the variance.
    """
    # Select only percentile columns
    x = copy(data);
    p_cols = names(x, Regex("^p_"));
    x_p = x[:, p_cols];

    # Compute standard deviation
    stdev = zeros(size(x_p, 1))
    for i in 1:size(x_p, 1)
        stdev[i] = std(x_p[i, :])
    end

    # Select only non-zero standard deviation
    stdev_nonzero = stdev[stdev .!= 0, :]
    max_sd = maximum(stdev_nonzero)
    #annotation = string(length(stdev_nonzero)) * " data points"
    range_x = range(1, size(stdev_nonzero, 1))

    @rput stdev_nonzero
    R"""
    library(ggplot2)
    library(extrafont)
    font_import(pattern = "lmroman*") 
    loadfonts()

    y = stdev_nonzero
    y_len = length(y)
    subtitle = paste0(y_len, " data points")
    # Create dataframe
    
    df = data.frame(x = (1:nrow(y)), y = y)

    # Scatter plot
    ggplot(
        df, 
        aes(x, y)
    ) + 
    geom_point(color = "deepskyblue3") + 
    labs(x = "Data points", 
         y = "Standard deviation [MW]",
         title = $title,
         subtitle = subtitle)+
    theme_bw() + 
    theme(text = element_text(family = "LM Roman 10"))
    """

    # # Scatterplot
    # plot(fontfamily="Computer Modern", grid=false, framestyle=:box, legend=false, xformatter = :plain);
    # scatter!(range_x, stdev_nonzero, ylimits = (0, max_sd), ms=0.7, ma=0.9, markerstrokewidth=0, mc= "deepskyblue3");
    # ylabel!("Standard deviation [MW]");
    # xlabel!("Number of data points");
    # title!(title, title_pos= :left);
    # #annotate!([],[]text((annotation),:bottom, :right, 15));
    # savefig("fig_" * string(title) * ".png")

    # # Histogram
    # plot(fontfamily="Computer Modern", grid=false, framestyle=:box, legend=false)
    # histogram!(stdev_nonzero, color="deepskyblue3")
    # title!("Distribution of standard deviation for " * title, title_pos=:left)
    # xlabel!("Standard deviation [MW]")
    # ylabel!("Frequency")
    # savefig("fig_" * string(title) * "_hist.png")
end;