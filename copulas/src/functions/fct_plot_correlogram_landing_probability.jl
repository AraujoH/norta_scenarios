function plot_correlogram_landing_probability(data, type)
    """
    Input:
        - 
    Output:
        -
    """
    # Get correlation matrix
    correl = cor(data);

    # Visualize only the upper triangle matrix
    nrow,ncol = size(correl);

    for i in 1:nrow
        for j in 1:ncol
            if i >= j
                correl[i,j] = NaN                
            end             
        end
    end

    # Prepare plot features
    title = "Landing Probability | " * type * " correlation";
    nticks = collect(1:3:48);


    plot(
        heatmap(
            correl, 
            yflip = true,
            fc = cgrad([:red, :white, :navyblue]),
            clims = (-1,1)        
            ), 
        fontfamily = "Computer Modern",
        title = title;
        xticks = nticks,
        yticks = nticks,
        bordercolor = :white,
        axiscolor = :gainsboro,
        size=(750,750)
        #tickfontcolor = :gainsboro
        #aspect_ratio = :equal
    )

    # Save fig
    filename = joinpath(pwd(),"output","figs", type*"_correlogram.png");
    png(filename)
end