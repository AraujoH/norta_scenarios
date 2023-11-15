function plot_correlation_heatmap(data)
    """
    compute_landing_probability:
    
    # Arguments
    
    # Returns
    
    """
    # Convert Julia DataFrame to R DataFrame
    @rput data
    
    # Load the required R packages
    R"""
    library(corrplot)
    
    # Create a correlation matrix in R
    correlation_matrix <- cor(data)
    
    # Plot the correlation matrix using corrplot in R
    corplot = corrplot::corrplot(correlation_matrix, method = "color")
    """
    @rget corrplot
    return corplot
end