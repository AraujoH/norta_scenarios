"""
    transform_landing_probability(data)

Transform landing probability data into a matrix format by filtering and reshaping.

This function processes landing probability data by:
1. Sorting the data by issue time
2. Grouping by issue time and counting occurrences
3. Filtering to keep only groups with exactly 48 entries
4. Reshaping the filtered landing probabilities into a matrix where each row represents
   an issue time and each column represents a time step

# Arguments
- `data`: A DataFrame containing at least the columns `:issue_time` and `:landing_probability`

# Returns
- A transposed matrix where each row corresponds to an issue time (with 48 entries) and 
  contains the landing probability values across time steps

# Example
```julia
transformed_matrix = transform_landing_probability(my_data)
```
"""
function transform_landing_probability(data)
    x = copy(data);
    # Sort data by issue time
    sort!(x, :issue_time);

    # Group data by issue time and count occurences in every group
    df = combine(groupby(x, [:issue_time]), DataFrames.nrow => :count);

    # Filter data by count. Only keep groups with 48 entries
    df_filtered = filter(:count => ==(48), df);
    
    issue_times_interest = df_filtered[!, :issue_time];
    landing_probability_filtered = innerjoin(x, df_filtered, on=:issue_time);
    landing_probability_filtered_matrix = reshape(landing_probability_filtered[!, :landing_probability], (48, size(df_filtered, 1)));
    landing_probability_filtered_matrix = transpose(landing_probability_filtered_matrix);
    return (landing_probability_filtered_matrix)
end