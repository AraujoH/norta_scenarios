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