"""
    write_scenarios(x::Matrix{Any}, type::String)

Write the scenarios matrix `x` to a CSV file. The file is saved in a directory
created under `datadir("exp_pro")` with a filename that includes the provided `type`.

# Arguments
- `x::Matrix{Any}`: The matrix containing the scenarios to be written to the CSV file.
- `type::String`: A string that will be included in the filename to distinguish different types of scenarios.

# Example
```julia
write_scenarios(scenarios_matrix, "example_type")
```
"""
function write_scenarios(x::AbstractMatrix, type::String)
    filepath = mkpath(datadir("exp_pro"))
    filename = "scenarios_jul_18_48_" * type * ".csv"
    CSV.write(joinpath(filepath, filename), Tables.table(x), header=false)
end