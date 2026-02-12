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

function write_scenarios_to_file(
    x::Array{Float64, 3},
    scenario_day::Int,
    scenario_month::Int,
    scenario_year::Int,
    data_type::String;
    filepath::Union{Nothing,String}=nothing,
    is_idm::Bool=false,
    intraday_hour::Union{Nothing,Int}=nothing,
    iteration::Union{Nothing,Int}=nothing
)

    # Reshape x to 2D
    reshaped_x = reduce(vcat, eachslice(x, dims=1))

    # Set output directory
    filepath = filepath === nothing ? mkpath(datadir("exp_pro")) : mkpath(filepath)

    # Build filename
    idm_tag = is_idm ? "IDM_" : ""
    iter_tag = iteration !== nothing ? "iter_$(iteration)_" : ""
    filename =
        data_type * "_" *
        "scenarios_" *
        idm_tag *
        iter_tag *
        string(intraday_hour) * "h_" *
        string(scenario_year) * "_" *
        lpad(scenario_month, 2, '0') * "_" *
        lpad(scenario_day, 2, '0') * "_" *
        ".csv"

    # Write to CSV
    CSV.write(joinpath(filepath, filename), Tables.table(reshaped_x), header=false)
end
