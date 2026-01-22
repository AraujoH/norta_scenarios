

"""
    generate_probability_IDM_scenarios_cube!(
        intraday_hour::Int,
        data::AbstractMatrix{<:Real},
        historical_normalVariate_4d_filepath::String,
        scenario_length::Int,
        number_of_scenarios::Int,
        number_of_iterations::Int,
        number_of_sheets::Int;
        seed::Int = 12345,
        save_path_scenarios_4d::AbstractString,
        save_path_W_4d::AbstractString,
    )

Generate 4D intraday market scenarios conditioned on realized values up to a
specified hour, using precomputed historical normal variates.

This function extends NORTA scenario generation for intraday markets by fixing
the first `intraday_hour` time steps using historical normal variates from a
previous run, then generating new scenarios for the remaining hours. This
conditioning reflects the intraday market context where early hours are already
realized and only future hours need stochastic simulation.

# Arguments
- `intraday_hour`: Hour index at which intraday market scenarios start (1-based).
  Hours 1 to `intraday_hour` are conditioned on historical variates.
- `data`: Historical observations matrix `(N, T)` where `N` is sample size and
  `T` is scenario length.
- `historical_normalVariate_4d_filepath`: Path to serialized 4D array of normal
  variates from non-intraday hours (from previous scenario generation).
- `scenario_length`: Number of time steps per scenario. Must equal `size(data, 2)`.
- `number_of_scenarios`: Number of scenarios per sheet.
- `number_of_iterations`: Number of independent iterations (outermost dimension).
- `number_of_sheets`: Number of sheets per iteration.

# Keyword Arguments
- `seed`: Master random seed for reproducibility. Default is 12345.
- `save_path_scenarios_4d`: File path for serializing the 4D scenario array.
- `save_path_W_4d`: File path for serializing the 4D normal variate array.

# Returns
- `scenarios_4d`: 4D array `(number_of_iterations, number_of_sheets, 
  number_of_scenarios, scenario_length)` containing generated scenarios.
- `W_4d`: 4D array `(number_of_iterations, number_of_sheets, scenario_length,
  number_of_scenarios)` containing all normal variates used in generation.

# Behavior
- Loads historical normal variates and extracts diagonal elements for hours
  1 to `intraday_hour`, conditioning scenarios on these realized values.
- Generates fresh random normal variates for hours `intraday_hour + 1` to `scenario_length`.
- Applies Gaussian copula transformation with empirical marginals from `data`.
- Automatically handles constant columns (zero variance) by excluding from
  correlation computation and filling with zeros in output.
- Serializes both output arrays to disk using Julia's `serialize()`.

# Example
```julia
# # Generate base scenarios and save normal variates
# generate_probability_scenarios_cube!(...)

# # Later, generate intraday scenarios conditioned on hour 12
# scenarios_4d, W_4d = generate_probability_IDM_scenarios_cube!(
#     12, data, "historical_W_4d.jls", 48, 48, 10, 50;
#     save_path_scenarios_4d = "idm_scenarios_4d.jls",
#     save_path_W_4d = "idm_W_4d.jls"
# )
```

# See Also
- `generate_probability_scenarios_cube!`: Base function for unconditional scenarios.
- `deserialize`: Load serialized arrays from disk.

"""
function generate_probability_IDM_scenarios_cube!(
    intraday_hour::Int,
    data::AbstractMatrix{<:Real},
    historical_normalVariate_4d_filepath::String,
    scenario_length::Int,
    number_of_scenarios::Int,
    number_of_iterations::Int,
    number_of_sheets::Int;
    seed::Int=12345,
    save_path_scenarios_4d::AbstractString,
    save_path_W_4d::AbstractString,
)
    # Distributions and randomness..............................................
    d = Normal(0, 1)
    rng_master = MersenneTwister(seed)

    # Prep: drop constant columns. This applies mostly to solar data as 
    # nighttime hours are constant zero.........................................
    x0 = Matrix{Float64}(data)
    T0 = size(x0, 2)

    if scenario_length != T0
        error("scenario_length must equal the number of columns in data.")
    end

    keep_mask = [!allequal(view(x0, :, j)) for j in axes(x0, 2)]
    keep_inds = findall(keep_mask)
    drop_inds = findall(.!keep_mask)

    x = x0[:, keep_inds]
    T = size(x, 2)

    if T == 0
        error("All columns are constant. Cannot generate scenarios.")
    end

    # Correlation and factorization ............................................
    C = Symmetric(cor(x))
    F = try
        cholesky(C)
    catch
        cholesky(C + 1e-10 * I)
        @info "Input matrix slightly perturbed."
    end
    M = F.L

    # Allocate outputs (full size)..............................................
    Nit = number_of_iterations
    nsheets = number_of_sheets
    Nscen = number_of_scenarios

    # Scenarios: full dimension (includes dropped columns filled with zeros)
    # scenarios_4d[i, s, n, t]
    scenarios_4d = Array{Float64}(undef, Nit, nsheets, Nscen, T0)

    # Normal variates: kept dimension only
    # W_4d[i, s, t, n]
    W_4d = Array{Float64}(undef, Nit, nsheets, T, Nscen)

    # Historical variates from the non-intraday hours
    historical_w_4d = deserialize(historical_normalVariate_4d_filepath)

    # Work buffers (reused)
    # W = Matrix{Float64}(undef, T, Nscen)
    Z = Matrix{Float64}(undef, T, Nscen)
    Y_red = Matrix{Float64}(undef, Nscen, T)
    Y_full = Matrix{Float64}(undef, Nscen, T0)
    cdf_Z = Vector{Float64}(undef, T)


    for i in 1:Nit
        for s in 1:nsheets
            # New RNG for this sheet
            sheet_seed = rand(rng_master, Int)
            rng = MersenneTwister(sheet_seed)

            # Fetch the Ws from hour 1 to intraday hour from the 
            # precomputed 4D array
            W_4d[i, s, 1:intraday_hour, :] .= diag(historical_w_4d[i, s, 1:intraday_hour, 1:intraday_hour])

            # Build the rest of W_4d for this (i, s) with new random values. 
            # This time we do not fix the past values.
            W_new = rand(rng, d, (T - intraday_hour, Nscen))

            # Fill in the rest of W_4d
            W_4d[i, s, intraday_hour+1:end, :] .= W_new

            # Build scenarios for this sheet
            for n in 1:Nscen
                @views Z[:, n] .= M * W_4d[i, s, :, n]

                @inbounds for t in 1:T
                    cdf_Z[t] = cdf(d, Z[t, n])
                end

                @inbounds for t in 1:T
                    Y_red[n, t] = quantile(view(x, :, t), cdf_Z[t])
                end
            end

            # Expand back to original dimension
            if T == T0
                Y_full .= Y_red
            else
                Y_full[:, drop_inds] .= 0.0
                Y_full[:, keep_inds] .= Y_red
            end

            # Store scenarios for this sheet
            @views scenarios_4d[i, s, :, :] .= Y_full
        end
    end

    serialize(save_path_W_4d, W_4d)
    serialize(save_path_scenarios_4d, scenarios_4d)

    return scenarios_4d, W_4d
end