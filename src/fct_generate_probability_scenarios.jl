"""
    generate_probability_scenarios(data, scenario_length, number_of_scenarios)

Generate correlated probability scenarios using NORTA (Normal-to-Anything) method.

Implements the NORTA approach with Cholesky decomposition to generate scenarios that preserve
the correlation structure of the input data. Automatically handles columns with zero variance
by excluding them from correlation matrix computation and re-expanding results.

# Arguments
- `data::DataFrame`: Input data with historical observations (typically 48 columns for hourly data)
- `scenario_length::Int`: Number of time periods per scenario (e.g., 48 for hourly)
- `number_of_scenarios::Int`: Number of scenarios to generate

# Returns
- `Y::Matrix{Float64}`: Generated scenarios with dimensions (number_of_scenarios × scenario_length)

# Details
- Removes columns with constant values (zero standard deviation) before computing correlations
- Uses Cholesky factorization: M such that MM' = Σ (correlation matrix)
- Generates independent standard normal vectors W and transforms via Z = MW
- Applies inverse CDF transformation to match marginal distributions of input data
- Re-expands results to original dimension with zero columns where constants were excluded

# Example
```julia
# scenarios = generate_probability_scenarios(df, 48, 1000)
```

"""
function generate_probability_scenarios(data, scenario_length, number_of_scenarios)

    # ==================================================================
    # COLUMNS TO KEEP
    # ==================================================================
    #=Here, we care only about the columns in x::DataFrame whose elements 
    are not all equal. If they are, the correlation is not defined b/c
    the standard deviation will be zero for columns whose elements
    are all the same =#
    x = copy(data)
    allequal_set = Set(findall(allequal, eachcol(x)))
    allequal_ind = sort(collect(allequal_set))
    allindex_set = Set(collect(1:48))
    alldifferent_ind = sort(collect(setdiff(allindex_set, allequal_set))) # Index for columns whose st. dev. isn't zero
    x = x[:, alldifferent_ind]


    # ==================================================================
    # CORRELATION MATRIX
    # ==================================================================
    #= Determine a lower-triangular, nonsingular factorization M of the 
        the correlation matrix for Z such that MM' = Sigma_Z. =#
    if ishermitian(cor(x))
        Σ_Z = LinearAlgebra.cholesky(cor(x))
    else
        Σ_Z = factorize(cor(x))
    end
    M = Σ_Z.L
    dim_M = size(M, 1)

    # ==================================================================
    # PROBABILITY GENERATION LOOP
    # ==================================================================
    #= In certain cases, as in solar power, not all columns will be 
    useful. Some will be discarded to avoid problems in the factorization
    of the correlation matrix. Here we check if the dimension n of the 
    matrix M (n x n) is equal to the scenario length stipulated as 48.
    If it is not, the vector W will have its length adjusted to n. 
    The variable allequal_ind stores the indices of the columns that 
    were discarded. After the scenarios Y are generated, Y column dimension
    will be expanded and all-zero columns will be introduced in the 
    location of the allequal_ind
    =#

    #Random.seed!(29031990)
    Random.seed!(12345)
    Y = Matrix{Float64}(undef, number_of_scenarios, scenario_length)

    need_expansion = 0 # This is specially important for solar data with several columns whose st. dev. is zero
    if dim_M != scenario_length
        original_scen_length = scenario_length
        scenario_length = dim_M
        Y = Matrix{Float64}(undef, number_of_scenarios, scenario_length)
        need_expansion = 1
    end

    for nscen in 1:number_of_scenarios
        # Set up normal distribution with mean 0 and sd equal to 1
        d = Normal(0, 1)

        #Generate vector W = (W_1, ..., W_k) ~ iid standard normal
        W = rand(d, scenario_length)

        # Create vector Z such that Z <- MW
        Z = M * W

        #Compute the CDF of Z
        #cdf_Z = sort(cdf.(d, Z));
        cdf_Z = cdf.(d, Z)

        for k in 1:scenario_length
            #Apply the inverse CDF for X_k
            # Y[nscen, k] = quantile(x[:, k], cdf_Z[k])
            Y[nscen, k] = quantile(x[:, k], cdf_Z[k])
        end
    end

    #= If there is the need for expansion, expansion is done in the 
    following block of code.
    =#
    if need_expansion == 1
        Y_aux = Matrix{}(undef, number_of_scenarios, original_scen_length)
        Y_aux[:, allequal_ind] .= 0
        Y_aux[:, alldifferent_ind] .= Y
        return Y_aux
    else
        return Y
    end
end


"""
Generate probability scenarios using a NORTA-type construction.

Arguments
---------
data : AbstractMatrix
    Historical data matrix with observations in rows and time periods in columns.
scenario_length : Int
    Number of time periods (must match number of columns in `data`).
number_of_scenarios : Int
    Number of scenarios to generate.
intraday_market_scenarios : Bool
    If true, scenarios share a growing prefix across scenarios.

Returns
-------
Matrix{Float64}
    Scenario matrix of size (number_of_scenarios × scenario_length).
"""
# function generate_probability_scenarios(
#     data::AbstractMatrix{<:Real},
#     scenario_length::Int,
#     number_of_scenarios::Int;
#     intraday_market_scenarios::Bool
# )

#     # ---------------------------------------------------------------
#     # INITIAL SETUP
#     # ---------------------------------------------------------------

#     Random.seed!(12345)                # reproducibility
#     d = Normal(0, 1)                   # standard normal for NORTA
#     T0 = size(data, 2)                 # original number of columns

#     # Sanity check
#     if scenario_length != T0
#         error("scenario_length must equal number of columns in data.")
#     end

#     # Ensure we work with a concrete numeric matrix
#     x0 = Matrix{Float64}(data)

#     # ---------------------------------------------------------------
#     # DROP CONSTANT COLUMNS (ZERO VARIANCE)
#     # ---------------------------------------------------------------
#     # Columns with all equal values have zero std dev and
#     # break correlation and Cholesky factorization.

#     keep_mask = [!allequal(view(x0, :, j)) for j in axes(x0, 2)]
#     keep_inds = findall(keep_mask)          # columns kept
#     drop_inds = findall(.!keep_mask)        # columns dropped

#     x = x0[:, keep_inds]                    # reduced data matrix
#     T = size(x, 2)                          # effective dimension

#     if T == 0
#         error("All columns are constant; cannot build correlation matrix.")
#     end

#     # ---------------------------------------------------------------
#     # CORRELATION MATRIX AND FACTORIZATION
#     # ---------------------------------------------------------------
#     # Correlation matrices are symmetric by construction.
#     # We wrap with Symmetric to enforce this mathematically
#     # and avoid numerical symmetry issues.

#     C = Symmetric(cor(x))

#     # Cholesky factorization may fail due to near-singularity.
#     # We add a tiny diagonal jitter only if needed.
#     F = try
#         cholesky(C)
#     catch
#         cholesky(C + 1e-10 * I)
#     end

#     M = F.L                               # lower-triangular factor

#     # ---------------------------------------------------------------
#     # ALLOCATE ARRAYS FOR SCENARIO GENERATION
#     # ---------------------------------------------------------------

#     Y = Matrix{Float64}(undef, number_of_scenarios, T)   # final scenarios (reduced)
#     W = Matrix{Float64}(undef, T, number_of_scenarios)   # latent Gaussian draws
#     Z = Matrix{Float64}(undef, T, number_of_scenarios)   # correlated Gaussian vars

#     # ---------------------------------------------------------------
#     # SCENARIO GENERATION
#     # ---------------------------------------------------------------
#     # Two modes:
#     # 1) Intraday: scenarios share a growing prefix
#     # 2) Independent: each scenario is fully independent

#     if intraday_market_scenarios
#         for nscen in 1:number_of_scenarios

#             # Fresh iid Gaussian draw
#             W[:, nscen] .= rand(d, T)

#             # Enforce growing prefix across scenarios
#             if nscen > 1
#                 p = min(nscen - 1, T)
#                 W[1:p, nscen] .= W[1:p, nscen - 1]
#             end

#             # Correlate Gaussian variables
#             Z[:, nscen] .= M * W[:, nscen]

#             # Map to uniform via Gaussian CDF
#             cdf_Z = cdf.(d, Z[:, nscen])

#             # Apply inverse empirical CDF marginally
#             for k in 1:T
#                 Y[nscen, k] = quantile(view(x, :, k), cdf_Z[k])
#             end
#         end
#     else
#         for nscen in 1:number_of_scenarios
#             W[:, nscen] .= rand(d, T)
#             Z[:, nscen] .= M * W[:, nscen]
#             cdf_Z = cdf.(d, Z[:, nscen])

#             for k in 1:T
#                 Y[nscen, k] = quantile(view(x, :, k), cdf_Z[k])
#             end
#         end
#     end

#     # ---------------------------------------------------------------
#     # EXPAND BACK TO ORIGINAL DIMENSION IF NEEDED
#     # ---------------------------------------------------------------
#     # If some columns were dropped (e.g. solar night hours),
#     # we reinsert them as all-zero columns in the original positions.

#     if T == T0
#         return Y
#     else
#         Y_aux = Matrix{Float64}(undef, number_of_scenarios, T0)
#         Y_aux[:, drop_inds] .= 0.0
#         Y_aux[:, keep_inds] .= Y
#         return Y_aux
#     end
# end


"""
    generate_probability_scenarios(
        data::AbstractMatrix{<:Real},
        scenario_length::Int,
        number_of_scenarios::Int;
        intraday_market_scenarios::Bool,
        num_batches::Int = 10,
        seed::Int = 12345
    )

Generate correlated probabilistic scenarios from historical data using a
Gaussian copula and empirical marginal distributions.

This function constructs scenarios by:
1. Estimating the empirical correlation structure of the input data,
2. Sampling correlated standard normal vectors,
3. Mapping them to empirical marginals via the probability integral transform.

Columns of `data` are interpreted as time steps within a scenario, and rows
correspond to historical realizations.

Constant columns (zero variance) are automatically removed before correlation
estimation and reinserted as zeros in the output.

# Arguments
- `data`: Matrix of historical observations with size `(N, T)`, where `N` is
  the number of samples and `T` is the scenario length.
- `scenario_length`: Number of time steps per scenario. Must equal `size(data, 2)`.
- `number_of_scenarios`: Number of scenarios to generate (ignored when
  `intraday_market_scenarios = true`).

# Keyword Arguments
- `intraday_market_scenarios`: If `true`, generates `scenario_length` scenarios
  per batch with a rolling-dependence structure suitable for intraday markets.
- `num_batches`: Number of independent batches to generate when
  `intraday_market_scenarios = true`. Default is 10.
- `seed`: Master random seed ensuring reproducibility across batches.

# Behavior
- When `intraday_market_scenarios = false`, the function returns a matrix of
  size `(number_of_scenarios, scenario_length)`.
- When `intraday_market_scenarios = true`, the function returns a matrix of
  size `(num_batches * scenario_length, scenario_length)`, stacking batches
  vertically.

In intraday mode, scenarios are generated sequentially such that earlier
time steps are shared across later scenarios, inducing temporal consistency
within each batch.

# Returns
A matrix of generated scenarios with the same column dimension as the input
data. Columns corresponding to constant input data are filled with zeros.

# Errors
Throws an error if:
- `scenario_length` does not match the number of columns in `data`,
- all columns in `data` are constant.

# Notes
- Correlation factorization is stabilized by adding a small diagonal jitter
  if the empirical correlation matrix is not numerically positive definite.
- A single master RNG is used to derive deterministic per-batch seeds.

"""
function generate_probability_scenarios(
    data::AbstractMatrix{<:Real},
    scenario_length::Int,
    number_of_scenarios::Int;
    intraday_market_scenarios::Bool,
    num_batches::Int = 10,
    seed::Int = 12345
)

    d = Normal(0, 1)

    # Master RNG: the single source of randomness for reproducibility.
    # We only use it to generate per batch seeds.
    rng_master = MersenneTwister(seed)

    # ==============================================================
    # PREP: drop constant columns (zero variance) so cor and cholesky work
    # ==============================================================
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

    # ==============================================================
    # CORRELATION MATRIX AND FACTORIZATION
    # ==============================================================
    C = Symmetric(cor(x))

    F = try
        cholesky(C)
    catch
        cholesky(C + 1e-10 * I)
        @info "Input matrix slightly perturbed."
    end

    M = F.L

    # ---------------------------------------------------------------
    # Decide how many scenarios we generate, and how big the output is
    # ---------------------------------------------------------------
    scenarios_per_batch = intraday_market_scenarios ? scenario_length : number_of_scenarios

    if intraday_market_scenarios && number_of_scenarios != scenario_length
        @warn "intraday_market_scenarios=true ignores number_of_scenarios; using scenario_length=$scenario_length per batch"
    end

    total_rows = intraday_market_scenarios ? (num_batches * scenarios_per_batch) : number_of_scenarios

    # Allocate the big output container for ALL batches, stacked vertically.
    Y_all_reduced = Matrix{Float64}(undef, total_rows, T)

    row0 = 0

    # ---------------------------------------------------------------
    # Batch loop
    # ---------------------------------------------------------------
    for b in 1:(intraday_market_scenarios ? num_batches : 1)

        # Deterministic per batch seed from master seed
        batch_seed = rand(rng_master, UInt)
        rng_batch = MersenneTwister(batch_seed)

        # Allocate arrays for THIS batch only
        Y = Matrix{Float64}(undef, scenarios_per_batch, T)
        W = Matrix{Float64}(undef, T, scenarios_per_batch)
        Z = Matrix{Float64}(undef, T, scenarios_per_batch)

        if intraday_market_scenarios
            for nscen in 1:scenarios_per_batch
                W[:, nscen] .= rand(rng_batch, d, T)

                if nscen > 1
                    p = min(nscen - 1, T)
                    W[1:p, nscen] .= W[1:p, nscen - 1]
                end

                Z[:, nscen] .= M * W[:, nscen]
                cdf_Z = cdf.(d, Z[:, nscen])

                for k in 1:T
                    Y[nscen, k] = quantile(view(x, :, k), cdf_Z[k])
                end
            end

            Y_all_reduced[row0 + 1 : row0 + scenarios_per_batch, :] .= Y
            row0 += scenarios_per_batch

        else
            for nscen in 1:scenarios_per_batch
                W[:, nscen] .= rand(rng_batch, d, T)
                Z[:, nscen] .= M * W[:, nscen]
                cdf_Z = cdf.(d, Z[:, nscen])

                for k in 1:T
                    Y[nscen, k] = quantile(view(x, :, k), cdf_Z[k])
                end
            end

            # Non intraday mode returns just one block
            return (T == T0) ? Y : begin
                Y_aux = Matrix{Float64}(undef, number_of_scenarios, T0)
                Y_aux[:, drop_inds] .= 0.0
                Y_aux[:, keep_inds] .= Y
                Y_aux
            end
        end
    end

    # ==============================================================
    # EXPAND BACK TO ORIGINAL DIMENSION IF WE DROPPED ANY COLUMNS
    # ==============================================================
    if T == T0
        return Y_all_reduced
    else
        Y_all = Matrix{Float64}(undef, total_rows, T0)
        Y_all[:, drop_inds] .= 0.0
        Y_all[:, keep_inds] .= Y_all_reduced
        return Y_all
    end
end


# HSD36 on Jan. 18, 2026 ===========================================
"""
    generate_probability_scenarios_cube!(
        data::AbstractMatrix{<:Real},
        scenario_length::Int,
        number_of_scenarios::Int,
        number_of_iterations::Int,
        number_of_sheets::Int;
        seed::Int = 12345,
        save_path_4d::AbstractString,
        save_path_mean::AbstractString,
    )

Generate a 4D array of correlated probabilistic scenarios with multiple iterations
and sheets, using a Gaussian copula and empirical marginal distributions.

This function extends the NORTA scenario generation approach to produce a large
ensemble of scenarios organized in a 4D structure suitable for advanced stochastic
analysis. Each iteration contains multiple sheets, where each sheet contains a set
of scenarios with a growing-prefix (intraday-style) dependency structure.

The function automatically handles constant columns (zero variance) and saves both
the complete 4D array and iteration-wise mean scenarios to disk using Julia's
serialization format.

# Arguments
- `data`: Matrix of historical observations with size `(N, T)`, where `N` is
  the number of samples and `T` is the scenario length.
- `scenario_length`: Number of time steps per scenario. Must equal `size(data, 2)`.
- `number_of_scenarios`: Number of scenarios to generate per sheet.
- `number_of_iterations`: Number of independent iterations (outermost dimension).
- `number_of_sheets`: Number of sheets per iteration (second dimension).

# Keyword Arguments
- `seed`: Master random seed ensuring reproducibility across all iterations and sheets.
  Default is 12345.
- `save_path_4d`: File path where the 4D scenario array will be serialized.
- `save_path_mean`: File path where the mean scenario array will be serialized.

# Output Structure
The function generates and saves two arrays:

1. **scenarios_4d**: 4D array with dimensions `(number_of_iterations, number_of_sheets, 
   number_of_scenarios, scenario_length)`. Each element `[i, s, n, t]` represents the
   value at time `t` for scenario `n` in sheet `s` of iteration `i`.

2. **scenarios_mean**: 3D array with dimensions `(number_of_iterations, 
   number_of_scenarios, scenario_length)`. Contains the mean across all sheets for
   each iteration, computed as `mean(scenarios_4d[i, :, :, :])` over the sheet dimension.

# Returns
Returns the `scenarios_mean` array (3D) with dimensions `(number_of_iterations, 
number_of_scenarios, scenario_length)`.

# Behavior
- Uses a hierarchical seeding structure: the master RNG generates deterministic
  per-sheet seeds, ensuring full reproducibility.
- Within each sheet, scenarios are generated with a rolling-dependence structure
  where earlier time steps are shared across later scenarios (intraday market style).
- Constant columns in the input data are automatically removed before correlation
  estimation and reinserted as zeros in the output.
- Correlation factorization is stabilized by adding a small diagonal jitter if needed.
- Both output arrays are serialized to disk and can be loaded with `deserialize()`.

# Errors
Throws an error if:
- `scenario_length` does not match the number of columns in `data`,
- all columns in `data` are constant.

# Example
```julia
#= Load historical data
data = Matrix(CSV.read("historical_data.csv", DataFrame))

# Generate 10 iterations × 50 sheets × 48 scenarios × 48 time steps
mean_scenarios = generate_probability_scenarios_cube!(
    data, 48, 48, 10, 50;
    seed = 12345,
    save_path_4d = "scenarios_4d.jls",
    save_path_mean = "scenarios_mean.jls"
)

# Load saved data later
scenarios_4d = deserialize("scenarios_4d.jls")
scenarios_mean = deserialize("scenarios_mean.jls")
=#
```
"""
function generate_probability_scenarios_cube!(
    data::AbstractMatrix{<:Real},
    scenario_length::Int,
    number_of_scenarios::Int,
    number_of_iterations::Int,
    number_of_sheets::Int;
    seed::Int = 12345,
    save_path_4d::AbstractString,
    save_path_mean::AbstractString,
)
    d = Normal(0, 1)

    rng_master = MersenneTwister(seed)

    # -----------------------------
    # Prep: drop constant columns
    # -----------------------------
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

    # -----------------------------
    # Correlation and factorization
    # -----------------------------
    C = Symmetric(cor(x))
    F = try
        cholesky(C)
    catch
        cholesky(C + 1e-10 * I)
        @info "Input matrix slightly perturbed."
    end
    M = F.L

    # -----------------------------
    # Allocate outputs (full size)
    # scenarios_4d[i, s, n, t]
    # scenarios_mean[i, n, t]
    # -----------------------------
    Nit = number_of_iterations
    nsheets = number_of_sheets
    Nscen = number_of_scenarios

    scenarios_4d = Array{Float64}(undef, Nit, nsheets, Nscen, T0)
    scenarios_mean = Array{Float64}(undef, Nit, Nscen, T0)

    # -----------------------------
    # Work buffers (reused)
    # Reduced-dimension buffers operate on kept columns only
    # -----------------------------
    W = Matrix{Float64}(undef, T, Nscen)
    Z = Matrix{Float64}(undef, T, Nscen)
    Y_red = Matrix{Float64}(undef, Nscen, T)
    Y_full = Matrix{Float64}(undef, Nscen, T0)
    cdf_Z = Vector{Float64}(undef, T)

    for i in 1:Nit
        # mean over sheets for this iteration (full dimension)
        @views scenarios_mean[i, :, :] .= 0.0

        for s in 1:nsheets
            # deterministic per (iteration, sheet) seed from master seed
            sheet_seed = rand(rng_master, Int)
            rng = MersenneTwister(sheet_seed)

            # build reduced scenarios for this sheet
            for n in 1:Nscen
                @views W[:, n] .= rand(rng, d, T)

                # remember the past across scenario paths
                if n > 1
                    p = min(n - 1, T)
                    @views W[1:p, n] .= W[1:p, n - 1]
                end

                @views Z[:, n] .= M * W[:, n]

                @inbounds for t in 1:T
                    cdf_Z[t] = cdf(d, Z[t, n])
                end

                @inbounds for t in 1:T
                    Y_red[n, t] = quantile(view(x, :, t), cdf_Z[t])
                end
            end

            # expand back to original dimension (put zeros in dropped columns)
            if T == T0
                Y_full .= Y_red
            else
                Y_full[:, drop_inds] .= 0.0
                Y_full[:, keep_inds] .= Y_red
            end

            # store sheet block and accumulate mean
            @views scenarios_4d[i, s, :, :] .= Y_full
            @views scenarios_mean[i, :, :] .+= Y_full
        end

        @views scenarios_mean[i, :, :] ./= nsheets
    end

    serialize(save_path_4d, scenarios_4d)
    serialize(save_path_mean, scenarios_mean)

    return scenarios_mean
end
