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
    generate_probability_scenarios_cube!(
        data::AbstractMatrix{<:Real},
        scenario_length::Int,
        number_of_scenarios::Int,
        number_of_iterations::Int,
        number_of_sheets::Int;
        seed::Int = 12345,
        save_path_scenarios_4d::AbstractString,
        save_path_W_4d::AbstractString,
    )

Generate a 4D ensemble of correlated scenarios with intraday-style rolling
dependencies using NORTA method and save to disk.

Creates `(number_of_iterations, number_of_sheets, number_of_scenarios, scenario_length)`
scenarios where each sheet uses a growing-prefix structure: later scenarios
share earlier time steps with preceding scenarios.

# Arguments
- `data`: Historical observations `(N, T)` for estimating marginals and correlations
- `scenario_length`: Time steps per scenario (must equal `size(data, 2)`)
- `number_of_scenarios`: Scenarios per sheet
- `number_of_iterations`: Independent iterations (outermost dimension)
- `number_of_sheets`: Sheets per iteration
- `seed`: Master random seed for reproducibility
- `save_path_scenarios_4d`: Output path for scenario array
- `save_path_W_4d`: Output path for normal variate array

# Returns
- `scenarios_4d`: 4D scenario array `[i, s, n, t]`
- `W_4d`: 4D normal variate array `[i, s, t, n]` (kept dimensions only)

Automatically handles constant columns and serializes outputs for later use with
`deserialize()`.
"""
function generate_probability_scenarios_cube!(
    data::AbstractMatrix{<:Real},
    scenario_length::Int,
    number_of_scenarios::Int,
    number_of_iterations::Int,
    number_of_sheets::Int;
    seed::Int = 12345,
    save_path_scenarios_4d::AbstractString,
    save_path_W_4d::AbstractString,
)
    d = Normal(0, 1)
    rng_master = MersenneTwister(seed)

    # Prep: drop constant columns
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

    # Correlation and factorization
    C = Symmetric(cor(x))
    F = try
        cholesky(C)
    catch
        cholesky(C + 1e-10 * I)
        @info "Input matrix slightly perturbed."
    end
    M = F.L

    # Allocate outputs (full size)
    Nit = number_of_iterations
    nsheets = number_of_sheets
    Nscen = number_of_scenarios

    # Scenarios: full dimension (includes dropped columns filled with zeros)
    # scenarios_4d[i, s, n, t]
    scenarios_4d = Array{Float64}(undef, Nit, nsheets, Nscen, T0)

    # Normal variates: kept dimension only
    # W_4d[i, s, t, n]
    W_4d = Array{Float64}(undef, Nit, nsheets, T, Nscen)

    # Work buffers (reused)
    W = Matrix{Float64}(undef, T, Nscen)
    Z = Matrix{Float64}(undef, T, Nscen)
    Y_red = Matrix{Float64}(undef, Nscen, T)
    Y_full = Matrix{Float64}(undef, Nscen, T0)
    cdf_Z = Vector{Float64}(undef, T)

    for i in 1:Nit
        for s in 1:nsheets
            sheet_seed = rand(rng_master, Int)
            rng = MersenneTwister(sheet_seed)

            # Build W for this (i, s)
            for n in 1:Nscen
                @views W[:, n] .= rand(rng, d, T)

                # Remember the past across scenario paths
                if n > 1
                    p = min(n - 1, T)
                    @views W[1:p, n] .= W[1:p, n - 1]
                end
            end

            # Store W for this sheet
            @views W_4d[i, s, :, :] .= W

            # Build scenarios for this sheet
            for n in 1:Nscen
                @views Z[:, n] .= M * W[:, n]

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