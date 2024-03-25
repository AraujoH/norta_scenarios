function generate_probability_scenarios(data, scenario_length, number_of_scenarios, q_landing_probability, start_time) #,
    #scenario_hour,
    #scenario_day,
    #scenario_month,
    #scenario_year)

    """
    generate_probability_scenarios: 

    # Arguments:
    - 

    # Returns:
    - 
    """
    # ==================================================================
    # COLUMNS TO KEEP
    # ==================================================================
    #=Here, we care only about the columns in x::DataFrame whose elements 
    are not all equal. If they are, the correlation is not defined b/c
    the standard deviation will be zero for columns whose elements
    are all the same =#
    x = copy(data)
    allequal_set = Set(findall(allequal, eachcol(x)));
    allequal_ind = sort(collect(allequal_set));
    allindex_set = Set(collect(1:48));
    alldifferent_ind = sort(collect(setdiff(allindex_set, allequal_set))); # Index for columns whose st. dev. isn't zero
    x = x[:, alldifferent_ind];


    # ==================================================================
    # CORRELATION MATRIX
    # ==================================================================
    #= Determine a lower-triangular, nonsingular factorization M of the 
        the correlation matrix for Z such that MM' = Sigma_Z. =#
    if ishermitian(cor(x))
        Σ_Z = LinearAlgebra.cholesky(cor(x));
    else
        Σ_Z = factorize(cor(x));
    end
    M = Σ_Z.L;
    dim_M = size(M, 1);

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
        d = Normal(0,1);

        #Generate vector W = (W_1, ..., W_k) ~ iid standard normal
        W = rand(d, scenario_length);

        # Create vector Z such that Z <- MW
        Z = M * W;

        #Compute the CDF of Z
        #cdf_Z = sort(cdf.(d, Z));
        cdf_Z = cdf.(d, Z);
        
        for k in 1:scenario_length
            #Apply the inverse CDF for X_k
            # Y[nscen, k] = quantile(x[:, k], cdf_Z[k])
            Y[nscen, k] = quantile(cdf_Z, q_landing_probability); 
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