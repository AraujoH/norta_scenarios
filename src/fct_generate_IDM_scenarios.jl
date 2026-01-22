number_of_iterations = 10
number_of_sheets = 50
Nscen = number_of_scenarios

@assert number_of_iterations == size(scenario_means, 1) "number_of_iterations must match the first dimension of scenario_means."

data = copy(lp_load)
d = Normal(0, 1)
seed = 29031990
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
n_iter = number_of_iterations
n_sheets = number_of_sheets


w_means = Array{Float64}(undef, n_iter, T0, Nscen)
w_means_transp = Array{Float64}(undef, n_iter, Nscen, T0)

w_means_path = joinpath(pwd(), "results", "w_means.jls")
w_means = deserialize(w_means_path)

for i in 1:n_iter
    w_means_transp[i, :, :] .= permutedims(w_means[i, :, :])
end

intraday_hours = [6]

W = Matrix{Float64}(undef, T, Nscen)
Z = Matrix{Float64}(undef, T, Nscen)
scenarios_idm_4d = Array{Float64}(undef, n_iter, n_sheets, Nscen, T0)
scnearios_idm_mean = Array{Float64}(undef, n_iter, Nscen, T0)

Y_red = Matrix{Float64}(undef, Nscen, T)
Y_full = Matrix{Float64}(undef, Nscen, T0)
cdf_Z = Vector{Float64}(undef, T)

for iter in 1:1 #number_of_iterations
    for hour in intraday_hours
        past_values = reshape(w_means_transp[iter, hour, 1:hour], 1, :)
        
        for sheet in 1:number_of_sheets
            sheet_seed = rand(rng_master, Int)
            rng = MersenneTwister(sheet_seed)

            # Fix the past values into W
            W[:, 1:hour] .= past_values    
            for n in (hour+1):Nscen
                # Generate normal random variables ~ N(0,1)
                W[n , hour+1:end] = rand(rng, d, T - hour)

                # Apply correlation
                @views Z[:, n] = M * W[:, n]

                # Transform to uniform
                for t in hour+1:T
                    cdf_Z[t] = cdf(d, Z[t, n])
                    Y_red[n, t] = quantile(x[iter, hour, t], cdf_Z[t])
                end
            end
            # Expand back to original dimension (put zeros in dropped columns)
            # Apply to solar generation
            if T == T0
                Y_full .= Y_red
            else
                Y_full[:, drop_inds] .= 0.0
                Y_full[:, keep_inds] .= Y_red
            end
        end
    end
end
