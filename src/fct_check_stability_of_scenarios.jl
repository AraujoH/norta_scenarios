"""
    check_stability_of_scenarios_mean(scenarios_4d, i;
        tol_rel = 1e-3,
        k_consecutive = 10,
        norm_kind = :fro
    )

Check whether the mean of scenario sheets stabilizes for a given iteration.

For a fixed iteration `i`, this function computes the running mean of the
`48×48` scenario matrices across sheets and monitors the relative change
between successive running means. Stability is declared when the relative
change stays below `tol_rel` for `k_consecutive` sheets in a row.

Returns:
- `(s, deltas)` if stability is reached at sheet `s`
- `(nothing, deltas)` otherwise

`deltas[s]` stores the relative change in the mean when adding sheet `s`.
"""
function check_stability_of_scenarios_mean(scenarios_4d, i;
    tol_rel::Float64 = 1e-3,
    k_consecutive::Int = 10,
)
    nsheets = size(scenarios_4d, 2)
    N       = size(scenarios_4d, 3)
    T       = size(scenarios_4d, 4)

    mean_prev = zeros(Float64, N, T)
    mean_curr = zeros(Float64, N, T)

    consec = 0
    deltas = fill(Inf, nsheets)

    for s in 1:nsheets
        Y = @view scenarios_4d[i, s, :, :]

        if s == 1
            mean_curr .= Y
            continue
        end

        mean_prev .= mean_curr
        @. mean_curr = mean_prev + (Y - mean_prev) / s

        Δ = mean_curr .- mean_prev
        num = opnorm(Δ, 2)                       # spectral norm
        den = max(opnorm(mean_curr, 2), eps())   # spectral norm
        deltas[s] = num / den

        if deltas[s] < tol_rel
            consec += 1
            if consec >= k_consecutive
                return s, deltas
            end
        else
            consec = 0
        end
    end

    return nothing, deltas
end
