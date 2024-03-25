<<<<<<< Updated upstream:copulas/src/test.jl
println("at least this works")
=======
findall(x -> x == active_issues[1], solar_forecast_2dayahead[!,:issue_time]) 

load_data[all_hours,:]


all_12_indices = findall(x -> x in  active_issues, load_data[!,:issue_time])
all_hours = load_data[all_12_indices, :forecast_time]

active_hours = all_hours[all_hours .>= start_date]

load_data[all_hours,:]

unique(load_data[all_hours,:forecast_time])

unique_active_hours = unique(active_hours)
# XXX so it goes from an 82 element vector to 58, what is going on?

# convert DateTime format to something readable...
Dates.format.(active_hours, "yyyy-mm-dd HH:MM:SS")
Dates.format.(unique_active_hours, "yyyy-mm-dd HH:MM:SS")
# print out these into csv files to check them out...
filepath = mkpath(datadir("test"))
CSV.write(joinpath(filepath,"active_hours.csv"), Tables.table(active_hours))
CSV.write(joinpath(filepath,"unique_active_hours.csv"), Tables.table(unique_active_hours))


function nonunique2(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatedvector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], xs[i])))
            push!(duplicatedvector,xs[i])
        end
    end
    duplicatedvector
end


# check the nonunique between the two
setdiff(unique_active_hours, active_hours) # doesn't work becasue the ones are repeating!
# use the above function instead
nonunique2(active_hours)
# oh but we already knew this
# this could be used for the whole 8764 to check for any duplicates through


# actually it is not weird that the 58 + 24 would give 82, there is overlap of the forecasts. 
### does the 2 day ahead include the exact same forecast for day ahead... or no?
# choose a forecast and issue time to test

# extract the probability quantiles from the day-ahead and the 2day-ahead

# boolean for TRUE or FALSE

#= this tells us if all we need to do is use the 2day ahead because it already has the day ahead data =#

### also try extract 2day ahead from solar at active issues to see its length

solar_12_indices = findall(x -> x in active_issues, solar_forecast_2dayahead[!,:issue_time])
# no it is only 48... i think this forecast skips a day
# and the reason the load has more is because it combines day ahead and 2day ahead
findall(x -> x in nonunique2(load_times[!,:forecast_time]), load_times[!,:forecast_time])
# note that load_data has 17528 rows... why the 8 more than 17520 which is what some raw data has


percentile_column_index = startswith.(names(load_data), "p_")
percentile_names = 

groupby(load_data, :forecast_time)

groupby()
>>>>>>> Stashed changes:scripts/test.jl
