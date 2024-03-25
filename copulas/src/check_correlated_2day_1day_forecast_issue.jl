




correlated_all_hours = load_data[all_12_indices, :]

correlated_active_indices = findall(x -> x .>= start_date, correlated_all_hours[!,:forecast_time])

# active_correlated_hours = correlated_all_hours[correlated_all_hours.forecast_time .>= start_date]

active_correlated_hours = correlated_all_hours[correlated_active_indices, :]

CSV.write(joinpath(filepath,"corr_active_hours.csv"), active_correlated_hours)

quantile_names = names(load_data)[percentile_column_index]

forecast_groups = groupby(correlated_all_hours, :forecast_time)

combine(forecast_groups, quantile_names .=> mean; renamecols=false)
# done!!
### actually, if you want to keep the other columns what do you do?
### do you need to keep issue times?
### XXX also you might need to keep actual?
combine(forecast_groups, "avg_actual" .=> mean; renamecols=false)

# leftjoin()

### XXX but then how do you know the issue time when that becomes important? or is it no longer important
### XXX now that we have determined the active_correlated_hours.... its




