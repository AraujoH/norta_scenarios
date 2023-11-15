
function plot_cdf_forecast_data(path)
    x = collect(load_forecast[1, 3:101]) / 1000
    p1 = plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation=:left,
        legend=false,
        legendcolumns=3,
        labelfontsize=12,
        background_color_legend=nothing,
        dpi=500)
    plot!([1:99], x, lc=:seagreen, linealpha=0.5)
    scatter!(([1:99], x), ms=1.5, mc=:darkgreen, markerstrokewidth=0)
    ylabel!("Load forecast [GW]")

    y = collect(solar_forecast_dayahead[18, 3:101]) / 1000
    p2 = plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation=:left,
        legend=false,
        legendcolumns=3,
        labelfontsize=12,
        background_color_legend=nothing,
        dpi=500)
    plot!([1:99], y, lc=:darkorange2, linealpha=0.5)
    scatter!(([1:99], y), ms=1.5, mc=:darkorange2, markerstrokewidth=0)
    ylabel!("Solar forecast [GW]")

    z = collect(wind_forecast_dayahead[18, 3:101]) / 1000
    p3 = plot(fontfamily="Computer Modern",
        linewidth=1.0,
        grid=false,
        framestyle=:box,
        titlelocation=:left,
        legend=false,
        legendcolumns=3,
        labelfontsize=12,
        background_color_legend=nothing,
        dpi=500)
    plot!([1:99], z, lc=:brown2, linealpha=0.5)
    scatter!(([1:99], z), ms=1.5, mc=:brown2, markerstrokewidth=0)
    ylabel!("Wind forecast [GW]")

    l = @layout [a b c]
    plot(p1, p2, p3, layout=l)
    savefig(joinpath(path,"some_forecast_data.png"))
end