using PlotAndSave
using Plots

function qubit_comparison()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_10"))
    d_10 = plotinfo.lines["|10>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_10"))
    dc_10 = plotinfo.lines["|100>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_10"))
    c200_10 = plotinfo.lines["|100>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_10"))
    pumped_10 = plotinfo.lines["|100>, n1+nc"]

    t_units = 2 * π * 5
    pl = plot(d_10.x / t_units, d_10.y, label = "Direct Coupling", xlabel = "t (μs)", ylabel = "<n1>")
    plot!(dc_10.x / t_units, dc_10.y, label = "Direct with Coupler", xlabel = "t (μs)", ylabel = "<n1+nc>")
    plot!(c200_10.x / t_units, c200_10.y, label = "Coupler200", xlabel = "t (μs)", ylabel = "<n1+nc>")
    plot!(pumped_10.x / t_units, pumped_10.y, label = "Pumped201", xlabel = "t (μs)", ylabel = "<n1+nc>")
    savefig(joinpath(@__DIR__, "qubit_comparison.png"))
    return pl
end

function qubit_comparison_short()
    pl = qubit_comparison()
    plot!(pl, xlims = (0, 1), ylims = (0.95, 1.01))
    savefig(joinpath(@__DIR__, "qubit_comparison_short.png"))
end

function leakage_comparison()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_20"))
    d_20 = plotinfo.lines["|20>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_20"))
    dc_20 = plotinfo.lines["|200>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_20"))
    c200_20 = plotinfo.lines["|200>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_20"))
    pumped_20 = plotinfo.lines["|200>, n1+nc"]

    t_units = 2 * π * 5

    pl = plot(ylims = (0, 2))

    plot!(d_20.x / t_units, d_20.y, label = "Direct Coupling", color = 1)
    plot!(dc_20.x / t_units, dc_20.y, label = "Direct with Coupler", color = 2)
    plot!(c200_20.x / t_units, c200_20.y, label = "Coupler200", color = 3)
    plot!(pumped_20.x / t_units, pumped_20.y, label = "Pumped201", color = 4)

    plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_010"))
    dc_010 = plotinfo.lines["|010>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_010"))
    c200_010 = plotinfo.lines["|010>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_010"))
    pumped_010 = plotinfo.lines["|010>, n1+nc"]

    plot!(dc_010.x / t_units, dc_010.y, label = "Direct with Coupler", linestyle = :dash, color = 2, alpha = 0.66)
    plot!(c200_010.x / t_units, c200_010.y, label = "Coupler200", linestyle = :dash, color = 3, alpha = 0.66)
    plot!(pumped_010.x / t_units, pumped_010.y, label = "Pumped201", linestyle = :dash, color = 4, alpha = 0.66)

    xlabel!("t (μs)")
    ylabel!("<n1+nc>")

    savefig(joinpath(@__DIR__, "leakage_comparison.png"))
end

function leakage_qubit_comparison()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_10_1mus"))
    d_10 = plotinfo.lines["|10>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_10_1mus"))
    dc_10 = plotinfo.lines["|100>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_10_20"))
    c200_10 = plotinfo.lines["|100>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_10_1mus"))
    pumped_10 = plotinfo.lines["|100>, n1+nc"]

    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_20"))
    d_20 = plotinfo.lines["|20>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_20"))
    dc_20 = plotinfo.lines["|200>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_20_longer"))
    c200_20 = plotinfo.lines["|200>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_20"))
    pumped_20 = plotinfo.lines["|200>, n1+nc"]

    function descending(x, y)
        out_x = [x[1]]
        out_y = [y[1]]
        for i in 2:length(x)
            x_i = x[i]
            if x_i < out_x[end]
                push!(out_x, x_i)
                push!(out_y, y[i])
            end
        end
        return out_x, out_y
    end

    d_x, d_y = descending(d_20.y, d_10.y)
    dc_x, dc_y = descending(dc_20.y, dc_10.y)
    c200_x, c200_y = descending(c200_20.y, c200_10.y)
    pumped_x, pumped_y = descending(pumped_20.y, pumped_10.y)

    pl = plot(xlims = (0, 1), ylims = (0.995, 1.01), xflip = true)

    plot!(d_x .- 1, d_y, label = "Direct Coupling L", color = 1)
    plot!(dc_x .- 1, dc_y, label = "Direct with Coupler L", color = 2)
    plot!(c200_x ./ 2, c200_y, label = "Coupler200 L", color = 3)
    plot!(pumped_x .- 1, pumped_y, label = "Pumped201 L", color = 4)


    xlabel!("leakage population")
    ylabel!("qubit population")

    savefig(joinpath(@__DIR__, "leakage_qubit_comparison.png"))
end

function qubit_comparison_10mus()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_10_10mus"))
    d_10 = plotinfo.lines["|10>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_10_10mus"))
    dc_10 = plotinfo.lines["|100>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_10_10mus"))
    c200_10 = plotinfo.lines["|100>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_10_10mus"))
    pumped_10 = plotinfo.lines["|100>, n1+nc"]

    t_units = 2 * π * 5
    xtick_pos = [10^-3, 10^-2, 10^-1, 10^0, 10^1]
    xtick_labels = ["0.001", "0.01", "0.1", "1", "10"]
    pl = plot(xlims = (10^-3, 10^1), ylims = (0.99, 1.01), xaxis = :log10, xticks = (xtick_pos, xtick_labels), size = (1200, 400))
    plot!(pl, d_10.x / t_units, d_10.y, label = "Direct Coupling", xlabel = "t (μs)", ylabel = "<n1>")
    plot!(pl, dc_10.x / t_units, dc_10.y, label = "Direct with Coupler", xlabel = "t (μs)", ylabel = "<n1+nc>")
    plot!(pl, c200_10.x / t_units, c200_10.y, label = "Coupler200", xlabel = "t (μs)", ylabel = "<n1+nc>")
    plot!(pl, pumped_10.x / t_units, pumped_10.y, label = "Pumped201", xlabel = "t (μs)", ylabel = "<n1+nc>")
    savefig(joinpath(@__DIR__, "qubit_comparison_10mus.png"))
    return pl
end

function leakage_comparison_10mus()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_20_10mus"))
    d_20 = plotinfo.lines["|20>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_20_10mus"))
    dc_20 = plotinfo.lines["|200>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_20_10mus"))
    c200_20 = plotinfo.lines["|200>, n1+nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_20_10mus"))
    pumped_20 = plotinfo.lines["|200>, n1+nc"]

    t_units = 2 * π * 5
    xtick_pos = [10^-3, 10^-2, 10^-1, 10^0, 10^1]
    xtick_labels = ["0.001", "0.01", "0.1", "1", "10"]
    pl = plot(ylims = (0, 2), xlims = (10^-3, 10^1), xaxis = :log10, xticks = (xtick_pos, xtick_labels), size = (1200, 400))

    plot!(d_20.x / t_units, d_20.y, label = "Direct Coupling", color = 1)
    plot!(dc_20.x / t_units, dc_20.y, label = "Direct with Coupler", color = 2)
    plot!(c200_20.x / t_units, c200_20.y, label = "Coupler200", color = 3)
    plot!(pumped_20.x / t_units, pumped_20.y, label = "Pumped201", color = 4)

    # plotinfo = loaddata(joinpath(@__DIR__, "DirectWithCoupler", "result_010"))
    # dc_010 = plotinfo.lines["|010>, n1+nc"]
    # plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_010"))
    # c200_010 = plotinfo.lines["|010>, n1+nc"]
    # plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_010"))
    # pumped_010 = plotinfo.lines["|010>, n1+nc"]

    # plot!(dc_010.x / t_units, dc_010.y, label = "Direct with Coupler", linestyle = :dash, color = 2, alpha = 0.66)
    # plot!(c200_010.x / t_units, c200_010.y, label = "Coupler200", linestyle = :dash, color = 3, alpha = 0.66)
    # plot!(pumped_010.x / t_units, pumped_010.y, label = "Pumped201", linestyle = :dash, color = 4, alpha = 0.66)

    xlabel!("t (μs)")
    ylabel!("<n1+nc>")

    savefig(joinpath(@__DIR__, "leakage_comparison_10mus.png"))
end

# qubit_comparison()
# qubit_comparison_short()
# leakage_comparison()
# leakage_qubit_comparison()

qubit_comparison_10mus()
leakage_comparison_10mus()