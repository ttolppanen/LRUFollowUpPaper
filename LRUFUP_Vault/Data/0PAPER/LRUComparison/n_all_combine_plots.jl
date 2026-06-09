using PlotAndSave
using Plots

function qubit_comparison_10mus()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_10_10mus"))
    d_10 = plotinfo.lines["|10>, n_all"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_10_10mus"))
    c200_10 = plotinfo.lines["|100>, n_all"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_10_10mus"))
    pumped_10 = plotinfo.lines["|100>, n_all"]

    t_units = 2 * π * 5
    xtick_pos = [10^-3, 10^-2, 10^-1, 10^0, 10^1]
    xtick_labels = ["0.001", "0.01", "0.1", "1", "10"]
    pl = plot(xlims = (10^-3, 10^1), ylims = (0.5, 2), xaxis = :log10, xticks = (xtick_pos, xtick_labels), size = (600, 400))
    plot!(pl, d_10.x / t_units, d_10.y, label = "Direct Coupling", xlabel = "t (μs)")
    plot!(pl, c200_10.x / t_units, c200_10.y, label = "Coupler200", xlabel = "t (μs)")
    plot!(pl, pumped_10.x / t_units, pumped_10.y, label = "Pumped201", xlabel = "t (μs)")
    
    ylabel!("<n_all>")
    savefig(joinpath(@__DIR__, "n_all_qubit_comparison_10mus.png"))
    return pl
end

function leakage_comparison_10mus()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_20_10mus"))
    d_20 = plotinfo.lines["|20>, n_all"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_20_10mus"))
    c200_20 = plotinfo.lines["|200>, n_all"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_20_10mus"))
    pumped_20 = plotinfo.lines["|200>, n_all"]

    t_units = 2 * π * 5
    xtick_pos = [10^-3, 10^-2, 10^-1, 10^0, 10^1]
    xtick_labels = ["0.001", "0.01", "0.1", "1", "10"]
    pl = plot(ylims = (0, 3), xlims = (10^-3, 10^1), xaxis = :log10, xticks = (xtick_pos, xtick_labels), size = (600, 400))

    plot!(d_20.x / t_units, d_20.y, label = "Direct Coupling", color = 1)
    plot!(c200_20.x / t_units, c200_20.y, label = "Coupler200", color = 2)
    plot!(pumped_20.x / t_units, pumped_20.y, label = "Pumped201", color = 3)

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
    ylabel!("<n_all>")

    savefig(joinpath(@__DIR__, "n_all_leakage_comparison_10mus.png"))
end


qubit_comparison_10mus()
leakage_comparison_10mus()
