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

    # Define consistent colors
    colors = [:blue, :orange, :green, :red]

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


qubit_comparison()
leakage_comparison()