using PlotAndSave
using Plots

function f()
    transmon_proj = loaddata(joinpath(@__DIR__, "single_transmon"), "proj_data.jld2")
    coupler_proj = loaddata(joinpath(@__DIR__, "coupler_LRU"), "proj_data.jld2")
    transmon_anh = loaddata(joinpath(@__DIR__, "single_transmon"), "anh_data.jld2")
    coupler_anh = loaddata(joinpath(@__DIR__, "coupler_LRU"), "anh_data.jld2")

    t_line = transmon_proj.lines["proj -> 0"]
    c_line = coupler_proj.lines["proj -> 0"]
    x = t_line.x
    pl = plot(x, t_line.y, label = "Single transmon", ylabel = "<Psi|0><0|Psi>", xlabel = "Number of gates")
    plot!(x, c_line.y, label = "Coupler LRU")
    savefig(pl, joinpath(@__DIR__, "proj_comparison.png"))

    t_line = transmon_anh.lines["|2> population"]
    c_line = coupler_anh.lines["Qubit subspace population"]
    x = t_line.x
    pl = plot(x, t_line.y, label = "Single transmon", ylabel = "<Psi|I - anh/2|Psi>", xlabel = "Number of gates")
    plot!(x, c_line.y, label = "Coupler LRU")
    savefig(pl, joinpath(@__DIR__, "anh_comparison.png"))
end

f()