using HDF5
using Plots

function f()
    file = h5open(joinpath(@__DIR__, "J_effs.h5"), "r")
    J_100 = vec((read(file["J_100"])))
    J_200 = vec((read(file["J_200"])))
    J_201 = vec((read(file["J_201"])))
    close(file)

    wc = 3100:0.01:3575

    zz = [calc_zz(wc_i, [], 4000, 4000, [300, 300, 300], [50, 50, 5]) for wc_i in wc]

    pl = plot(wc, J_100, label = "|100>", xlabel = "wc (MHz)", ylabel = "J_eff", title = "Effective coupling strength")
    plot!(pl, wc, J_200, label = "|200>", xlabel = "wc (MHz)")
    plot!(pl, wc, J_201, label = "|201>", xlabel = "wc (MHz)")
    plot!(pl, wc, 50 *zz, label = "ZZ * 50", xlabel = "wc (MHz)")

    savefig(joinpath(@__DIR__, "coupling_strengths.png"))
end

f()