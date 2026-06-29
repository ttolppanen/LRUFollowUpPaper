using PlotAndSave
using Plots
using HDF5
using LaTeXStrings
using Plots.Measures
using LocalPackages

function make_fig()
    plots = []
    annotation_font_size = 28
    annotation_pos = (-0.27, 1.1)

    # ----------------------------------------------------
    plotinfo = loaddata(@__DIR__, "wc_off.jld2")
    display(keys(plotinfo.lines))

    pl = plot(xlabel = L"\Delta_{12} / 2\pi ~(\textrm{MHz})", ylabel = L"\Delta_{1c} / 2\pi~(\textrm{MHz})")

    g_x = plotinfo.lines["g0"].x
    g_y = plotinfo.lines["g0"].y

    zz_x = plotinfo.lines["zz min"].x
    zz_y = plotinfo.lines["zz min"].y

    plot!(pl, g_x .- 4000, 4000 .- g_y, label = L"g_{\textrm{eff}} = 0")
    plot!(pl, zz_x .- 4000, 4000 .- zz_y, label = L"\zeta = 0")

    annotate!(pl, ((-0.24, 1.1), text("(a)", annotation_font_size, "Computer Modern")))
    plot!(pl, topmargin = 16mm, leftmargin = 16mm, bottommargin = 17mm)

    xtick_pos = [-300, -150, 0, 150, 300]
    xtick_labels = ["-300", "-150", "0", "150", "300"]

    plot!(pl, xticks = (xtick_pos, xtick_labels))
    # yflip!(pl)
    push!(plots, pl)

    # ----------------------------------------------------

    file = h5open(joinpath(@__DIR__, "J_effs.h5"), "r")
    J_100 = vec((read(file["J_100"])))
    J_200 = vec((read(file["J_200"])))
    J_201 = vec((read(file["J_201"])))
    close(file)

    wc = 3100:0.01:3575

    zz = [calc_zz(wc_i, [], 4000, 4000, [300, 300, 300], [50, 50, 5]) for wc_i in wc]

    pl = plot()
    plot!(pl, 4000 .- wc, J_100, label = L"g_{\textrm{eff}}", xlabel = L"\Delta_{1c} / 2\pi~(\textrm{MHz})", ylabel = L"J / 2\pi~(\textrm{MHz})", legend = :right)
    hline!(pl, [0], lc=:grey, lw=1, label="")
    plot!(pl, [NaN], [NaN], label = L"\zeta", lc = 2)
    plot!(pl, 4000 .- wc, J_200, label = L"J_{20}", lc = 3)
    plot!(pl, 4000 .- wc, J_201, label = L"J_{21}", lc = 4)
    plot!(twinx(), 4000 .- wc, zz, label = "", ylabel = L"\zeta / 2\pi~(\textrm{MHz})", ylims = (-0.05, 0.0107), lc = 2, legend = :right)

    vline!([495], ls=:dash, lc=1, label="")
    # vline!([-750], ls=:dash, lc=2, label="")
    vline!([793], ls=:dash, lc=4, label="")
    vline!([754], ls=:dash, lc=2, label="")

    plot!(pl, rightmargin = 20mm)
    annotate!(pl, (annotation_pos, text("(b)", annotation_font_size, "Computer Modern")))

    push!(plots, pl)

    # ----------------------------------------------------

    file = h5open(joinpath(@__DIR__, "J_effs_tunneling.h5"), "r")
    J_201 = vec((read(file["J_201"])))
    J_211 = vec((read(file["J_211"])))
    w1 = vec((read(file["w_1"])))

    w1_l = length(w1)
    sl = Int(round(w1_l / 14)) # segment length

    close(file)

    pl = plot( 4000 .- w1, J_201, label = L"J_{201}", xlabel = L"\Delta_{q} / 2\pi~(\textrm{MHz})", ylabel = L"J / 2\pi~(\textrm{MHz})", ylims = (-0.101, 0.101))
    w1_start = w1[1:Int(round(w1_l / 2 - sl / 2))]
    w1_mid = w1[Int(round(w1_l / 2 - sl / 2)):Int(round(w1_l / 2 + sl / 2))]
    w1_end = w1[Int(round(w1_l / 2 + sl / 2)):end]

    J_211_start = J_211[1:Int(round(w1_l / 2 - sl / 2))]
    J_211_mid = J_211[Int(round(w1_l / 2 - sl / 2)):Int(round(w1_l / 2 + sl / 2))]
    J_211_end = J_211[Int(round(w1_l / 2 + sl / 2)):end]

    mid_n = div(length(w1_mid), 2)

    plot!(pl, 4000 .- w1_start, J_211_start, label = "", lc = 2)
    plot!(pl, 4000 .- w1_mid[1:(mid_n)], J_211_mid[1:(mid_n)], label = "", lc = 2, ls = :dash)
    plot!(pl, 4000 .- w1_mid[(mid_n + 2):end], J_211_mid[(mid_n + 2):end], label = "", lc = 2, ls = :dash)
    plot!(pl, 4000 .- w1_end, J_211_end, label = L"J_{211}", lc = 2)

    annotate!(pl, (annotation_pos, text("(c)", annotation_font_size, "Computer Modern")))

    plot!(pl, xticks = (xtick_pos, xtick_labels))
    push!(plots, pl)

    # ----------------------------------------------------

    layout = @layout [a b c]
    pl = plot(plots..., layout = layout, size = 2 .* (3 * 400, 300), margin = 15mm, guidefont = 24, tickfont = 20, legendfont = 20)
    display(pl)
    savefig(pl, joinpath(@__DIR__, "effective_coupling.pdf"))
end

make_fig()