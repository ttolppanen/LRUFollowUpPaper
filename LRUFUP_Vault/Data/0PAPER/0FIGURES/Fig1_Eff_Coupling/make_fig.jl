using PlotAndSave
using Plots
using HDF5
using LaTeXStrings
using Plots.Measures
using LocalPackages

function make_fig()
    plots = []
    annotation_font_size = 36
    annotation_pos = (-0.25, 1.15)

    # ----------------------------------------------------
    plotinfo = loaddata(@__DIR__, "wc_off.jld2")
    display(keys(plotinfo.lines))

    pl = plot(xlabel = L"\Delta_{21} / 2\pi ~(\textrm{MHz})", ylabel = L"\Delta_{c1} / 2\pi~(\textrm{MHz})")

    g_x = plotinfo.lines["g0"].x
    g_y = plotinfo.lines["g0"].y

    zz_x = plotinfo.lines["zz min"].x
    zz_y = plotinfo.lines["zz min"].y

    plot!(pl, g_x .- 4000, g_y .- 4000, label = L"g_{\textrm{eff}} = 0")
    plot!(pl, zz_x .- 4000, zz_y .- 4000, label = L"\zeta = 0")

    annotate!(pl, (annotation_pos, text("(g)", annotation_font_size, "Computer Modern")))
    plot!(pl, topmargin = 17mm)

    xtick_pos = [-300, -150, 0, 150, 300]
    xtick_labels = ["-300", "-150", "0", "150", "300"]

    plot!(pl, xticks = (xtick_pos, xtick_labels))
    push!(plots, pl)

    # ----------------------------------------------------

    file = h5open(joinpath(@__DIR__, "J_effs.h5"), "r")
    J_100 = vec((read(file["J_100"])))
    J_200 = vec((read(file["J_200"])))
    J_201 = vec((read(file["J_201"])))
    close(file)

    wc = 3100:0.01:3575

    zz = [calc_zz(wc_i, [], 4000, 4000, [300, 300, 300], [50, 50, 5]) for wc_i in wc]

    pl = plot(wc .- 4000, J_100, label = L"g_{\textrm{eff}}", xlabel = L"\Delta_{c1} / 2\pi~(\textrm{MHz})", ylabel = L"J / 2\pi")
    plot!(pl, wc .- 4000, J_200, label = L"J_{20}")
    plot!(pl, wc .- 4000, J_201, label = L"J_{21}")
    plot!(twinx(), wc .- 4000, zz, label = L"\zeta", ylabel = L"\zeta / 2\pi", ylims = (-0.05, 0.0107), lc = 4, legend = :right)

    vline!([-495], ls=:dash, lc=1, label="")
    # vline!([-750], ls=:dash, lc=2, label="")
    vline!([-795], ls=:dash, lc=3, label="")
    vline!([-755], ls=:dash, lc=4, label="")


    annotate!(pl, (annotation_pos, text("(h)", annotation_font_size, "Computer Modern")))

    push!(plots, pl)

    # ----------------------------------------------------

    file = h5open(joinpath(@__DIR__, "J_effs_tunneling.h5"), "r")
    J_201 = vec((read(file["J_201"])))
    J_211 = vec((read(file["J_211"])))
    w1 = vec((read(file["w_1"])))

    w1_l = length(w1)
    sl = Int(round(w1_l / 14)) # segment length

    close(file)

    pl = plot(w1 .- 4000, J_201, label = L"J_{201}", xlabel = L"\Delta_{q} / 2\pi~(\textrm{MHz})", ylabel = L"J / 2\pi", ylims = (-0.1, 0.1))
    w1_start = w1[1:Int(round(w1_l / 2 - sl / 2))]
    w1_mid = w1[Int(round(w1_l / 2 - sl / 2)):Int(round(w1_l / 2 + sl / 2))]
    w1_end = w1[Int(round(w1_l / 2 + sl / 2)):end]

    J_211_start = J_211[1:Int(round(w1_l / 2 - sl / 2))]
    J_211_mid = J_211[Int(round(w1_l / 2 - sl / 2)):Int(round(w1_l / 2 + sl / 2))]
    J_211_end = J_211[Int(round(w1_l / 2 + sl / 2)):end]

    mid_n = div(length(w1_mid), 2)

    plot!(pl, w1_start .- 4000, J_211_start, label = "", lc = 2)
    plot!(pl, w1_mid[1:(mid_n)] .- 4000, J_211_mid[1:(mid_n)], label = "", lc = 2, ls = :dash)
    plot!(pl, w1_mid[(mid_n + 2):end] .- 4000, J_211_mid[(mid_n + 2):end], label = "", lc = 2, ls = :dash)
    plot!(pl, w1_end .- 4000, J_211_end, label = L"J_{211}", lc = 2)

    annotate!(pl, (annotation_pos, text("(j)", annotation_font_size, "Computer Modern")))

    plot!(pl, xticks = (xtick_pos, xtick_labels))
    push!(plots, pl)

    # ----------------------------------------------------

    layout = @layout [a; b; c]
    pl = plot(plots..., layout = layout, size = 1.5 .* (400, 3 * 300), leftmargin = 13mm, rightmargin = 13mm, margin = 2mm, guidefont = 20, tickfont = 16, legendfont = 16)
    display(pl)
    savefig(pl, joinpath(@__DIR__, "Fig1_without_scheme.pdf"))
end

make_fig()