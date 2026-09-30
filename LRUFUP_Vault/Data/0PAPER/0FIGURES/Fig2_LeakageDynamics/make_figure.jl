using Plots
using HDF5
using LaTeXStrings
using Plots.Measures

function make_fig()
    filename = joinpath(@__DIR__, "Fig2_data.h5")
    plots = []
    annotation_font_size = 17
    annotation_pos = (-0.36, 1.15)
    annotation_pos_h = (-0.3, 1.15)
    h5open(filename, "r") do file
        # Panel (a) ------------------------------------------------------
        plot_a = plot(ylims = (0, 2.02), xlabel = L"t~" * "(μs)", ylabel = L"\left\langle \hat{n}_1 \right\rangle")
        x_h100 = read(file["a"]["100"], "x")
        y_h100 = read(file["a"]["100"], "y")
        plot!(plot_a, x_h100, y_h100, label = "|10>")

        x_h200 = read(file["a"]["200"], "x")
        y_h200 = read(file["a"]["200"], "y")
        plot!(plot_a, x_h200, y_h200, label = "|20>")

        plot!(plot_a, leftmargin = 10mm, topmargin = 10mm)
        annotate!(plot_a, (annotation_pos, text("(a)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_a)

        # Panel (b) ------------------------------------------------------
        x = read(file["b"], "x")
        y = read(file["b"], "y")
        z = read(file["b"], "z")

        plot_b = heatmap(x, y, z, xlabel = L"\Delta_{12} / 2\pi~\textrm{(MHz)}", ylabel = L"t~" * "(μs)", colorbar_title = L"\left\langle \hat{n}_2 / 2 \right\rangle")
        plot!(plot_b, clims = (0, 2), colorbar = false)
        annotate!(plot_b, (annotation_pos_h, text("(b)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_b)

        # Panel (c) ------------------------------------------------------
        plot_c = plot(xlabel = L"\Delta_{12} / 2\pi~\textrm{(MHz)}", ylabel = L"\max\left(\left\langle \hat{n}_2 (t) \right\rangle\right)", ylims = (-0.02, 2))

        x_g_0 = read(file["c"]["g_0"], "x")
        y_g_0 = read(file["c"]["g_0"], "y")
        plot!(plot_c, x_g_0, y_g_0, label = L"g_{\textrm{eff}} = 0")

        x_zz_0 = read(file["c"]["zz_0"], "x")
        y_zz_0 = read(file["c"]["zz_0"], "y")
        plot!(plot_c, x_zz_0, y_zz_0, label = L"\zeta = 0")

        annotate!(plot_c, (annotation_pos, text("(c)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_c)

        # Panel (d) ------------------------------------------------------
        x = read(file["d"], "x")
        y = read(file["d"], "y")
        z = read(file["d"], "z")

        plot_d = heatmap(x, y, z, xlabel = L"\Delta U_{12} / 2\pi~\textrm{(MHz)}", ylabel = L"t~"* "(μs)", colorbar_title = L"\left\langle \hat{n}_2 - 1 \right\rangle", clims = (0, 2))
        plot!(plot_d, colorbar = false, xlims = (-50, 50))
        annotate!(plot_d, (annotation_pos_h, text("(d)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_d)

        # Colorbar
        limits = collect(range(0, 2, length = 100))
        colorbar = heatmap([1], limits, [limits;;], colorbar=false, xaxis=false, ymirror=true, ylabel = L"\left\langle \hat{n}_2 \right\rangle")
        plot!(colorbar, rightmargin = 12mm)
        push!(plots, colorbar)

        # Panel (e) ------------------------------------------------------
        plot_e = plot(ylims = (0, 2.02), xlabel = L"t~\textrm{(ns)}", ylabel = L"\left\langle \hat{n}_1 \right\rangle")
        x_h100 = read(file["e"]["100"], "x")
        y_h100 = read(file["e"]["100"], "y")
        plot!(plot_e, x_h100, y_h100, label = "|10>")

        x_h201 = read(file["e"]["201"], "x")
        y_h201 = read(file["e"]["201"], "y")
        plot!(plot_e, x_h201, y_h201, label = "|21>")
        plot!(plot_e, bottommargin = 8mm)

        annotate!(plot_e, (annotation_pos, text("(e)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_e)

        # Panel (f) ------------------------------------------------------
        x = read(file["f"], "x")
        y = read(file["f"], "y")
        z = read(file["f"], "z")

        plot_f = heatmap(x, y, z, xlabel = L"\Delta_{12} / 2\pi~\textrm{(MHz)}", ylabel = L"t~\textrm{(ns)}", colorbar_title = L"\left\langle \hat{n}_2 - 1 \right\rangle")
        plot!(plot_f, clims = (1, 2), colorbar = false)
        annotate!(plot_f, (annotation_pos_h, text("(f)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_f)

        # Panel (g) ------------------------------------------------------
        plot_g = plot(xlabel = L"\Delta_{12} / 2\pi~\textrm{(MHz)}", ylabel = L"\max\left(\left\langle \hat{n}_2 (t) \right\rangle\right)", ylims = (0.99, 2))

        x_g_0 = read(file["g"]["g_0"], "x")
        y_g_0 = read(file["g"]["g_0"], "y")
        plot!(plot_g, x_g_0, y_g_0, label = L"g_{\textrm{eff}} = 0")

        x_zz_0 = read(file["g"]["zz_0"], "x")
        y_zz_0 = read(file["g"]["zz_0"], "y")
        plot!(plot_g, x_zz_0, y_zz_0, label = L"\zeta = 0")

        annotate!(plot_g, (annotation_pos, text("(g)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_g)

        # Panel (h) ------------------------------------------------------
        x = read(file["h"], "x")
        y = read(file["h"], "y")
        z = read(file["h"], "z")

        plot_h = heatmap(x, y, z, xlabel = L"\Delta U_{12} / 2\pi~\textrm{(MHz)}", ylabel = L"t~\textrm{(ns)}", colorbar_title = L"\left\langle \hat{n}_2 - 1 \right\rangle", clims = (1, 2))
        plot!(plot_h, colorbar = false, xlims = (-50, 50))
        annotate!(plot_h, (annotation_pos_h, text("(h)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_h)

        # Colorbar
        limits = collect(range(0, 1, length = 100))
        colorbar = heatmap([1], limits, [limits;;], colorbar=false, xaxis=false, ymirror=true, ylabel = L"\left\langle \hat{n}_2 \right\rangle")
        plot!(colorbar, rightmargin = 12mm)
        push!(plots, colorbar)
    end
    layout = @layout [a b c d e{0.02w}; f g h i j{0.02w}]
    pl = plot(plots..., layout = layout, size = 1 .* (4 * 400, 2 * 300), margin = 4mm, guidefont = 18, tickfont = 13, legendfont = 13)
    display(pl)
    savefig(pl, joinpath(@__DIR__, "Fig2.pdf"))
end

function compare_coupling()
    filename = joinpath(@__DIR__, "Fig2_data.h5")
    plots = []
    h5open(filename, "r") do file
        # Panel (c) ------------------------------------------------------
        plot_c = plot(xlabel = L"\Delta_{12} / 2\pi~\textrm{(MHz)}", ylabel = L"\max\left(\left\langle \hat{n}_2 (t) \right\rangle\right)")
        plot!(plot_c, ylims = (0, 2.1))

        x_g_0 = read(file["c"]["g_0"], "x")
        y_g_0 = read(file["c"]["g_0"], "y")
        plot!(plot_c, x_g_0 .- 4000, 1.0138 * y_g_0, label = L"g_{\textrm{eff}} = 0")

        x_zz_0 = read(file["c"]["zz_0"], "x")
        y_zz_0 = read(file["c"]["zz_0"], "y")
        # plot!(plot_c, x_zz_0 .- 4000, y_zz_0, label = L"\zeta = 0")

        plot!(plot_c, x -> 2 / sqrt(1 + (x / 0.2)^2), label = "ref")

        display(plot_c)
    end
end


make_fig()
# compare_coupling()