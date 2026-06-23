using Plots
using HDF5
using LaTeXStrings
using Plots.Measures

function make_fig()
    filename = joinpath(@__DIR__, "Fig3_data.h5")
    plots = []
    annotation_font_size = 28
    # colorbar_lims_201 = (1, 2)
    # colorbar_lims_211 = (1, 2)
    annotation_pos = (-0.12, 1.15)
    h5open(filename, "r") do file
        # Panel (a) ------------------------------------------------------
        plot_a = plot(ylims = (-0.02, 2.02), xlabel = L"t~" * "(μs)")
        x_n1 = read(file["a"]["n1"], "x")
        y_n1 = read(file["a"]["n1"], "y")
        plot!(plot_a, x_n1, y_n1, label = L"\left\langle \hat{n}_1 \right\rangle")

        x_n2 = read(file["a"]["n2"], "x")
        y_n2 = read(file["a"]["n2"], "y")
        plot!(plot_a, x_n2, y_n2, label = L"\left\langle \hat{n}_2 \right\rangle")

        x_n3 = read(file["a"]["n3"], "x")
        y_n3 = read(file["a"]["n3"], "y")
        plot!(plot_a, x_n3, y_n3, label = L"\left\langle \hat{n}_3 \right\rangle")

        # x_nc1 = read(file["a"]["nc1"], "x")
        # y_nc1 = read(file["a"]["nc1"], "y")
        # plot!(plot_a, x_nc1, y_nc1, label = L"|n_{c_1}>")

        # x_nc2 = read(file["a"]["nc2"], "x")
        # y_nc2 = read(file["a"]["nc2"], "y")
        # plot!(plot_a, x_nc2, y_nc2, label = L"|n_{c_2}>")

        plot!(plot_a, leftmargin = 15mm, legend = :bottomright)
        annotate!(plot_a, (annotation_pos, text("(a)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_a)

        # Panel (b) ------------------------------------------------------
        x = read(file["b"], "x")
        y = read(file["b"], "y")
        z = read(file["b"], "z")

        plot_b = heatmap(x, y, z, xlabel = L"\omega_1 / 2\pi = \omega_3 / 2\pi~\textrm{(MHz)}", ylabel = L"t~" * "(μs)", colorbar_title = L"\left\langle \hat{n}_3 \right\rangle", clims = (0, 2))
        annotate!(plot_b, (annotation_pos, text("(b)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_b)

        # Panel (c) ------------------------------------------------------
        x = read(file["c"], "x")
        y = read(file["c"], "y")
        z = read(file["c"], "z")

        plot_c = heatmap(x, y, z, xlabel = L"\omega_1 / 2\pi~\textrm{(MHz)}", ylabel = L"\left( \omega_3 - \omega_1 \right) / 2\pi ~\textrm{(MHz)}", colorbar_title = L"\textrm{max}\left( \left\langle \hat{n}_3 \right\rangle \right)", clims = (1, 2))
        plot!(plot_c, rightmargin = 10mm)
        annotate!(plot_c, (annotation_pos, text("(c)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_c)

        # Panel (d) ------------------------------------------------------
        plot_d = plot(ylims = (-0.02, 2.02), xlabel = L"t~" * "(μs)")
        x_n1 = read(file["d"]["n1"], "x")
        y_n1 = read(file["d"]["n1"], "y")
        plot!(plot_d, x_n1, y_n1, label = L"\left\langle \hat{n}_1 \right\rangle")

        x_n2 = read(file["d"]["n2"], "x")
        y_n2 = read(file["d"]["n2"], "y")
        plot!(plot_d, x_n2, y_n2, label = L"\left\langle \hat{n}_2 \right\rangle")

        x_n3 = read(file["d"]["n3"], "x")
        y_n3 = read(file["d"]["n3"], "y")
        plot!(plot_d, x_n3, y_n3, label = L"\left\langle \hat{n}_3 \right\rangle")

        # x_nc1 = read(file["d"]["nc1"], "x")
        # y_nc1 = read(file["d"]["nc1"], "y")
        # plot!(plot_d, x_nc1, y_nc1, label = L"|n_{c_1}>")

        # x_nc2 = read(file["d"]["nc2"], "x")
        # y_nc2 = read(file["d"]["nc2"], "y")
        # plot!(plot_d, x_nc2, y_nc2, label = L"|n_{c_2}>")

        plot!(plot_d, leftmargin = 10mm)
        annotate!(plot_d, (annotation_pos, text("(d)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_d)

        # Panel (e) ------------------------------------------------------
        x = read(file["e"], "x")
        y = read(file["e"], "y")
        z = read(file["e"], "z")

        plot_e = heatmap(x, y, z, xlabel = L"\omega_1 / 2\pi = \omega_3 / 2\pi~\textrm{(MHz)}", ylabel = L"t~" * "(μs)", colorbar_title = L"\left\langle \hat{n}_3 \right\rangle", clims = (0, 2))
        annotate!(plot_e, (annotation_pos, text("(e)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_e)

        # Panel (f) ------------------------------------------------------
        x = read(file["f"], "x")
        y = read(file["f"], "y")
        z = read(file["f"], "z")

        plot_f = heatmap(x, y, z, xlabel = L"\omega_1 / 2\pi~\textrm{(MHz)}", ylabel = L"\left( \omega_3 - \omega_1 \right) / 2\pi ~\textrm{(MHz)}", colorbar_title = L"\textrm{max}\left( \left\langle \hat{n}_3 \right\rangle \right)", clims = (1, 2))
        plot!(plot_f, rightmargin = 10mm)
        annotate!(plot_f, (annotation_pos, text("(f)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_f)
    end
    layout = @layout [a b{0.37w} c{0.37w}; d e f]
    pl = plot(plots..., layout = layout, size = 2 .* (3 * 400, 2 * 300), margin = 12mm, bottommargin = 18mm, topmargin = 18mm, guidefont = 28, tickfont = 22, legendfont = 24, colorbar_titlefontsize = 28)
    display(pl)
    savefig(pl, joinpath(@__DIR__, "Fig3.pdf"))
end

make_fig()