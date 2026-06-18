using Plots
using HDF5
using LaTeXStrings
using Plots.Measures

function make_fig()
    filename = joinpath(@__DIR__, "Fig4_data.h5")
    plots = []
    annotation_font_size = 17

    scheme_200 = L"\textrm{Coupler}_{20}"
    scheme_pumped = L"\textrm{Pumped}_{21}"
    scheme_direct = L"\textrm{Direct}"

    annotation_pos = ((-0.18, 1.15))

    h5open(filename, "r") do file
        xtick_pos = [10^-3, 10^-2, 10^-1, 10^0, 10^1]
        xtick_labels = ["0.001", "0.01", "0.1", "1", "10"]
        # Panel (a) ------------------------------------------------------
        plot_a = plot(ylims = (-0.0, 2.02), xaxis = :log10, xticks = (xtick_pos, xtick_labels),xlabel = L"t~\textrm{(\mu s)}", ylabel = L"\left\langle n_1 \right\rangle")
        plot!(plot_a, legend = :bottomleft)
        x_c200 = read(file["a"]["Coupler_200"], "x")
        y_c200 = read(file["a"]["Coupler_200"], "y")
        plot!(plot_a, x_c200, y_c200, label = scheme_200)

        x_pumped = read(file["a"]["Pumped_201"], "x")
        y_pumped = read(file["a"]["Pumped_201"], "y")
        plot!(plot_a, x_pumped, y_pumped, label = scheme_pumped)

        x_direct = read(file["a"]["Direct"], "x")
        y_direct = read(file["a"]["Direct"], "y")
        plot!(plot_a, x_direct, y_direct, label = scheme_direct)

        plot!(plot_a, leftmargin = 10mm)
        annotate!(plot_a, (annotation_pos, text("(a)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_a)

        # Panel (c) ------------------------------------------------------
        plot_c = plot(ylims = (-0.0, 2.02), xaxis = :log10, xticks = (xtick_pos, xtick_labels),xlabel = L"t~\textrm{(\mu s)}", ylabel = L"\left\langle n_c \right\rangle")
        x_c200 = read(file["c"]["Coupler_200"], "x")
        y_c200 = read(file["c"]["Coupler_200"], "y")
        plot!(plot_c, x_c200, y_c200, label = scheme_200)

        x_pumped = read(file["c"]["Pumped_201"], "x")
        y_pumped = read(file["c"]["Pumped_201"], "y")
        plot!(plot_c, x_pumped, y_pumped, label = scheme_pumped)

        plot!(plot_c, leftmargin = 10mm)
        annotate!(plot_c, (annotation_pos, text("(c)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_c)

        # Panel (b) ------------------------------------------------------
        plot_b = plot(ylims = (-0.0, 2.02), xaxis = :log10, xticks = (xtick_pos, xtick_labels),xlabel = L"t~\textrm{(\mu s)}", ylabel = L"\left\langle n_1 \right\rangle")
        x_c200 = read(file["b"]["Coupler_200"], "x")
        y_c200 = read(file["b"]["Coupler_200"], "y")
        plot!(plot_b, x_c200, y_c200, label = scheme_200)

        x_pumped = read(file["b"]["Pumped_201"], "x")
        y_pumped = read(file["b"]["Pumped_201"], "y")
        plot!(plot_b, x_pumped, y_pumped, label = scheme_pumped)

        x_direct = read(file["b"]["Direct"], "x")
        y_direct = read(file["b"]["Direct"], "y")
        plot!(plot_b, x_direct, y_direct, label = scheme_direct)

        plot!(plot_b, leftmargin = 10mm)
        annotate!(plot_b, (annotation_pos, text("(b)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_b)

        # Panel (d) ------------------------------------------------------
        plot_d = plot(ylims = (-0.0, 2.02), xaxis = :log10, xticks = (xtick_pos, xtick_labels),xlabel = L"t~\textrm{(\mu s)}", ylabel = L"\left\langle n_1 \right\rangle")
        x_c200 = read(file["d"]["Coupler_200"], "x")
        y_c200 = read(file["d"]["Coupler_200"], "y")
        plot!(plot_d, x_c200, y_c200, label = scheme_200)

        x_pumped = read(file["d"]["Pumped_201"], "x")
        y_pumped = read(file["d"]["Pumped_201"], "y")
        plot!(plot_d, x_pumped, y_pumped, label = scheme_pumped)

        plot!(plot_d, leftmargin = 10mm)
        annotate!(plot_d, (annotation_pos, text("(d)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_d)
    end
    layout = @layout [a c; b d]
    pl = plot(plots..., layout = layout, size = 3 .* (400, 300), margin = 0mm, bottommargin = 2mm, topmargin = 14mm, guidefont = 18, tickfont = 13, legendfont = 13)
    display(pl)
    savefig(pl, joinpath(@__DIR__, "Fig4.pdf"))
end

function make_b_zoom()
        filename = joinpath(@__DIR__, "Fig4_data.h5")
    plots = []
    annotation_font_size = 17

    scheme_200 = L"\textrm{Coupler}_{20}"
    scheme_pumped = L"\textrm{Pumped}_{21}"
    scheme_direct = L"\textrm{Direct}"

    annotation_pos = ((0.1, 1))

    h5open(filename, "r") do file
        xtick_pos = [10^-2, 10^-1, 10^0, 10^1]
        xtick_labels = ["", "0.1", "1", "10"]

        ytick_pos = [0.9, 0.95, 1]
        ytick_labels = ["0.90", "0.95", "1.00"]

        # Panel (b) ------------------------------------------------------
        plot_b = plot(ylims = (0.9, 1.01), xaxis = :log10, xticks = (xtick_pos, xtick_labels), yticks = (ytick_pos, ytick_labels), xlabel = L"t~\textrm{(\mu s)}", ylabel = L"\left\langle n_1 \right\rangle")
        plot!(plot_b, legend = :none)
        x_c200 = read(file["b"]["Coupler_200"], "x")
        y_c200 = read(file["b"]["Coupler_200"], "y")
        plot!(plot_b, x_c200, y_c200, label = scheme_200, linewidth = 0.5)

        x_pumped = read(file["b"]["Pumped_201"], "x")
        y_pumped = read(file["b"]["Pumped_201"], "y")
        plot!(plot_b, x_pumped, y_pumped, label = scheme_pumped, linewidth = 0.5)

        x_direct = read(file["b"]["Direct"], "x")
        y_direct = read(file["b"]["Direct"], "y")
        plot!(plot_b, x_direct, y_direct, label = scheme_direct, linewidth = 0.5)

        plot!(plot_b, leftmargin = 0mm)
        # annotate!(plot_b, (annotation_pos, text("(b.i)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_b)
    end
    pl = plot(plots..., size = 0.6 .* (400, 300), margin = 0mm, bottommargin = 0mm, topmargin = 0mm, guidefont = 18, tickfont = 13, legendfont = 13)
    display(pl)
    savefig(pl, joinpath(@__DIR__, "b_zoom.pdf"))
end

make_fig()
make_b_zoom()
