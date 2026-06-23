using Plots
using HDF5
using LaTeXStrings
using Plots.Measures

function make_fig()
    filename = joinpath(@__DIR__, "Fig5_data.h5")
    plots = []
    annotation_font_size = 17
    annotation_pos = ((-0.17, 1.15))

    t1_legend = L"\textrm{W/o~LRU}"
    lru_legend = L"\textrm{With~LRU}"

    h5open(filename, "r") do file
        # Panel (a) ------------------------------------------------------
        plot_a = plot(ylims = (0, 2.02), xlabel = L"t~\textrm{(ns)}", ylabel = L"\left\langle \hat{n}_1 \right\rangle")
        x_leak = read(file["a"], "x")
        y_leak = read(file["a"], "y")
        plot!(plot_a, x_leak, y_leak, legend = :none)

        plot!(plot_a, leftmargin = 10mm)
        annotate!(plot_a, (annotation_pos, text("(a)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_a)

        # Panel (b) ------------------------------------------------------
        plot_b = plot(framestyle = :none, axis = nothing, legend = false)

        annotate!(plot_b, ((-0.0, 1.15), text("(b)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_b)

        # Panel (c) ------------------------------------------------------
        plot_b = plot(ylims = (0, 1.02), xlabel = L"t~" * "(μs)", ylabel = L"\left\langle \hat{n}_1 \right\rangle")
        
        x_leak = read(file["c"]["T1 and LRU"], "x")
        y_leak = read(file["c"]["T1 and LRU"], "y")
        plot!(plot_b, x_leak, y_leak, label = lru_legend)
        
        x_leak = read(file["c"]["T1 only"], "x")
        y_leak = read(file["c"]["T1 only"], "y")
        plot!(plot_b, x_leak, y_leak, label = t1_legend)

        annotate!(plot_b, (annotation_pos, text("(c)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_b)

        # Panel (d) ------------------------------------------------------
        plot_d = plot(ylims = (-0.001, 0.082), xlabel = L"t~" * "(μs)", ylabel = L"\left\langle \hat{n}_1 \right\rangle")
        
        x_lru_00 = read(file["d"]["T1 and LRU"], "x")
        y_lru_00 = read(file["d"]["T1 and LRU"], "y")
        plot!(plot_d, x_lru_00, y_lru_00, label = lru_legend)
        
        x_t1_00 = read(file["d"]["T1 only"], "x")
        y_t1_00 = read(file["d"]["T1 only"], "y")
        plot!(plot_d, x_t1_00, y_t1_00, label = t1_legend)

        annotate!(plot_d, (annotation_pos, text("(d)", annotation_font_size, "Computer Modern")))
        push!(plots, plot_d)
    end
    layout = @layout [a c; b d]
    pl = plot(plots..., layout = layout, size = 3 .* (400, 300), margin = 0mm, bottommargin = 2mm, topmargin = 14mm, guidefont = 18, tickfont = 13, legendfont = 13)
    display(pl)
    savefig(pl, joinpath(@__DIR__, "Fig5_without_circuit.pdf"))
end

make_fig()