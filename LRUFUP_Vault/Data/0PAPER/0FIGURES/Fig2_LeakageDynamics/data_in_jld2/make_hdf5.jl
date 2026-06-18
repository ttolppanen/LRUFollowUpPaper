using PlotAndSave
using JLD2
using HDF5

function save_hdf5()
    filename = joinpath(@__DIR__, "Fig2_data.h5")
    h5open(filename, "w") do file
        attrs(file)["DOI"] = "10.5281/zenodo.20733347"
        attrs(file)["Description"] = "Data for Fig. 2 of the paper Leakage Mobility and Passive Leakage Removal in Transmons with Tunable Couplers"

        # Panel (a) ------------------------------------------------------
        g = create_group(file, "a")
        attrs(g)["Description"] = "Data for panel (a) of Fig. 2, leakage propagation in the g_eff = 0 tuning. The x-axis shows the time in microseconds. The y-axis shows the population of the first transmon <n1>."

        plotinfo = loaddata(@__DIR__, "200_hopping_data.jld2")
        J = parse(Float64, plotinfo.parameters[:J12][1:(end - 3)])

        h200 = plotinfo.lines["below, s_0 = |200>, <n1>"]
        h100 = plotinfo.lines["below, s_0 = |100>, <n1>"]

        sg = create_group(g, "100")
        attrs(sg)["Description"] = "Data for initial state 100"
        write(sg, "x", collect(h100.x) / (2 * pi * J))
        write(sg, "y", h100.y)

        sg = create_group(g, "200")
        attrs(sg)["Description"] = "Data for initial state 200"
        write(sg, "x", collect(h200.x) / (2 * pi * J))
        write(sg, "y", h200.y)

        # Panel (b) ------------------------------------------------------
        g = create_group(file, "b")
        attrs(g)["Description"] = "Data for panel (b) of Fig. 2, leakage propagation when the second transmon is detuned from the first one. The x-axis shows the frequency of the second transmon in megahertz. The y-axis shows the time in microseconds. The color of the plot indicates the population of the second transmon <n2 / 2> where we have divided by 2 to normalize the data, for which the values here are under the z-label."

        @load joinpath(@__DIR__, "200_off_resonance_data.jld2") x y z p

        write(g, "x", collect(x))
        write(g, "y", collect(y))
        write(g, "z", collect(z))

        # Panel (c) ------------------------------------------------------
        g = create_group(file, "c")
        attrs(g)["Description"] = "Data for panel (c) of Fig. 2, maximum leakage propagation for different detuning of the second transmon. We show the maximum leakage for the g_eff = 0 tuning, and the zz = 0 tuning of the coupler. The x-axis shows the frequency of the second transmon in megahertz. The y-axis shows the maximum population of the second transmon max(<n2 / 2>) during a time of 10 microseconds."

        plotinfo = loaddata(@__DIR__, "200_off_resonance_max_data.jld2")

        g_0 = plotinfo.lines["g = 0, max(<n2 / 2>)"]
        zz_0 = plotinfo.lines["zz = 0 (min), max(<n2 / 2>)"]
        sg = create_group(g, "g_0")
        attrs(sg)["Description"] = "Data for g_eff = 0 tuning"
        write(sg, "x", g_0.x)
        write(sg, "y", Vector{Float64}(g_0.y))

        sg = create_group(g, "zz_0")
        attrs(sg)["Description"] = "Data for zz = 0 tuning"
        write(sg, "x", zz_0.x)
        write(sg, "y", Vector{Float64}(zz_0.y))

        # Panel (d) ------------------------------------------------------
        g = create_group(file, "d")
        attrs(g)["Description"] = "Data for panel (d) of Fig. 2, leakage propagation when the second transmon does not have equal anharmonicity U as the first one. Even though the anharmonicities are different, the second transmon is tuned so that the second excited states of the first and second transmon have equal energy. The x-axis shows the anharmonicity of the second transmon in megahertz. The y-axis shows the time in microseconds. The color of the plot indicates the population of the second transmon <n2 / 2> where we have divided by 2 to normalize the data, for which the values here are under the z-label."

        @load joinpath(@__DIR__, "200_U_data.jld2") x y z p

        write(g, "x", collect(x))
        write(g, "y", collect(y))
        write(g, "z", collect(z))

        # Panel (e) ------------------------------------------------------
        g = create_group(file, "e")
        attrs(g)["Description"] = "Data for panel (e) of Fig. 2, leakage propagation in the g_eff = 0 tuning with an excitation in the second transmon. The x-axis shows the time in nanoseconds. The y-axis shows the population of the first transmon <n1>."

        plotinfo = loaddata(@__DIR__, "201_hopping_data.jld2")
        J = parse(Float64, plotinfo.parameters[:J12][1:(end - 3)])

        h201 = plotinfo.lines["below, s_0 = |201>, <n1>"]
        h100 = plotinfo.lines["below, s_0 = |100>, <n1>"]

        sg = create_group(g, "100")
        attrs(sg)["Description"] = "Data for initial state 100"
        write(sg, "x", collect(h100.x) / (2 * pi * J) * 1e3)
        write(sg, "y", h100.y)

        sg = create_group(g, "201")
        attrs(sg)["Description"] = "Data for initial state 201"
        write(sg, "x", collect(h201.x) / (2 * pi * J) * 1e3)
        write(sg, "y", h201.y)

        # Panel (f) ------------------------------------------------------
        g = create_group(file, "f")
        attrs(g)["Description"] = "Data for panel (f) of Fig. 2, leakage propagation when the second transmon is detuned from the first one while the second transmon is excited. The x-axis shows the frequency of the second transmon in megahertz. The y-axis shows the time in nanoseconds. The color of the plot indicates the population of the second transmon <n2 - I> where we have used n2 - I to normalize the data between 0 and 1, for which the values here are under the z-label."

        @load joinpath(@__DIR__, "201_off_resonance_data.jld2") x y z p

        write(g, "x", collect(x))
        write(g, "y", collect(y) .* 1e3)
        write(g, "z", collect(z))

        # Panel (g) ------------------------------------------------------
        g = create_group(file, "g")
        attrs(g)["Description"] = "Data for panel (g) of Fig. 2, maximum leakage propagation for different detuning of the second transmon. We show the maximum leakage for the g_eff = 0 tuning, and the zz = 0 tuning of the coupler. The x-axis shows the frequency of the second transmon in megahertz. The y-axis shows the maximum population of the second transmon max(<n2 - I>) during a time of 3 microseconds."

        plotinfo = loaddata(@__DIR__, "201_off_resonance_max_data.jld2")

        g_0 = plotinfo.lines["g = 0, max(<n2 - I>)"]
        zz_0 = plotinfo.lines["zz = 0 (min), max(<n2 - I>)"]
        sg = create_group(g, "g_0")
        attrs(sg)["Description"] = "Data for g_eff = 0 tuning"
        write(sg, "x", g_0.x)
        write(sg, "y", Vector{Float64}(g_0.y))

        sg = create_group(g, "zz_0")
        attrs(sg)["Description"] = "Data for zz = 0 tuning"
        write(sg, "x", zz_0.x)
        write(sg, "y", Vector{Float64}(zz_0.y))

        # Panel (h) ------------------------------------------------------
        g = create_group(file, "h")
        attrs(g)["Description"] = "Data for panel (h) of Fig. 2, leakage propagation when the second transmon does not have equal anharmonicity U as the first one. Even though the anharmonicities are different, the second transmon is tuned so that the second excited states of the first and second transmon have equal energy. The x-axis shows the anharmonicity of the second transmon in megahertz. The y-axis shows the time in nanoseconds. The color of the plot indicates the population of the second transmon <n2 - I> where we have used n2 - I to normalize the data between 0 and 1, for which the values here are under the z-label."

        @load joinpath(@__DIR__, "201_U_data.jld2") x y z p

        write(g, "x", collect(x))
        write(g, "y", collect(y) .* 1e3)
        write(g, "z", collect(z))
    end
end

save_hdf5()