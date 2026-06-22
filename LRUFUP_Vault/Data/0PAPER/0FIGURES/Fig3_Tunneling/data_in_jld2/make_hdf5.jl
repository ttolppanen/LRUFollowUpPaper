using PlotAndSave
using JLD2
using HDF5

function save_hdf5()
    filename = joinpath(@__DIR__, "Fig3_data.h5")
    h5open(filename, "w") do file
        attrs(file)["DOI"] = "10.5281/zenodo.20733347"
        attrs(file)["Description"] = "Data for Fig. 3 of the paper Leakage Mobility and Passive Leakage Removal in Transmons with Tunable Couplers"

        # Panel (a) ------------------------------------------------------
        g = create_group(file, "a")
        attrs(g)["Description"] = "Data for panel (a) of Fig. 3, leakage tunneling time dynamics for initial state 201. The x-axis shows the time in microseconds. The y-axis shows the populations of the transmon and couplers."

        plotinfo = loaddata(@__DIR__, "201_time_data.jld2")
        J = parse(Float64, plotinfo.parameters[:J12][1:(end - 3)])

        n1_201 = plotinfo.lines["s_0 = |201>, <n1>"]
        n2_201 = plotinfo.lines["s_0 = |201>, <n2>"]
        n3_201 = plotinfo.lines["s_0 = |201>, <n3>"]
        nc1_201 = plotinfo.lines["s_0 = |201>, <nc1>"]
        nc2_201 = plotinfo.lines["s_0 = |201>, <nc2>"]

        sg = create_group(g, "n1")
        attrs(sg)["Description"] = "Data for initial state 201, population of the first transmon qubit <n1>"
        write(sg, "x", collect(n1_201.x) / (2 * pi * J))
        write(sg, "y", n1_201.y)

        sg = create_group(g, "n2")
        attrs(sg)["Description"] = "Data for initial state 201, population of the second transmon qubit <n2>"
        write(sg, "x", collect(n2_201.x) / (2 * pi * J))
        write(sg, "y", n2_201.y)

        sg = create_group(g, "n3")
        attrs(sg)["Description"] = "Data for initial state 201, population of the third transmon qubit <n3>"
        write(sg, "x", collect(n3_201.x) / (2 * pi * J))
        write(sg, "y", n3_201.y)

        sg = create_group(g, "nc1")
        attrs(sg)["Description"] = "Data for initial state 201, population of the first tunable coupler transmon <nc1>"
        write(sg, "x", collect(nc1_201.x) / (2 * pi * J))
        write(sg, "y", nc1_201.y)

        sg = create_group(g, "nc2")
        attrs(sg)["Description"] = "Data for initial state 201, population of the second tunable coupler transmon <nc2>"
        write(sg, "x", collect(nc2_201.x) / (2 * pi * J))
        write(sg, "y", nc2_201.y)

        # Panel (b) ------------------------------------------------------
        g = create_group(file, "b")
        attrs(g)["Description"] = "Data for panel (b) of Fig. 3, leakage tunneling with different tunings related to the middle transmon for initial state 201. The x-axis shows the frequency of the first and third transmon in megahertz. The y-axis shows the time in microseconds. The color of the plot indicates the population of the third transmon <n3> for which the values here are under the z-label."

        @load joinpath(@__DIR__, "201_detuned_data.jld2") x y z p

        write(g, "x", collect(x))
        write(g, "y", collect(y))
        write(g, "z", collect(z))

        # Panel (c) ------------------------------------------------------
        g = create_group(file, "c")
        attrs(g)["Description"] = "Data for panel (c) of Fig. 3, leakage tunneling with off-resonant initial (first) and target (third) sites for initial state 201. The x-axis shows the frequency of the first transmon in megahertz. The y-axis shows the detuning of the first and third transmons w3 - w1 in megahertz. The color of the plot indicates the maximal population of the third transmon <n3> over a period of 10 microseconds, for which the values here are under the z-label."

        @load joinpath(@__DIR__, "201_relative_freq_data.jld2") x y max_amps

        write(g, "x", collect(x))
        write(g, "y", collect(y))
        write(g, "z", collect(max_amps))

        # Panel (d) ------------------------------------------------------
        g = create_group(file, "d")
        attrs(g)["Description"] = "Data for panel (d) of Fig. 3, leakage tunneling time dynamics for initial state 211. The x-axis shows the time in microseconds. The y-axis shows the populations of the transmon and couplers."

        plotinfo = loaddata(@__DIR__, "211_time_data.jld2")
        J = parse(Float64, plotinfo.parameters[:J12][1:(end - 3)])

        n1_211 = plotinfo.lines["s_0 = |211>, <n1>"]
        n2_211 = plotinfo.lines["s_0 = |211>, <n2>"]
        n3_211 = plotinfo.lines["s_0 = |211>, <n3>"]
        nc1_211 = plotinfo.lines["s_0 = |211>, <nc1>"]
        nc2_211 = plotinfo.lines["s_0 = |211>, <nc2>"]

        sg = create_group(g, "n1")
        attrs(sg)["Description"] = "Data for initial state 211, population of the first transmon qubit <n1>"
        write(sg, "x", collect(n1_211.x) / (2 * pi * J))
        write(sg, "y", n1_211.y)

        sg = create_group(g, "n2")
        attrs(sg)["Description"] = "Data for initial state 211, population of the second transmon qubit <n2>"
        write(sg, "x", collect(n2_211.x) / (2 * pi * J))
        write(sg, "y", n2_211.y)

        sg = create_group(g, "n3")
        attrs(sg)["Description"] = "Data for initial state 211, population of the third transmon qubit <n3>"
        write(sg, "x", collect(n3_211.x) / (2 * pi * J))
        write(sg, "y", n3_211.y)

        sg = create_group(g, "nc1")
        attrs(sg)["Description"] = "Data for initial state 211, population of the first tunable coupler transmon <nc1>"
        write(sg, "x", collect(nc1_211.x) / (2 * pi * J))
        write(sg, "y", nc1_211.y)

        sg = create_group(g, "nc2")
        attrs(sg)["Description"] = "Data for initial state 211, population of the second tunable coupler transmon <nc2>"
        write(sg, "x", collect(nc2_211.x) / (2 * pi * J))
        write(sg, "y", nc2_211.y)

        # Panel (e) ------------------------------------------------------
        g = create_group(file, "e")
        attrs(g)["Description"] = "Data for panel (e) of Fig. 3, leakage tunneling with different tunings related to the middle transmon for initial state 211. The x-axis shows the frequency of the first and third transmon in megahertz. The y-axis shows the time in microseconds. The color of the plot indicates the population of the third transmon <n3> for which the values here are under the z-label."

        @load joinpath(@__DIR__, "211_detuned_data.jld2") x y z p

        write(g, "x", collect(x))
        write(g, "y", collect(y))
        write(g, "z", collect(z))

        # Panel (f) ------------------------------------------------------
        g = create_group(file, "f")
        attrs(g)["Description"] = "Data for panel (f) of Fig. 3, leakage tunneling with off-resonant initial (first) and target (third) sites for initial state 211. The x-axis shows the frequency of the first transmon in megahertz. The y-axis shows the detuning of the first and third transmons w3 - w1 in megahertz. The color of the plot indicates the maximal population of the third transmon <n3> over a period of 10 microseconds, for which the values here are under the z-label."

        @load joinpath(@__DIR__, "211_relative_freq_data.jld2") x y z p

        write(g, "x", collect(x))
        write(g, "y", collect(y))
        write(g, "z", collect(z))
    end
    println("Data saved!")
end

save_hdf5()