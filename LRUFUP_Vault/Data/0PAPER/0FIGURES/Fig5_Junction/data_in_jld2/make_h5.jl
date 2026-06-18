using PlotAndSave
using JLD2
using HDF5

function save_hdf5()
    filename = joinpath(@__DIR__, "Fig5_data.h5")
    t_to_nanoseconds(times) = times ./ (2 * π * 1) # Assuming J = 1GHz!!!
    h5open(filename, "w") do file
        attrs(file)["DOI"] = "10.5281/zenodo.20733347"
        attrs(file)["Description"] = "Data for Fig. 5 of the paper Leakage Mobility and Passive Leakage Removal in Transmons with Tunable Couplers"

        # Panel (a) ------------------------------------------------------
        g = create_group(file, "a")
        attrs(g)["Description"] = "Data for panel (a) of Fig. 5, leakage removal for the junction readout leakage removal scheme. The x-axis shows the time in nanoseconds. The y-axis shows the population of the transmon <n1>."

        plotinfo = loaddata(@__DIR__, "leakage_data.jld2")

        leakage = plotinfo.lines["LRU dissipation, n1"]

        write(g, "x", t_to_nanoseconds(collect(leakage.x)))
        write(g, "y", leakage.y)

        # Panel (b) ------------------------------------------------------
        g = create_group(file, "b")
        attrs(g)["Description"] = "Data for panel (b) of Fig. 5, qubit decay with and withouth leakage removal on, for a qubit with T1 = 50mus. The x-axis shows the time in microseconds. The y-axis shows the population of the transmon <n1>."

        plotinfo = loaddata(@__DIR__, "qubit_data.jld2")

        qubit_t1 = plotinfo.lines["T1 only, n1"]
        qubit_t1_lru = plotinfo.lines["T1 and LRU dissipation, n1"]

        sg = create_group(g, "T1 only")
        attrs(sg)["Description"] = "Data for the qubit decay without leakage removal."
        write(sg, "x", t_to_nanoseconds(collect(qubit_t1.x)) * 1e-3) # Convert to microseconds
        write(sg, "y", qubit_t1.y)

        sg = create_group(g, "T1 and LRU")
        attrs(sg)["Description"] = "Data for the qubit decay with leakage removal."
        write(sg, "x", t_to_nanoseconds(collect(qubit_t1_lru.x)) * 1e-3) # Convert to microseconds
        write(sg, "y", qubit_t1_lru.y)

        # Panel (d) ------------------------------------------------------
        g = create_group(file, "d")
        attrs(g)["Description"] = "Data for panel (d) of Fig. 5, groundstate dynamics following from the two-mode squeezing. The x-axis shows the time in microseconds. The y-axis shows the population of the transmon <n1>."

        plotinfo = loaddata(@__DIR__, "00_data.jld2")

        ground_t1 = plotinfo.lines["T1 only, n1"]
        ground_t1_lru = plotinfo.lines["T1 and LRU dissipation, n1"]

        sg = create_group(g, "T1 only")
        attrs(sg)["Description"] = "Data for the qubit decay without leakage removal."
        write(sg, "x", t_to_nanoseconds(collect(ground_t1.x)) * 1e-3) # Convert to microseconds
        write(sg, "y", ground_t1.y)

        sg = create_group(g, "T1 and LRU")
        attrs(sg)["Description"] = "Data for the qubit decay with leakage removal."
        write(sg, "x", t_to_nanoseconds(collect(ground_t1_lru.x)) * 1e-3) # Convert to microseconds
        write(sg, "y", ground_t1_lru.y)
    end
end

save_hdf5()