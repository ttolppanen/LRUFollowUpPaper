using PlotAndSave
using HDF5

function make_h5()
    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_10_10mus"))
    d_10 = plotinfo.lines["|10>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_100_10mus"))
    c200_10 = plotinfo.lines["|100>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_100_10mus"))
    pumped_10 = plotinfo.lines["|100>, n1"]

    plotinfo = loaddata(joinpath(@__DIR__, "DirectCoupling", "result_20_10mus"))
    d_20 = plotinfo.lines["|20>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_20_10mus"))
    c200_20 = plotinfo.lines["|200>, n1"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_20_10mus"))
    pumped_20 = plotinfo.lines["|200>, n_1"]

    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_010_10mus"))
    c200_010 = plotinfo.lines["|010>, nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_010_10mus"))
    pumped_010 = plotinfo.lines["|010>, nc"]

    plotinfo = loaddata(joinpath(@__DIR__, "Coupler200", "result_110_10mus"))
    c200_110 = plotinfo.lines["|110>, n1"]
    # c200_10_nc = plotinfo.lines["|110>, nc"]
    plotinfo = loaddata(joinpath(@__DIR__, "Pumped201", "result_110_10mus"))
    pumped_110 = plotinfo.lines["|110>, n1"]
    # pumped_10_nc = plotinfo.lines["|110>, nc"]

    t_to_mus(times) = times ./ (2 * π * 5) # ASSUMING IN UNITS OF J = 5MHZ !!!!
    plots = []

    scheme_200 = "Coupler_200"
    scheme_pumped = "Pumped_201"
    scheme_direct = "Direct"
    
    filepath = joinpath(@__DIR__, "Fig4_data.h5")
    h5open(filepath, "w") do file
        attrs(file)["DOI"] = "10.5281/zenodo.20733347"
        attrs(file)["Description"] = "Data for Fig. 4 of the paper Leakage Mobility and Passive Leakage Removal in Transmons with Tunable Couplers"

        # Panel (a) ------------------------------------------------------
        g = create_group(file, "a")
        attrs(g)["Description"] = "Data for panel (a) of Fig. 4, leakage removal for the three schemes $(scheme_200), $(scheme_pumped) and $(scheme_direct). The x-axis shows the time in microseconds. The time-steps are separeted so that on a log scale they have even distribution. The y-axis shows the population of the first transmon <n1>, starting from a leaked initial state |200> (|20> for $(scheme_direct))."

        sg = create_group(g, scheme_200)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_200)"
        write(sg, "x", t_to_mus(collect(c200_20.x)))
        write(sg, "y", c200_20.y)

        sg = create_group(g, scheme_pumped)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_pumped)"
        write(sg, "x", t_to_mus(collect(pumped_20.x)))
        write(sg, "y", pumped_20.y)

        sg = create_group(g, scheme_direct)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_direct)"
        write(sg, "x", t_to_mus(collect(d_20.x)))
        write(sg, "y", d_20.y)

        # Panel (b) ------------------------------------------------------
        g = create_group(file, "b")
        attrs(g)["Description"] = "Data for panel (b) of Fig. 4, qubit decay for the three schemes $(scheme_200), $(scheme_pumped) and $(scheme_direct). The x-axis shows the time in microseconds. The time-steps are separeted so that on a log scale they have even distribution. The y-axis shows the population of the first transmon <n1>, starting from a leaked initial state |100> (|10> for $(scheme_direct)."

        sg = create_group(g, scheme_200)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_200)"
        write(sg, "x", t_to_mus(collect(c200_10.x)))
        write(sg, "y", c200_10.y)

        sg = create_group(g, scheme_pumped)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_pumped)"
        write(sg, "x", t_to_mus(collect(pumped_10.x)))
        write(sg, "y", pumped_10.y)

        sg = create_group(g, scheme_direct)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_direct)"
        write(sg, "x", t_to_mus(collect(d_10.x)))
        write(sg, "y", d_10.y)

        # Panel (c) ------------------------------------------------------
        g = create_group(file, "c")
        attrs(g)["Description"] = "Data for panel (c) of Fig. 4, leakage removal from the coupler for the two coupled enabled schemes $(scheme_200) and $(scheme_pumped). The x-axis shows the time in microseconds. The time-steps are separeted so that on a log scale they have even distribution. The y-axis shows the population of the coupler <nc>, starting from a leaked initial state |010>."

        sg = create_group(g, scheme_200)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_200)"
        write(sg, "x", t_to_mus(collect(c200_010.x)))
        write(sg, "y", c200_010.y)

        sg = create_group(g, scheme_pumped)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_pumped)"
        write(sg, "x", t_to_mus(collect(pumped_010.x)))
        write(sg, "y", pumped_010.y)

        # Panel (d) ------------------------------------------------------
        g = create_group(file, "d")
        attrs(g)["Description"] = "Data for panel (d) of Fig. 4, qubit decay if there is leakage in the coupler for the two coupler enabled schemes $(scheme_200) and $(scheme_pumped). The x-axis shows the time in microseconds. The time-steps are separeted so that on a log scale they have even distribution. The y-axis shows the population of the first transmon <n1>, starting from a leaked initial state |110>."

        sg = create_group(g, scheme_200)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_200)"
        write(sg, "x", t_to_mus(collect(c200_110.x)))
        write(sg, "y", c200_110.y)

        sg = create_group(g, scheme_pumped)
        attrs(sg)["Description"] = "Data for the leakage removal scheme $(scheme_pumped)"
        write(sg, "x", t_to_mus(collect(pumped_110.x)))
        write(sg, "y", pumped_110.y)
    end
end

make_h5()