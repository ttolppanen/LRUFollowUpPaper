using Graphs
using GraphPlot
using Distributions


function f()
    # --- Create the graph ---
    width = 7
    height = 7
    frequency_threshold = 0.0000 # 50MHz
    n_sites = width * height

    # from Table 8: ibm_brisbane
    frequencies = rand(Normal(4.789996, 0.109104), n_sites)
    anharmonicities = rand(Normal(-0.310393, 0.00548), n_sites)

    g = grid((width, height))

    remove_edges!(g, frequencies, frequency_threshold)

    x_locs = vcat([0:(width - 1) for _ in 1:height]...)
    y_locs = vcat([[i for _ in 1:width] for i in 0:(height - 1)]...)

    gplot(g, x_locs, y_locs, nodelabel=1:n_sites, NODESIZE = 0.5 / sqrt(n_sites), plot_size = (13cm, 13cm))
end


function remove_edges!(g, frequencies, threshold)
    for transmon in vertices(g)
        for neighbor in neighbors(g, transmon)
            if abs(frequencies[transmon] - frequencies[neighbor]) > threshold
                rem_edge!(g, transmon, neighbor)
            end
        end
    end
end


f()