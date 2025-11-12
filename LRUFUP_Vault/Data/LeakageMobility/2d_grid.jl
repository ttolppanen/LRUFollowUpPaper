using Graphs
using GraphPlot
using Distributions
using Statistics


function f()
    # --- Create the graph ---
    width = 7
    height = 7
    frequency_threshold = 40 * 1e-3 # 50MHz
    tunneling_threshold = 10 * 1e-3 # 50MHz
    n_sites = width * height

    # from Table 8: ibm_brisbane
    frequencies = rand(Normal(4.789996, 0.109104), n_sites)
    # anharmonicities = rand(Normal(-0.310393, 0.00548), n_sites)

    g = grid((width, height))

    remove_edges!(g, frequencies, frequency_threshold)
    add_tunneling_edges!(g, frequencies, tunneling_threshold, width, height)

    x_locs = vcat([0:(width - 1) for _ in 1:height]...)
    y_locs = vcat([[i for _ in 1:width] for i in 0:(height - 1)]...)

    pl = gplot(g, x_locs, y_locs, nodelabel=1:n_sites, NODESIZE = 0.5 / sqrt(n_sites), plot_size = (13cm, 13cm))

    ratios = []
    for _ in 1:1000
        frequencies = rand(Normal(4.789996, 0.109104), n_sites)
        g = grid((width, height))
        remove_edges!(g, frequencies, frequency_threshold)
        add_tunneling_edges!(g, frequencies, tunneling_threshold, width, height)

        mobile_sites = connected_to_edge(g, width, height)
        n_inner_mobile = length(mobile_sites) - 2 * width - 2 * height + 4
        n_i_m_ratio = n_inner_mobile / ((width - 2) * (height - 2))
        push!(ratios, n_i_m_ratio)
    end

    println("Average inner mobile site ratio: ", round(mean(ratios), digits = 2))
    return pl
end


function remove_edges!(g, frequencies, threshold)
    for site in vertices(g)
        for neighbor in copy(neighbors(g, site))
            if abs(frequencies[site] - frequencies[neighbor]) > threshold
                rem_edge!(g, site, neighbor)
            end
        end
    end
end


function add_tunneling_edges!(g, frequencies, threshold, width, height)
    for y in 1:height
        for x in 1:width
            site = x + (y - 1) * width
            for next_point in [(-2, 0), (-1, 1), (-1, -1), (0, -2), (0, 2), (1, 1), (1, -1), (2, 0)]
                new_x = x + next_point[1]
                new_y = y + next_point[2]
                if is_point_inside(new_x, new_y, width, height)
                    new_site = new_x + (new_y - 1) * width
                    if abs(frequencies[site] - frequencies[new_site]) < threshold
                        add_edge!(g, site, new_site)
                    end
                end
            end
        end
    end
end


function is_point_inside(x, y, width, height)
    return x >= 1 && x <= width && y >= 1 && y <= height
end


function connected_to_edge(g, width, height)
    visited = Set{Int}()

    # BFS queue initialized with outer nodes
    outer_nodes = Set(1:width 
    ) ∪ Set(
        (height-1)*width+1 : width*height  # top row
    ) ∪ Set(
        1:width:width*(height-1)+1        # left column
    ) ∪ Set(
        width:width:width*height           # right column
    )
    queue = collect(outer_nodes)

    while !isempty(queue)
        v = popfirst!(queue)
        if v ∉ visited
            push!(visited, v)
            for n in neighbors(g, v)
                if n ∉ visited
                    push!(queue, n)
                end
            end
        end
    end

    return visited # this is the set of nodes connected to the edge
end

f()