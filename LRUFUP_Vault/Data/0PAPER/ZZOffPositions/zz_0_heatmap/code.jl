using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using OrderedCollections
using SparseArrays
using LocalPackages
using QuantumToolbox
using JLD2
using HDF5

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3

    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1s] = ValueInfo(collect(range(3702, 4298, length = 100)) ./ param[:J12].val, "3.702 - 4.298GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wc_min] = ValueInfo(2500 / param[:J12].val, "2.5GHz")
    param[:wc_max] = ValueInfo(4000 / param[:J12].val, "4")
    param[:Us] = ValueInfo(collect(range(200, 300, length = 10)) ./ param[:J12].val, "250 - 300MHz")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()

    gs = [param[:Jc].val, param[:Jc].val, param[:J12].val]

    x_N = length(param[:w1s].val)
    y_N = length(param[:Us].val)
    wcs = zeros(y_N, x_N)
    for I in CartesianIndices((y_N, x_N))
        i, j = Tuple(I)
        w1 = param[:w1s].val[j]
        U = param[:Us].val[i]
            
        Us = [U, U, U]
        if abs(w1 - 4000 / 5) > U
            wcs[i, j] = 0
            continue
        end
        try
            wc_pos, zz_val = zz_0_with_roots(w1, param[:w2].val, param[:wc_min].val, param[:wc_max].val, Us, gs)
            wcs[i, j] = wc_pos[1]
        catch
            wcs[i, j] = 0
        end
    end
    x = param[:w1s].val .* 5
    y = param[:Us].val .* 5
    # t = t_list ./ (2 * π * 5)
    
    z = wcs .* 5
    p = param
    @save joinpath(@__DIR__, "data.jld2") x y z p
end

function make_heatmap()
    @load joinpath(@__DIR__, "data.jld2") x y z
    heatmap(x, y, z, xlabel = "w1 (MHz)", ylabel = "U (MHz)", colorbar_title = "wc", clims = (3000, 3500))
    savefig(joinpath(@__DIR__, "heatmap.png"))
end

function data_to_h5()
    filename = joinpath(@__DIR__, "data.h5")
    h5open(filename, "w") do file
        @load joinpath(@__DIR__, "data.jld2") x y z p

        write(file, "x", collect(x))
        write(file, "y", collect(y))
        write(file, "z", collect(z))
    end
    println("Data saved!")
end

# @time f()
# make_heatmap()
data_to_h5()