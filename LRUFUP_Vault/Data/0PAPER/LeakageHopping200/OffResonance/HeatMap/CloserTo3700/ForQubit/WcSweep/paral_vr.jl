using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays
using LocalPackages
using QuantumToolbox
using QuantumTimes
using OrdinaryDiffEq

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:alg] = "toolbox Vern9"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(0 * 50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    
    N = 100
    w2s = collect(range(3700, 4000, length=N)) ./ param[:J12].val
    param[:w2] = ValueInfo(w2s, "(3.702-3.750)GHz")
    wc_range = collect(range(3300, 3500, length=N)) ./ param[:J12].val
    param[:wc] = ValueInfo(wc_range, "(3.600-3.800)GHz")

    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t_below = 0.1
    param[:t_below] = ValueInfo(t_below * 2 * π * param[:J12].val, "$(t_below)μs")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    t_list_below = LinRange(0, param[:t_below].val, param[:t_points])

    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    n2 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n2]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    state = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))


    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

    H = qt_tunable_coupler_paper_param(d, U, g)
    w1 = [param[:w1].val]
    wc = param[:wc].val
    w2 = param[:w2].val

    out = mesolve_map(H, state, t_list_below; params = (w1, wc, w2), e_ops = e_ops, progress_bar= false, abstol = param[:abstol], reltol = param[:reltol], alg = Vern9(), ensemblealg = EnsembleSplitThreads())
 
    max_amps = zeros(length(param[:w2].val), length(param[:wc].val))
    for I in CartesianIndices(max_amps)
        w2_i, wc_i = Tuple(I)
        max_amps[wc_i, w2_i] = maximum(real(out[1, 1, wc_i, w2_i].expect[1, :]))
    end

    x = param[:w2].val .* 5
    y = param[:wc].val .* 5
    
    heatmap(x, y, max_amps, xlabel = "w2 (MHz)", ylabel = "wc (MHz)", colorbar_title = "max(<n2>)")
    savefig(joinpath(@__DIR__, "heatmap.png"))
end

default_threads = BLAS.get_num_threads()
BLAS.set_num_threads(1)

@time f()

BLAS.set_num_threads(default_threads)