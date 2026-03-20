using Distributed

addprocs(7)

@everywhere begin
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
end

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

@everywhere function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:alg] = "toolbox Vern9"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    
    N = 100
    w2s = collect(range(3900, 4100, length=N)) ./ param[:J12].val
    param[:w2] = ValueInfo(w2s, "(3.702-3.750)GHz")
    wc_range = collect(range(3400, 3800, length=N)) ./ param[:J12].val
    param[:wc] = ValueInfo(wc_range, "(3.600-3.800)GHz")

    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t_below = 0.15
    param[:t_below] = ValueInfo(t_below * 2 * π * param[:J12].val, "$(t_below)μs")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

@everywhere function calc_max_amps(param, state, t_list, e_ops, wc, w2)
    w = [param[:w1].val, wc, w2]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

    H = qt_tunable_coupler_paper(param[:d], w, U, g)

    out = mesolve(H, state, t_list, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9())
    return maximum(real(out.expect[1, :]))
end

@everywhere function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    t_list = LinRange(0, param[:t_below].val, param[:t_points])

    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    n2 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n2]

    state = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))

    param_sweep = [(wc, w2) for w2 in param[:w2].val for wc in param[:wc].val]
    max_amps = pmap(p -> calc_max_amps(param, state, t_list, e_ops, p[1], p[2]), param_sweep)

    max_amps = reshape(max_amps, length(param[:wc].val), length(param[:w2].val))

    x = param[:w2].val .* 5
    y = param[:wc].val .* 5
    
    heatmap(x, y, max_amps, xlabel = "w2 (MHz)", ylabel = "wc (MHz)", colorbar_title = "max(<n2>)")
    savefig(joinpath(@__DIR__, "heatmap.png"))
end

@everywhere BLAS.set_num_threads(1)

@time f()

rmprocs(workers())
