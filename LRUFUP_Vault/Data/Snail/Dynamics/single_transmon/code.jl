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

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function single_transmon_test(d, phi_c, E_c, E_j)
    phi_zpf = (2 * E_c / E_j)^(1/4)
    n_zpf = (E_j / (32 * E_c))^(1/4)
    c_op = create(d)
    d_op = destroy(d)
    n = -1im * n_zpf * (d_op - c_op)
    ec_part = 4 * E_c * n * n
    phi_1 = phi_zpf * (d_op + c_op)
    ej_part = -1 * (E_j * cos(phi_1 - phi_c))
    return ec_part + ej_part
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 1
    N = 0
    state = basis(param[:d], N)
    param[:state] = ValueInfo(state, "|$(N)>")
    param[:alg] = "toolbox"
    param[:H] = "transmon"
    param[:phi_c] = ValueInfo(-π * 3 / 2, "-3π/2")
    Ej = 20
    param[:Ej] = ValueInfo(Ej, "$(Ej)GHz")
    param[:Ec] = ValueInfo(0.3 / param[:Ej].val, "300MHz")
    t = 100
    param[:t] = ValueInfo(t * 2 * π * 10^-3 * param[:Ej].val, "$(t)ns")
    param[:diss_rate] = 0

    param[:Ej] = ValueInfo(1, param[:Ej].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d];
    state = param[:state].val

    n1 = num(d)
    e_ops = [n1]
    
    H = single_transmon_test(d, param[:phi_c].val, param[:Ec].val, param[:Ej].val)
    λ, psi, V = eigenstates(H)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    return λ
end


@time a = f()