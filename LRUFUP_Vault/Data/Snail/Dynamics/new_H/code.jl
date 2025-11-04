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

# this is calculated so that |2> -> |1> is in resonance
function Ec2(Ej1, Ej2, Ec1)
    c = 2 * (sqrt(8 * Ej1 * Ec1) - Ec1) - Ec1
    Ec2 = (sqrt(32 * Ej2 + 4 * c) - 2 * sqrt(8 * Ej2)) / 2
    return Ec2
end

function two_transmon_snail(d, phi_c, E_c, E_j, N_junc)
    beta = E_j[3] ./ E_j[1:2]
    k = (phi_c <= π/2 && phi_c >= -π/2) ? 1.0 : 0.0
    alpha = 1 / sqrt(sin(phi_c)^2 + N_junc^6 * cos(phi_c)^2)
    K2 = alpha * (cos(phi_c) + (-1)^k * N_junc^2 * abs(cos(phi_c)))
    K3 = -alpha * sin(phi_c) * (1 - 1 / N_junc^2)
    println("alpha = $alpha \n K2 = $K2 \n K3 = $K3")

    phi_zpf(i) = ((2 * E_c[i]) / (E_j[i] * (1 + beta[i] * K2)))^(1/4)
    n_zpf(i) = ((E_j[i] * (1 + beta[i] * K2)) / (32 * E_c[i]))^(1/4)
    c_op = create(d)
    d_op = destroy(d)

    n_1 = -1im * n_zpf(1) * (d_op - c_op)
    phi_1 = phi_zpf(1) * (d_op + c_op)

    H1 = 4 * E_c[1] * n_1^2 - E_j[1] * (cos(phi_1) - beta[1] * K2 / 2 * phi_1^2 - beta[1] * K3 / 6 * phi_1^3)
    λ, psi, U1 = eigenstates(H1)
    U1 = QuantumObject(U1, dims = d)
    
    n_2 = -1im * n_zpf(2) * (d_op - c_op)
    phi_2 = phi_zpf(2) * (d_op + c_op)
    H2 = 4 * E_c[2] * n_2^2 - E_j[2] * (cos(phi_2) - beta[2] * K2 / 2 * phi_2^2 + beta[2] * K3 / 6 * phi_2^3)
    λ, psi, U2 = eigenstates(H2)
    U2 = QuantumObject(U2, dims = d)

    U = tensor(U1, U2)
    coupling = -E_j[3] * (K2 * tensor(phi_1, phi_2) - K3 / 2 * (tensor(phi_1, phi_2^2) - tensor(phi_1^2, phi_2)))
    return U' * (tensor(H1, H2) + coupling) * U
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 11
    param[:L] = 2
    N = 2
    state = QuantumObject(tensor(basis(param[:d], N), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|$(N)0>")
    param[:alg] = "toolbox"
    param[:H] = "snail"
    param[:phi_c] = ValueInfo(-π * 2 / 3, "-2π/3")
    Ej1 = 20
    param[:Ej1] = ValueInfo(Ej1, "$(Ej1)GHz")
    param[:Ej2] = ValueInfo(15 / 2, "15/2*EJ1")
    param[:Ec1] = ValueInfo(0.3 / param[:Ej1].val, "300MHz")
    param[:Ec2] = ValueInfo(0.3 / param[:Ej1].val, "300MHz")
    param[:Ej_snail] = ValueInfo(3 / 2 , "3/2*EJ1")
    param[:N_junc] = 2
    t = 10000
    param[:t] = ValueInfo(t * 2 * π * 10^-3 * param[:Ej1].val, "$(t)ns")
    param[:diss_rate] = 0

    param[:Ej1] = ValueInfo(1, param[:Ej1].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    state = param[:state].val

    n1 = qt_single_site_op(d, L, num(d), 1)
    n_all = qt_n_all(d, L)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n_all]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    E_c = [param[:Ec1].val, param[:Ec2].val]
    E_j = [param[:Ej1].val, param[:Ej2].val, param[:Ej_snail].val]
    H = two_transmon_snail(d, param[:phi_c].val, E_c, E_j, param[:N_junc])
    # drive_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    # drive_op = drive_op / 2
    # pulse = gaussian_pulse(param[:gate_time].val, π)
    # Ht = H + QobjEvo(drive_op, (p, t) -> pulse(t))

    # t_list = LinRange(0, param[:t].val, param[:t_points])
    # out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))
    # push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n_all"))

    # path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f();