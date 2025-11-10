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

function two_transmon_snail(d, phi_c, E_c, E_j, N_junc; D = 15)
    if d > D
        D = d
    end
    beta = E_j[3] ./ E_j[1:2]
    k = (phi_c <= π/2 && phi_c >= -π/2) ? 1.0 : 0.0
    alpha = 1 / sqrt(sin(phi_c)^2 + N_junc^6 * cos(phi_c)^2)
    K2 = alpha * (cos(phi_c) + (-1)^k * N_junc^2 * abs(cos(phi_c)))
    K3 = -alpha * sin(phi_c) * (1 - 1 / N_junc^2)
    # println("alpha = $alpha \n K2 = $K2 \n K3 = $K3")

    phi_zpf(i) = ((2 * E_c[i]) / (E_j[i] * (1 + beta[i] * K2)))^(1/4)
    n_zpf(i) = ((E_j[i] * (1 + beta[i] * K2)) / (32 * E_c[i]))^(1/4)
    c_op = create(D)
    d_op = destroy(D)

    n_1 = -1im * n_zpf(1) * (d_op - c_op)
    phi_1 = phi_zpf(1) * (d_op + c_op)

    H1 = 4 * E_c[1] * n_1^2 - E_j[1] * (cos(phi_1) - beta[1] * K2 / 2 * phi_1^2 - beta[1] * K3 / 6 * phi_1^3)
    λ, psi, U1 = eigenstates(H1)
    U1 = QuantumObject(U1, dims = D)
    println(λ[1:3])
    
    n_2 = -1im * n_zpf(2) * (d_op - c_op)
    phi_2 = phi_zpf(2) * (d_op + c_op)
    H2 = 4 * E_c[2] * n_2^2 - E_j[2] * (cos(phi_2) - beta[2] * K2 / 2 * phi_2^2 + beta[2] * K3 / 6 * phi_2^3)
    λ, psi, U2 = eigenstates(H2)
    U2 = QuantumObject(U2, dims = D)
    println(λ[1:3])

    U = tensor(U1, U2)
    coupling = -E_j[3] * (K2 * tensor(phi_1, phi_2) - K3 / 2 * (tensor(phi_1, phi_2^2) - tensor(phi_1^2, phi_2)))
    H = U' * (tensor(H1, qeye(D)) + tensor(qeye(D), H2) + coupling) * U
    
    H_row = []
    for i in 0:(d-1)
        i_ind = (i * D + 1):(i * D + d)
        push!(H_row, H[(1:d), i_ind])
    end

    H_out = hcat(H_row...)
    for j in 1:(d-1)
        j_ind = (j * D + 1):(j * D + d)
        H_row = []
        for i in 0:(d-1)
            i_ind = (i * D + 1):(i * D + d)
            push!(H_row, H[j_ind, i_ind])
        end
        H_out = vcat(H_out, hcat(H_row...))
    end
    
    # first = 1:d # indeces of differenct blocks, look at tensor product for matrix to make sense hehe
    # second = (D + 1):(D + d) # length is d, each block is between of D, that is the big dimension for accuracy
    # third = (2 * D + 1):(2 * D + d)
    # H = [H[first, first] H[first, second] H[first, third]; 
    #     H[second, first] H[second, second] H[second, third];
    #     H[third, first] H[third, second] H[third, third]] # cutted dimensions
    H = QuantumObject(H_out, dims = (d, d))
    λ, psi, U2 = eigenstates(H)
    E = labeled_eigenvalue_dict(λ, d)
    println("00 -> 20 = $(E["20"] - E["00"])")
    # println("00 -> 10 = $(E["10"] - E["00"])")
    println("00 -> 01 = $(E["01"] - E["00"])")
    return QuantumObject(H_out, dims = (d, d))
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 2
    N = 2
    state = QuantumObject(tensor(basis(param[:d], N), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|$(N)0>")
    param[:alg] = "toolbox"
    param[:H] = "snail"
    param[:phi_c] = ValueInfo(-π * 2 / 3, "-2π/3")
    Ej1 = 2.5
    param[:Ej1] = ValueInfo(Ej1, "$(Ej1)GHz")
    param[:Ej2] = ValueInfo(7.1, "7.1 * EJ1")
    param[:Ec1] = ValueInfo(0.3 / param[:Ej1].val, "300MHz")
    param[:Ec2] = ValueInfo(0.195 / param[:Ej1].val, "195MHz")
    param[:Ej_snail] = ValueInfo(3 / 2 , "3/2*EJ1")
    param[:N_junc] = 2
    t = 10
    param[:t] = ValueInfo(t * 2 * π * param[:Ej1].val, "$(t)ns")
    param[:diss_rate] = 0 * 30 * 5 * 10^-3 / Ej1

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
    n2 = qt_single_site_op(d, L, num(d), 2)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n2]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    E_c = [param[:Ec1].val, param[:Ec2].val]
    E_j = [param[:Ej1].val, param[:Ej2].val, param[:Ej_snail].val]
    H = two_transmon_snail(d, param[:phi_c].val, E_c, E_j, param[:N_junc])

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))
    push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n2"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)

    # drive_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    # drive_op = drive_op / 2
    # pulse = gaussian_pulse(param[:gate_time].val, π)
    # Ht = H + QobjEvo(drive_op, (p, t) -> pulse(t))
end


@time f();