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

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end


function two_transmon_snail_with_k(d, K2, K3, E_c, E_j, phi; D = 15)
    if d > D
        D = d
    end
    beta = E_j[1:2] ./ E_j[3]

    phi_zpf(i) = ((2 * E_c[i]) / (E_j[3] * (beta[i] * cos(phi[i]) + K2)))^(1/4)
    n_zpf(i) = ((E_j[3] * (beta[i] * cos(phi[i]) + K2)) / (32 * E_c[i]))^(1/4)
    c_op = create(D)
    d_op = destroy(D)

    n_1 = -1im * n_zpf(1) * (d_op - c_op)
    phi_1 = phi_zpf(1) * (d_op + c_op)

    H1 = 4 * E_c[1] * n_1^2 + E_j[3] / 2 * (beta[1] * cos(phi[1]) + K2) * phi_1^2 - E_j[3] / 6 * (beta[1] * sin(phi[1]) - K3) * phi_1^3 - E_j[3]/24 * beta[1]*cos(phi[1]) * phi_1^4
    λ, psi, U1 = eigenstates(H1)
    U1 = QuantumObject(U1, dims = D)
    # println(λ[1:3])
    # println((λ[1:3] .- λ[1]))
    # a1 = λ[3] - λ[1]
    
    n_2 = -1im * n_zpf(2) * (d_op - c_op)
    phi_2 = phi_zpf(2) * (d_op + c_op)
    H2 = 4 * E_c[2] * n_2^2 + E_j[3] / 2 * (beta[2] * cos(phi[2]) + K2) * phi_2^2 - E_j[3] / 6 * (beta[2] * sin(phi[2]) + K3) * phi_2^3 - E_j[3]/24 * beta[2]*cos(phi[2]) * phi_2^4
    λ, psi, U2 = eigenstates(H2)
    U2 = QuantumObject(U2, dims = D)
    # println(λ[1:3])
    # println((λ[1:3] .- λ[1]))
    # a2 = λ[2] - λ[1]
    # println("Resonance: $(a1 - a2)")

    U = tensor(U1, U2)
    coupling = 8 * E_c[3] * tensor(n_1, n_2) - E_j[3] * K2 * tensor(phi_1, phi_2) - E_j[3]/2 * K3 * (tensor(phi_1^2, phi_2) - tensor(phi_1, phi_2^2))
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
    
    H = QuantumObject(H_out, dims = (d, d))
    # λ, psi, U2 = eigenstates(H)
    # E = labeled_eigenvalue_dict(λ, d)
    # println("00 -> 20 = $(E["20"] - E["00"])")
    # # println("00 -> 10 = $(E["10"] - E["00"])")
    # println("00 -> 01 = $(E["01"] - E["00"])")
    return QuantumObject(H_out, dims = (d, d))
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 2
    N = 1
    state = QuantumObject(tensor(basis(param[:d], N), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|$(N)0>")
    param[:alg] = "toolbox"
    param[:H] = "snail"

    Ejs = 20
    param[:Ejs] = ValueInfo(1, "$(Ejs)GHz")
    param[:Units] = ValueInfo(Ejs, "Ejs")

    param[:Ej1] = ValueInfo(0.5, "b_1 = 0.5 -> Ej1 = $(0.5 * Ejs)GHz")
    param[:Ej2] = ValueInfo(4, "b_2 = 4 -> Ej2 = $(4 * Ejs)GHz")
    param[:Ec1] = ValueInfo(0.15 / param[:Units].val, "150MHz")
    param[:Ec2] = ValueInfo(0.082 / param[:Units].val, "82MHz")
    param[:Ec12] = ValueInfo(0.0087 / param[:Units].val, "8.7MHz")
    Ec12s = collect(range(0.0001, 0.01, length = 1000)) ./ param[:Units].val
    param[:Ec12] = ValueInfo(Ec12s, "0.1 - 10MHz")

    param[:K2] = 0.123
    param[:K3] = -0.177
    param[:phi1] = -0.089
    param[:phi2] = 0.011

    t = 10
    param[:t] = ValueInfo(t * 2 * π * param[:Units].val, "$(t)ns")
    # param[:diss_rate] = 0.01
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------
    lines = []
    d = param[:d]; L = param[:L]
    t_list = LinRange(0, param[:t].val, param[:t_points])
    state = param[:state].val

    n1 = qt_single_site_op(d, L, num(d), 1)
    n2 = qt_single_site_op(d, L, num(d), 2)
    anh1 = qt_single_site_op(d, L, qt_anh_op(d), 1)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n2, anh1]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    
    E_j = [param[:Ej1].val, param[:Ej2].val, param[:Ejs].val]
    phi = [param[:phi1], param[:phi2]]
    

    max_vals = []
    for ec12 in param[:Ec12].val
        E_c = [param[:Ec1].val, param[:Ec2].val, ec12]
        H = two_transmon_snail_with_k(d, param[:K2], param[:K3], E_c, E_j, phi)
        out = mesolve(H, state, t_list, e_ops = e_ops; abstol = 1e-14, reltol = 1e-8, maxiters = 10^7, progress_bar = false)
        push!(max_vals, maximum(real(out.expect[4, :])))
    end
    lines = [LineInfo(param[:Ec12].val .* param[:Units].val, max_vals, 1, "max(<n2>), min=$(round(minimum(max_vals), digits = 3))")]

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "Ec12", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)

    # drive_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    # drive_op = drive_op / 2
    # pulse = gaussian_pulse(param[:gate_time].val, π)
    # Ht = H + QobjEvo(drive_op, (p, t) -> pulse(t))
    # pl = plot(lines[1].x, lines[1].y, ylim = (0, 1.1))
    # plot!(pl, lines[1].x, t -> exp(-t / (2 * π * 20.1) / 30))
    # t_time_from_n(t_list, real(out.expect[1, :]), 10 * 1e9)
end

@time f()