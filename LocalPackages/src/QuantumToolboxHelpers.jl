# using QuantumToolbox
# using QuantumOperators

export qt_single_site_op
export qt_anh_op
export qt_op_all
export qt_anh_all
export qt_n_all
export qt_bosehubbard
export tunable_coupler
export qt_tunable_coupler
export tunable_coupler_below
export qt_tunable_coupler_below
export tunable_coupler_plus_coupler
export qt_tunable_coupler_plus_coupler
export tunnel_tunable_coupler
export qt_tunnel_tunable_coupler
export tunnel_tunable_coupler_below
export qt_tunnel_tunable_coupler_below
export cross_tunable_coupler
export qt_cross_tunable_coupler
export qt_two_transmon_snail
export qt_two_transmon_snail_test
export qt_single_transmon_snail
export print_labeled_eigenvalues
export labeled_eigenvalue_dict


function qt_single_site_op(d::Integer, L::Integer, op::QuantumObject, i::Integer)
    ops = []
    for site in 1:L
        if site != i
            push!(ops, qeye(d))
        else
            push!(ops, op)
        end
    end
    return tensor(ops...)
end

function qt_anh_op(d::Integer)
    return num(d) * (num(d) - qeye(d))
end

function qt_op_all(d::Integer, L::Integer, op::QuantumObject)
    all_op = qt_single_site_op(d, L, op, 1)
    for site in 2:L
        all_op += qt_single_site_op(d, L, op, site)
    end
    return all_op
end

function qt_anh_all(d::Integer, L::Integer)
    anh_op = qt_anh_op(d)
    return qt_op_all(d, L, anh_op)
end

function qt_n_all(d::Integer, L::Integer)
    return qt_op_all(d, L, num(d))
end

function qt_bosehubbard(d::Integer, L::Integer; kwargs...) # Order might be reversed here? The last site could be the first?
    return QuantumObject(bosehubbard(d, L; kwargs...), dims = tuple([d for _ in 1:L]...))
end
# coupling between not neighboring sites
function qt_bosehubbard(d::Integer, L::Integer, couplings::Vector{Tuple{Int64, Int64, T}}; kwargs...) where T <: Real # Order might be reversed here? The last site could be the first?
    return QuantumObject(bosehubbard(d, L, couplings; kwargs...), dims = tuple([d for _ in 1:L]...))
end

# 3-site tunable coupler, parameters as lists
# w = [w1, w_c, w2]
# U = [U1, U_c, U2]
# g = [g_c1, g_c2, g_direct]
function tunable_coupler(w, U, g; J_sign = -1)
    return bosehubbard(3, 3, [(1, 2, J_sign * g[1]), (2, 3, g[2]), (1, 3, g[3])]; w, U)
end
function qt_tunable_coupler(w, U, g; kwargs...)
    return QuantumObject(tunable_coupler(w, U, g; kwargs...), dims = (3, 3, 3))
end

# In Tunable Coupling Architecture for Fixed-Frequency Transmon Superconducting Qubits
# they have negative direct coupling
function tunable_coupler_below(w, U, g)
    return bosehubbard(3, 3, [(1, 2, g[1]), (2, 3, g[2]), (1, 3, -g[3])]; w, U)
end
function qt_tunable_coupler_below(w, U, g)
    return QuantumObject(tunable_coupler_below(w, U, g), dims = (3, 3, 3))
end

function tunable_coupler_plus_coupler(w, U, g; J_sign = -1)
    return bosehubbard(3, 4, [(1, 2, g[1]), (2, 3, g[2]), (3, 4, g[3]), (1, 3, J_sign * g[4])]; w, U)
end
function qt_tunable_coupler_plus_coupler(w, U, g; kwargs...)
    return QuantumObject(tunable_coupler_plus_coupler(w, U, g; kwargs...), dims = (3, 3, 3, 3))
end

function tunnel_tunable_coupler(w, U, g; J_sign = -1)
    return bosehubbard(3, 5, [(1, 2, J_sign * g[1]), (2, 3, g[2]), (3, 4, J_sign * g[3]), (4, 5, g[4]), (1, 3, g[5]), (3, 5, g[5])]; w, U)
end
function qt_tunnel_tunable_coupler(w, U, g; kwargs...)
    return QuantumObject(tunnel_tunable_coupler(w, U, g; kwargs...), dims = (3, 3, 3, 3, 3))
end

function tunnel_tunable_coupler_below(w, U, g)
    return bosehubbard(3, 5, [(1, 2, g[1]), (2, 3, g[2]), (3, 4, g[3]), (4, 5, g[4]), (1, 3, -g[5]), (3, 5, -g[5])]; w, U)
end
function qt_tunnel_tunable_coupler_below(w, U, g; kwargs...)
    return QuantumObject(tunnel_tunable_coupler_below(w, U, g; kwargs...), dims = (3, 3, 3, 3, 3))
end


# site structure
#        6
#        7
# 1   2  3  4   5
#        8
#        9

function cross_tunable_coupler(w, U, J_c, J_direct)
    return bosehubbard(3, 9, [
        (1, 2, J_c), (2, 3, J_c), 
        (3, 4, J_c), (4, 5, J_c),
        (6, 7, J_c), (7, 3, J_c),
        (3, 8, J_c), (8, 9, J_c),
        (1, 3, -J_direct), (3, 5, -J_direct),
        (6, 3, -J_direct), (3, 9, -J_direct)
        ]; w, U)
end
function qt_cross_tunable_coupler(w, U, J_c, J_direct)
    return QuantumObject(cross_tunable_coupler(w, U, J_c, J_direct), dims = (3, 3, 3, 3, 3, 3, 3, 3, 3))
end

# phi_c can be for example -2π/3 or -π/2 (for -π/2 the K_2 term vanishes)
# E_c should contain E_c for both transmons
# E_j should contain for all three elements
function qt_two_transmon_snail(d, phi_c, E_c, E_j, N_junc)
    beta = E_j[3] ./ E_j[1:2]
    k = (phi_c <= π/2 && phi_c >= -π/2) ? 0.0 : 1.0
    alpha = 1 / sqrt(sin(phi_c)^2 + N_junc^6 * cos(phi_c)^2)
    K2 = alpha * (cos(phi_c) + (-1)^k * N_junc^2 * abs(cos(phi_c)))
    phi_zpf(i) = ((2 * E_c[i]) / (E_j[i] * (1 + beta[i] * K2)))^(1/4)
    n_zpf(i) = ((E_j[i] * (1 + beta[i] * K2)) / (32 * E_c[i]))^(1/4)
    c_op = [qt_single_site_op(d, 2, create(d), 1), qt_single_site_op(d, 2, create(d), 2)]
    d_op = [qt_single_site_op(d, 2, destroy(d), 1), qt_single_site_op(d, 2, destroy(d), 2)]
    n2_1 = -1im * n_zpf(1) * (d_op[1] - c_op[1]) 
    n2_2 = -1im * n_zpf(2) * (d_op[2] - c_op[2])
    ec_part = 4 * E_c[1] * n2_1 * n2_1 + 4 * E_c[2] * n2_2 * n2_2
    phi_1 = phi_zpf(1) * (d_op[1] + c_op[1])
    phi_2 = phi_zpf(2) * (d_op[2] + c_op[2])
    phi_m = phi_1 - phi_2
    ej_part = -1 * (E_j[1] * cos(phi_1 - phi_c) + E_j[2] * cos(phi_2 + phi_c))
    phi_ext = phi_c + N_junc * ((-1)^k * asin(alpha * sin(phi_c)) + k * π)
    coupling = -1 * (E_j[3] * (alpha * cos(phi_m) - N_junc * cos((phi_ext - phi_m) / N_junc)))
    return ec_part + ej_part + coupling
    # H = ec_part + ej_part
    # λ, psi, U = eigenstates(H)
    # U = QuantumObject(U, dims = (d, d))
    # return U' * (H + coupling) * U # returns as not hermitian?
end

function qt_single_transmon_snail(d, phi_c, E_c, E_j, N_junc)
    beta = E_j[3] ./ E_j[1:2]
    k = (phi_c <= π/2 && phi_c >= -π/2) ? 0.0 : 1.0
    alpha = 1 / sqrt(sin(phi_c)^2 + N_junc^6 * cos(phi_c)^2)
    K2 = alpha * (cos(phi_c) + (-1)^k * N_junc^2 * abs(cos(phi_c)))
    phi_zpf(i) = ((2 * E_c[i]) / (E_j[i] * (1 + beta[i] * K2)))^(1/4)
    n_zpf(i) = ((E_j[i] * (1 + beta[i] * K2)) / (32 * E_c[i]))^(1/4)
    c_op = [qt_single_site_op(d, 1, create(d), 1), qt_single_site_op(d, 1, create(d), 2)]
    d_op = [qt_single_site_op(d, 1, destroy(d), 1), qt_single_site_op(d, 1, destroy(d), 2)]
    n2_1 = -1im * n_zpf(1) * (d_op[1] - c_op[1]) 
    n2_2 = -1im * n_zpf(2) * (d_op[2] - c_op[2])
    ec_part = 4 * E_c[1] * n2_1 * n2_1
    phi_1 = phi_zpf(1) * (d_op[1] + c_op[1])
    phi_2 = phi_zpf(2) * (d_op[2] + c_op[2])
    phi_m = phi_1 - phi_2
    ej_part = -1 * E_j[1] * cos(phi_1 - phi_c)
    H = ec_part + ej_part
    display(H - H')
    λ, psi, U = eigenstates(H)
    U = QuantumObject(U, dims = d)
    return U' * H * U
end

function qt_two_transmon_snail_test(d, phi_c, E_c, E_j, N_junc)
    beta = E_j[3] ./ E_j[1:2]
    k = (phi_c <= π/2 && phi_c >= -π/2) ? 0.0 : 1.0
    alpha = 1 / sqrt(sin(phi_c)^2 + N_junc^6 * cos(phi_c)^2)
    K2 = alpha * (cos(phi_c) + (-1)^k * N_junc^2 * abs(cos(phi_c)))
    phi_zpf(i) = ((2 * E_c[i]) / (E_j[i] * (1 + beta[i] * K2)))^(1/4)
    n_zpf(i) = ((E_j[i] * (1 + beta[i] * K2)) / (32 * E_c[i]))^(1/4)
    c_op = [qt_single_site_op(d, 2, create(d), 1), qt_single_site_op(d, 2, create(d), 2)]
    d_op = [qt_single_site_op(d, 2, destroy(d), 1), qt_single_site_op(d, 2, destroy(d), 2)]
    n2_1 = -1im * n_zpf(1) * (d_op[1] - c_op[1]) 
    n2_2 = -1im * n_zpf(2) * (d_op[2] - c_op[2])
    phi_1 = phi_zpf(1) * (d_op[1] + c_op[1])
    phi_2 = phi_zpf(2) * (d_op[2] + c_op[2])
    H1 = 4 * E_c[1] * n2_1' * n2_1 - E_j[1] * cos(phi_1 - phi_c)
    H2 = 4 * E_c[2] * n2_2' * n2_2 - E_j[2] * cos(phi_2 + phi_c)
    H1, U1 = in_diag_basis(H1)
    H2, U2 = in_diag_basis(H2)
    U = U1 * U2
    H = H1 * H2
    phi_m = phi_1 - phi_2
    phi_ext = phi_c + N_junc * ((-1)^k * asin(alpha * sin(phi_c)) + k * π)
    coupling = -1 * (E_j[3] * (alpha * cos(phi_m) - N_junc * cos((phi_ext - phi_m) / N_junc)))
    return H + U' * (coupling) * U # returns as not hermitian?
end

function in_diag_basis(H::QuantumObject)
    λ, psi, U = eigenstates(H)
    U = QuantumObject(U, dims = H.dims)
    return U' * H * U, U
end

function print_labeled_eigenvalues(λ, d; n = length(λ))
    println("Eigenvalues with computational-basis labels:")
    for idx in 0:(n-1)
        # recover computational-basis indices (i,j)
        i = div(idx, d)
        j = mod(idx, d)
        label = "|$(i)$(j)⟩"
        println(rpad(label, 6), "  ->  ", λ[idx+1])
    end
end

function labeled_eigenvalue_dict(λ, d)
    dict = Dict{String, ComplexF64}()
    N = length(λ)

    for idx in 0:(N-1)
        i = div(idx, d)
        j = mod(idx, d)
        label = "$(i)$(j)"          # or "[$i,$j]" if you prefer
        dict[label] = λ[idx+1]
    end

    return dict
end