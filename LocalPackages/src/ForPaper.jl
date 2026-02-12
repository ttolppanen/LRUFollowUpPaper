export tunable_coupler_paper, qt_tunable_coupler_paper
export tunable_coupler_paper_positive, qt_tunable_coupler_paper_positive
export tunnel_tunable_coupler_paper, qt_tunnel_tunable_coupler_paper
export four_qubit_coupler_paper, qt_four_qubit_coupler_paper
export five_qubit_coupler_paper, qt_five_qubit_coupler_paper
export w_off_position
export logical_state_index, zz_coupling_strength, zz_zero_point


# In Tunable Coupling Architecture for Fixed-Frequency Transmon Superconducting Qubits
# They show that for the setup where the coupler is below in frequency, the direct coupling
# and one of the coupligns to coupler have negative sign. But they don't use Bose-Hubbard model,
# in the BH model this should translate to g_2c being negative.
function tunable_coupler_paper(d, w, U, g)
    return bosehubbard(d, 3, [(1, 2, g[1]), (2, 3, -g[2]), (1, 3, g[3])]; w, U)
end
function qt_tunable_coupler_paper(d, w, U, g)
    return QuantumObject(tunable_coupler_paper(d, w, U, g), dims = (d, d, d))
end

function tunable_coupler_paper_positive(d, w, U, g)
    return bosehubbard(d, 3, [(1, 2, g[1]), (2, 3, g[2]), (1, 3, g[3])]; w, U)
end
function qt_tunable_coupler_paper_positive(d, w, U, g)
    return QuantumObject(tunable_coupler_paper_positive(d, w, U, g), dims = (d, d, d))
end

function tunnel_tunable_coupler_paper(d, w, U, g)
    return bosehubbard(d, 5, [(1, 2, g[1]), (2, 3, -g[2]), (3, 4, -g[3]), (4, 5, g[4]), (1, 3, g[5]), (3, 5, g[5])]; w, U)
end
function qt_tunnel_tunable_coupler_paper(d, w, U, g; kwargs...)
    return QuantumObject(tunnel_tunable_coupler_paper(d, w, U, g; kwargs...), dims = (d, d, d, d, d))
end

function four_qubit_coupler_paper(d, w, U, g)
    return bosehubbard(d, 7, [(1, 2, g[1]), (2, 3, -g[1]), (3, 4, -g[1]), (4, 5, g[1]), (5, 6, g[1]), (6, 7, -g[1]), (1, 3, g[2]), (3, 5, g[2]), (5, 7, g[2])]; w, U)
end
function qt_four_qubit_coupler_paper(d, w, U, g; kwargs...)
    return QuantumObject(four_qubit_coupler_paper(d, w, U, g; kwargs...), dims = (d, d, d, d, d, d, d))
end

function five_qubit_coupler_paper(d, w, U, g)
    return bosehubbard(d, 9, [(1, 2, g[1]), (2, 3, -g[1]), (3, 4, -g[1]), (4, 5, g[1]), (5, 6, g[1]), (6, 7, -g[1]), (7, 8, -g[1]), (8, 9, g[1]), (1, 3, g[2]), (3, 5, g[2]), (5, 7, g[2]), (7, 9, g[2])]; w, U)
end
function qt_five_qubit_coupler_paper(d, w, U, g; kwargs...)
    return QuantumObject(five_qubit_coupler_paper(d, w, U, g; kwargs...), dims = (d, d, d, d, d, d, d, d, d))
end

# WEIRD +- G12!
function w_off_position(w1, w2, g1c, g2c, g12; is_above = false)
    G = 2 * g12 / (g1c * g2c)
    if is_above
        return (w1 + w2) / 2 + 1 / G + sqrt(1 / G^2 + ((w1 - w2) / 2)^2) - g12
    else
        return (w1 + w2) / 2 - 1 / G - sqrt(1 / G^2 + ((w1 - w2) / 2)^2) + g12
    end
end

function logical_state_index(comparison_states, eigen_states, d)
    # state00 = QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 0)))
    # state10 = QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 0)))
    # state01 = QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 1)))
    # state11 = QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 1)))
    out = []
    for bare_state in comparison_states
        max_overlap = 0.0
        max_index = 0
        for (i, eigen_state) in enumerate(eigen_states)
            if i in out
                continue
            end
            overlap = abs(bare_state' * eigen_state)^2
            if overlap > max_overlap
                max_overlap = overlap
                max_index = i
            end
        end
        push!(out, max_index)
    end
    return out[1], out[2], out[3], out[4]
end

function zz_coupling_strength(dressed_states_i, E, d)
    i_00, i_10, i_01, i_11 = dressed_states_i[1], dressed_states_i[2], dressed_states_i[3], dressed_states_i[4]
    return (E[i_11] - E[i_10]) - (E[i_01] - E[i_00])
end

function zz_zero_point(w1, w2, wcs, Us, gs)
    zzs = []
    comparison_states = []
    d = 3
    for wc in wcs
        if isempty(comparison_states)
            push!(comparison_states, QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 0))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 0))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 1))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 1))))
        end
        w = [w1, wc, w2]
        U = Us
        g = gs
        H = qt_tunable_coupler_paper(d, w, U, g)

        E, psi, Uni = eigenstates(H)
        dressed_states_i = logical_state_index(comparison_states, psi, d)
        push!(zzs, zz_coupling_strength(dressed_states_i, E, d))
        comparison_states = psi[collect(dressed_states_i)]
    end
    wc_zeros = []
    prev_zz_sign = sign(zzs[1])
    for (i, zz) in enumerate(zzs)
        if sign(zz) != prev_zz_sign
            push!(wc_zeros, wcs[i])
            prev_zz_sign = sign(zz)
        end
    end
    if isempty(wc_zeros)
        error("No ZZ zero points found")
    end
    return wc_zeros
end