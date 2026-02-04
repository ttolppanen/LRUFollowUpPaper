export tunable_coupler_paper, qt_tunable_coupler_paper
export tunable_coupler_paper_positive, qt_tunable_coupler_paper_positive
export w_off_position


# In Tunable Coupling Architecture for Fixed-Frequency Transmon Superconducting Qubits
# They show that for the setup where the coupler is below in frequency, the direct coupling
# and one of the coupligns to coupler have negative sign. But they don't use Bose-Hubbard model,
# in the BH model this should translate to g_2c being negative.
function tunable_coupler_paper(w, U, g)
    return bosehubbard(3, 3, [(1, 2, g[1]), (2, 3, -g[2]), (1, 3, g[3])]; w, U)
end
function qt_tunable_coupler_paper(w, U, g)
    return QuantumObject(tunable_coupler_paper(w, U, g), dims = (3, 3, 3))
end

function tunable_coupler_paper_positive(w, U, g)
    return bosehubbard(3, 3, [(1, 2, g[1]), (2, 3, g[2]), (1, 3, g[3])]; w, U)
end
function qt_tunable_coupler_paper_positive(w, U, g)
    return QuantumObject(tunable_coupler_paper_positive(w, U, g), dims = (3, 3, 3))
end

function w_off_position(w1, w2, g1c, g2c, g12; is_above = false)
    G = 2 * g12 / (g1c * g2c)
    if is_above
        return (w1 + w2) / 2 + 1 / G + sqrt(1 / G^2 + ((w1 - w2) / 2)^2) - 1
    else
        return (w1 + w2) / 2 - 1 / G - sqrt(1 / G^2 + ((w1 - w2) / 2)^2) + 1
    end
end