export tunable_coupler_paper, qt_tunable_coupler_paper
export tunable_coupler_paper_positive, qt_tunable_coupler_paper_positive


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