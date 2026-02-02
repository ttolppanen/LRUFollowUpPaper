export tunable_coupler_paper, qt_tunable_coupler_paper


# In Tunable Coupling Architecture for Fixed-Frequency Transmon Superconducting Qubits
# They show that for the setup where the coupler is below in frequency, the direct coupling
# and one of the coupligns to coupler have negative sign. (?)
function tunable_coupler_paper(w, U, g)
    return bosehubbard(3, 3, [(1, 2, -g[1]), (2, 3, g[2]), (1, 3, -g[3])]; w, U)
end
function qt_tunable_coupler_paper(w, U, g)
    return QuantumObject(tunable_coupler_paper(w, U, g), dims = (3, 3, 3))
end