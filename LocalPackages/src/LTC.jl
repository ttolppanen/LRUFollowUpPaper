export lossy_tunable_coupler, qt_lossy_tunable_coupler

function lossy_tunable_coupler(d, w, U, g)
    return bosehubbard(d, 4, [(1, 2, g[1]), (2, 3, -g[2]), (1, 3, g[3]), (2, 4, g[4])]; w, U)
end
function qt_lossy_tunable_coupler(d, w, U, g)
    return QuantumObject(lossy_tunable_coupler(d, w, U, g), dims = (d, d, d, d))
end