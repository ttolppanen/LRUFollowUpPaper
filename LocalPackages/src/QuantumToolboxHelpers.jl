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
function tunable_coupler(w, U, g)
    return bosehubbard(3, 3, [(1, 2, -g[1]), (2, 3, g[2]), (1, 3, g[3])]; w, U)
end
function qt_tunable_coupler(w, U, g)
    return QuantumObject(tunable_coupler(w, U, g), dims = (3, 3, 3))
end
