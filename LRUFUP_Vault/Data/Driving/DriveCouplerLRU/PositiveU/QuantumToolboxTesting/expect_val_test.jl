using QuantumToolbox

function f()
    H = 0.5 * sigmay()
    ψ0 = basis(2, 0)
    e_ops = [
        proj(basis(2, 0)),
        proj(basis(2, 1)),
        basis(2, 0) * basis(2, 1)'
    ]
    tlist = LinRange(0, 10, 100)
    sol = sesolve(H, ψ0, tlist; saveat=[tlist[end]], progress_bar = Val(false))
    print(sol.states[1])
end

f()