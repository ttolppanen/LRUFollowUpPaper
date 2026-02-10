# using Plots
# using LinearAlgebra
# using QuantumStates
# using QuantumOperators
# using PlotAndSave
# using QuantumTimeEvolution
# using OrderedCollections
# using SparseArrays
using LocalPackages
# using QuantumToolbox
# using QuantumTimes

function initial_parameters()
    param = OrderedDict()

    param[:J12] = 5
    param[:Jc] = 50 / param[:J12]
    param[:w1] = 4000 / param[:J12]
    param[:w2] = 4000 / param[:J12]
    param[:U] = 300 / param[:J12]


    U = [param[:U], param[:U], param[:U]]
    gs = [param[:Jc], param[:Jc], 1]
    wcs = collect(range(3000, 3800, length=10000)) ./ param[:J12]
    wcs = zz_zero_point(param[:w1], param[:w2], wcs, U, gs)
    display(wcs .* param[:J12])
end

initial_parameters()