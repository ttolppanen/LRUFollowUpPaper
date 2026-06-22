using PlotAndSave
using JLD2
using HDF5

function save_hdf5()
    filename = joinpath(@__DIR__, "wc_zz0_data.h5")
    w1 = collect(range(3702, 4298, length = 1000))
    wcs = [zz_0_with_w_4GHZ(wi) for wi in w1]
    h5open(filename, "w") do file
        write(file, "w1", w1)
        write(file, "wc", wcs)
    end
    plot(w1, wcs, xlabel = L"\omega / 2\pi~\textrm{(MHz)}", ylabel = L"wc~\textrm{(MHz)}")
end

save_hdf5()