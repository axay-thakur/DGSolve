# --- 7. mesh_utility.m ---
export mesh_utility
function mesh_utility(xstart::Real, xend::Real, nelem::Int, porder::Int)
    xcg = collect(range(xstart, stop=xend, length=porder * nelem + 1))
    doubling_idx = 1:porder:(porder * nelem + 1)
    indexing = ones(Int, length(xcg))
    for idx in doubling_idx[2:end-1]
        indexing[idx] = 2
    end
    xrep = vcat([fill(xcg[i], indexing[i]) for i in eachindex(xcg)]...)
    e2vcg = reshape(xrep, porder + 1, nelem)
    return (xcg=xcg, e2vcg=e2vcg)
end