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

function mesh_utility(xstart::Real, xend::Real, nelem::Int, porder::Int,z::Union{AbstractArray,AbstractMatrix,AbstractVector})
    dx = (xend - xstart) / nelem
    ptspe = porder + 1
    T = eltype(xend)
    xcg = zeros(T, ptspe * nelem)
    for i in 1:nelem
        x0 = xstart + (i-1)*dx
        for j in 1:ptspe
            xcg[(i-1)*ptspe+j] = x0 + (1+z[j])*dx/2
        end
    end
    e2vcg = reshape(xcg, porder + 1, nelem)
    xcg = unique(xcg)
    return (xcg=xcg, e2vcg=e2vcg)
end
