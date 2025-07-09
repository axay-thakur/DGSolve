export element_mass_matrix, inv_mass_matrix
# --- 3. element_mass_matrix.m ---
function element_mass_matrix(we::AbstractVector, grad_g::AbstractVector, Tv::AbstractArray)
    nq = size(Tv, 4)
    w = we .* grad_g
    n_rows_Tv_slice = size(Tv, 1)
    Melem = zeros(Float64, n_rows_Tv_slice, n_rows_Tv_slice)
    for i in 1:nq
        Tv_slice_1 = Tv[:, :, 1, i]
        Melem += w[i] * Tv_slice_1 * transpose(Tv_slice_1)
    end
    return Melem
end

# --- 6. inv_mass_matrix.m ---
function inv_mass_matrix(w::AbstractVector, Gj::AbstractMatrix, Tv::AbstractArray)
    nelem, _ = size(Gj)
    ndof = size(Tv, 1)
    Minv = zeros(Float64, ndof, ndof, nelem)
    for i in 1:nelem
        Melem_i = element_mass_matrix(w, Gj[i,:], Tv)
        Minv[:, :, i] .= inv(Melem_i)
    end
    return Minv
end
