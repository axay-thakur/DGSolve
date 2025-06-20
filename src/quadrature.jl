# Helper function for define_quad_onedim_gaussleg
export quad_onedim_gauss_from_recursion, define_quad_onedim_gaussleg
function quad_onedim_gauss_from_recursion(a::AbstractVector, b::AbstractVector, c::AbstractVector, mu0::Real)
    N = length(a)
    d0 = -b ./ a
    d1 = sqrt.(c[2:end] ./ a[1:end-1] ./ a[2:end])
    A = diagm(0 => d0) + diagm(-1 => d1) + diagm(1 => d1)
    F = eigen(A)
    lam = F.values
    v = F.vectors
    z = lam
    w = mu0 .* (v[1, :].^2)
    return w, z
end

# --- 8. define_quad_onedim_gaussleg.m ---
function define_quad_onedim_gaussleg(N::Int)
    mu0 = 2.0
    j_vec = collect(1:N)
    a = (2 .* j_vec .- 1) ./ j_vec
    c = (j_vec .- 1) ./ j_vec
    b = zeros(N)
    w, z = quad_onedim_gauss_from_recursion(a, b, c, mu0)
    return w, z
end