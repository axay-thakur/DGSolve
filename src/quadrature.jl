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

"""
    legendre_poly(p::Int, x::Float64)

Computes Legendre polynomial of degree `p` and its first and second derivatives at `x`.

Returns a tuple: `(L, L1, L2)` where:
- `L` is the value of the polynomial
- `L1` is the first derivative
- `L2` is the second derivative
"""
function legendre_poly_rec(n::Int, x::T) where T
    L0 = 1.0
    L0_1 = 0.0
    L0_2 = 0.0

    L1 = 0.0
    L1_1 = 0.0
    L1_2 = 0.0

    for i in 1:n
        L2, L2_1, L2_2 = L1, L1_1, L1_2
        L1, L1_1, L1_2 = L0, L0_1, L0_2

        a = (2i - 1) / i
        b = (i - 1) / i

        L0    = a * x * L1 - b * L2
        L0_1  = a * (L1 + x * L1_1) - b * L2_1
        L0_2  = a * (2 * L1_1 + x * L1_2) - b * L2_2
    end

    return L0, L0_1, L0_2
end


"""
    legendre_gauss_lobatto(P::Int)

Computes the Legendre-Gauss-Lobatto points and weights for a given number of points `P`.

Returns two vectors:
- `xgl`: nodes (Float64)
- `wgl`: weights (Float64)
"""
function legendre_gauss_lobatto(N::Int)
    n = N - 1
    nh = fld(n + 1, 2)  # floor((p + 1) / 2)

    xgl = zeros(Float64, N)
    wgl = zeros(Float64, N)

    for i in 1:nh
        x = cos((2i - 1) * pi / (2n + 1))
        for _ in 1:20
            L0, L0_1, L0_2 = legendre_poly_rec(n, x)
            denom = -2x * L0_1 + (1 - x^2) * L0_2
            dx = -((1 - x^2) * L0_1) / denom
            x += dx
            if abs(dx) < 1e-20
                break
            end
        end
        L0, _, _ = legendre_poly_rec(n, x)
        j = n + 2 - i
        xgl[j] = x
        wgl[j] = 2 / (n * (n + 1) * L0^2)
    end

    # Handle midpoint for odd P
    if (n + 1) != 2 * nh
        x = 0.0
        L0, _, _ = legendre_poly_rec(n, x)
        xgl[nh + 1] = x
        wgl[nh + 1] = 2 / (n * (n + 1) * L0^2)
    end

    # Use symmetry to fill lower half
    for i in 1:nh
        j = n + 2 - i
        xgl[i] = -xgl[j]
        wgl[i] = wgl[j]
    end

    return xgl, wgl
end

