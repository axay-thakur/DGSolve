export compute_transf_quant, eval_interp_onedim_lagrange, create_elem_basis, Ilrbasis
export eval_modal_onedim_legendre, legendre_poly, legendre_poly_deriv
# --- 1. compute_transf_quant.m ---
function compute_transf_quant(xe::AbstractArray, Q::AbstractArray)
    xq = transpose(xe) * Q[:, 1, :]
    Gj = abs.(transpose(xe) * Q[:, 2, :])
    return xq, Gj
end

# --- 4. eval_interp_onedim_lagrange.m ---
function eval_interp_onedim_lagrange(xk::AbstractVector, x::AbstractVector)
    nv = length(xk)
    nx = length(x)
    xk_col = reshape(xk, :)
    x_row = reshape(x, :)
    Q_values = ones(Float64, nv, nx)
    Q_derivs = zeros(Float64, nv, nx)

    for j in 1:nv
        for i in 1:nv
            if i == j
                continue
            end
            Q_values[j, :] .*= (x_row[:] .- xk_col[i]) ./ (xk_col[j] - xk_col[i])
        end
        for i in 1:nv
            if i == j
                continue
            end
            tmp_prod = ones(Float64, nx)
            for m in 1:nv
                if m == j || m == i
                    continue
                end
                tmp_prod .*= (x_row[:] .- xk_col[m]) ./ (xk_col[j] - xk_col[m])
            end
            Q_derivs[j, :] .+= tmp_prod ./ (xk_col[j] - xk_col[i])
        end
    end
    return reshape(vcat(Q_values, Q_derivs), nv, 2, nx)
end

# --- 2. create_elem_basis.m ---
function create_elem_basis(nvar::Int, Qv::AbstractArray)
    ndim = size(Qv, 2) - 1
    nv, _, nq = size(Qv)
    Tv = zeros(nv * nvar, nvar, ndim + 1, nq)
    for k in 1:nvar
            Tv[k:nvar:end, k, :, :] .= Qv
    end
    return Tv
end

# --- 5. Ilrbasis.m ---
function Ilrbasis(::Val{:nodal},nvar::Int, nv::Int,zk::AbstractVector,ndim::Int=1)
    Tvb = zeros(Float64, nv * nvar, nvar,ndim+1, 2)
    Qder1 = eval_interp_onedim_lagrange(zk,zk[end:end])[:,2,1]
    Qder2 = eval_interp_onedim_lagrange(zk,zk[1:1])[:,2,1]
    Tvb[((nv*nvar - nvar) + 1):end, :,1, 1] .= Matrix{Float64}(I, nvar, nvar)
    Tvb[1:nvar, :,1,2] .= Matrix{Float64}(I, nvar, nvar)
    for k =1:nvar
        Tvb[k:nvar:end,k,2,1] .= Qder1
        Tvb[k:nvar:end,k,2,2] .= Qder2
    end
    return Tvb
end

function Ilrbasis(::Val{:modal},nvar::Int, porder::Int,ndim::Int=1)
    Tvb = zeros(Float64, (porder+1) * nvar, nvar,ndim+1, 2)
    Qbl = eval_modal_onedim_legendre(porder,[1.0])
    Qbr = eval_modal_onedim_legendre(porder, [-1.0])
    for k in 1:nvar
        Tvb[k:nvar:end, k,:, 1] .= Qbl
        Tvb[k:nvar:end, k,:, 2] .= Qbr
    end
    return Tvb
end

"""
    legendre_poly(n::Int, x::Real)

Evaluates the standard (non-normalized) Legendre polynomial of degree `n` at point `x`.
Uses the recursive definition.
"""
function legendre_poly(n::Int, x::Real)
    if n == 0
        return 1.0
    elseif n == 1
        return x
    else
        # Legendre recursion formula: P_n(x) = ( (2n-1)*x*P_{n-1}(x) - (n-1)*P_{n-2}(x) ) / n
        return ((2n - 1) * x * legendre_poly(n - 1, x) - (n - 1) * legendre_poly(n - 2, x)) / n
    end
end

"""
    legendre_poly_deriv(n::Int, x::Real)

Evaluates the derivative of the standard (non-normalized) Legendre polynomial of degree `n` at point `x`.
Uses a stable recursive relation.
"""
function legendre_poly_deriv(n::Int, x::Real)
    if n == 0
        return 0.0 # P_0(x) = 1, derivative is 0
    elseif n == 1
        return 1.0 # P_1(x) = x, derivative is 1
    else
        if abs(x) ≈ 1.0 # Handle endpoints specially to avoid division by zero
            if x ≈ 1.0
                return n * (n + 1) / 2.0
            else # x ≈ -1.0
                return n * (n + 1) / 2.0 * (-1.0)^(n + 1)
            end
        else
            return n / (1 - x^2) * (legendre_poly(n - 1, x) - x * legendre_poly(n, x))
        end
    end
end


"""
    eval_modal_onedim_legendre(porder::Int, x::AbstractVector)

Evaluates **orthonormal** Legendre polynomials up to `porder` and their derivatives
at points `x` on the reference interval [-1, 1].

The orthonormalization factor for Legendre polynomial P_n(x) is `sqrt((2n+1)/2)`.

# Input arguments
- `porder`: The maximum degree of Legendre polynomial (degree 0 to `porder`).
- `x`: Array of points where Legendre polynomials are evaluated (reference interval [-1, 1]).

# Output arguments
- `Q`: Array (PORDER+1, 2, NX).
    - `Q[k, 1, j]` = value of orthonormal Legendre polynomial L_{k-1}(x_j)
    - `Q[k, 2, j]` = derivative of orthonormal Legendre polynomial L'_{k-1}(x_j)
"""
function eval_modal_onedim_legendre(porder::Int, x::AbstractVector)
    nx = length(x)
    Q = zeros(Float64,porder + 1, 2, nx)

    for i in 1:(porder + 1) # Loop for polynomial degree k = 0 to porder
        n = i - 1 # Current polynomial degree (0-indexed, corresponding to P_n)
        
        # Calculate the orthonormalization factor
        # L_n(x) = sqrt((2n+1)/2) * P_n(x)
        # ortho_factor = sqrt((2n + 1) / 2.0)
        ortho_factor = Float64(1.0)
        for j in 1:nx # Loop for evaluation points
            # Apply normalization to both the value and the derivative
            Q[i, 1, j] = ortho_factor * legendre_poly(n, x[j])
            Q[i, 2, j] = ortho_factor * legendre_poly_deriv(n, x[j])
        end
    end
    return Q
end
