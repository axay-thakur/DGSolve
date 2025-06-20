export ConservationLaw, ScalarAdvectionDiffusion
export inviscid_flux, viscous_flux
export numerical_flux_inviscid, numerical_flux_viscous

abstract type ConservationLaw end

@inline inviscid_flux(::ConservationLaw,u) = throw("Not implemented")
@inline viscous_flux(::ConservationLaw, u, ∇u) = throw("Not implemented")
@inline  numerical_flux_viscous(::ConservationLaw, ul, ur, ∇ul, ∇ur) = throw("Not impelemented")
@inline  numerical_flux_inviscid(::ConservationLaw,ul,ur) =  throw("Not implemented")
@inline source_term(::ConservationLaw,u) = throw("Not implemented")

struct ScalarAdvectionDiffusion{T,S} <: ConservationLaw
    a::T         # advection speed
    ν::T         # viscosity
    dim::S
    pen::T
    viscous::Bool
end

function ScalarAdvectionDiffusion(a::T,dim::S) where{T,S}
    return ScalarAdvectionDiffusion{T, S}(a, 0, 0,dim, false)
end

function ScalarAdvectionDiffusion(a::T,ν::T,dim::S,Cip::T,porder::S,h::T) where {T,S}
    pen = Cip*ν*(porder^2)/h
    return ScalarAdvectionDiffusion{T, S}(a, ν, dim,pen, true) 
end
function inviscid_flux(claw::ScalarAdvectionDiffusion,u)
    return (claw.a * u)
end

function viscous_flux(claw::ScalarAdvectionDiffusion, u, ∇u)
        return claw.ν * ∇u
end

function numerical_flux_viscous(claw::ScalarAdvectionDiffusion, ul, ur, ∇ul, ∇ur)
        return claw.ν * 0.5 * (∇ul + ∇ur)
end

function numerical_flux_inviscid(claw::ScalarAdvectionDiffusion,ul,ur)
    if claw.a < 0
        fs = claw.a * ur
    else
        fs = claw.a * ul
    end
    return fs
end
