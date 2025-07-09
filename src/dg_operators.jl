export eval_vel_dg
# --- 9. eval_vel_dg.m ---
function eval_vel_dg(t::Real, Q::AbstractMatrix,claw::ConservationLaw, T::AbstractArray, Tvb::AbstractArray, w::AbstractVector, Mi::AbstractArray,Gj::AbstractMatrix, Gjbnd::AbstractMatrix,Qlr::Union{AbstractMatrix, Nothing})
    nv = size(Q, 1)
    nvar = size(T, 2)
    nelem = size(Q,2)
    nq = size(T, 4)
    F = zeros(Float64, size(Q))

    Qlr = if isnothing(Qlr) || isempty(Qlr)
        hcat(Q[end:end, end],Q[1:1, 1])
    else
        Qlr
    end
    Tvb_left_face_val = Tvb[:, :,1, 1]
    Tvb_right_face_val = Tvb[:, :,1, 2]
    Tvb_left_face_deriv = Tvb[:,:,2,1]
    Tvb_right_face_deriv = Tvb[:,:,2,2]
    if claw.viscous
        if isempty(Qlr[nvar+1:end,:])
            catted = hcat((1/Gjbnd[1,2]).*transpose(Tvb_right_face_deriv) * Q[:, 1:1],(1/Gjbnd[end,1]).*transpose(Tvb_left_face_deriv)*Q[:,end:end])
            Qlr = vcat(Qlr,catted)
        end
    end
    for i in 1:nelem
        volume_term = zeros(Float64, nv)
        Fright = zeros(Float64, nv)
        Fleft = zeros(Float64,nv)
        

        for j in 1:nq
            T_val_slice_j = T[:, :, 1, j]
            T_deriv_slice_j = T[:, :, 2, j]
            interpolated_state_j = transpose(T_val_slice_j) * Q[:, i]
            vol_flux = zeros(Float64,nvar)
            if claw.inviscid
                vol_flux .+= inviscid_flux(claw,interpolated_state_j)
            end
            if claw.viscous
                vol_flux .-= viscous_flux(claw,nothing,(1/Gj[i,j]).*transpose(T_deriv_slice_j) * Q[:, i])

            end
            volume_term .+= w[j] .* (T_deriv_slice_j * vol_flux)
        end
        F_full = volume_term 
        if claw.inviscid
            if i == nelem
                Fright_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q[:, i],Qlr[1:nvar, 2])
                Fright .= Tvb_left_face_val * Fright_numeric
            else
                Fright_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q[:, i], transpose(Tvb_right_face_val) * Q[:, i+1])
                Fright .= Tvb_left_face_val * Fright_numeric
            end
            if i ==1
                Fleft_numeric = numerical_flux_inviscid(claw,Qlr[1:nvar, 1], transpose(Tvb_right_face_val) * Q[:, i])
                Fleft .= Tvb_right_face_val * Fleft_numeric
            else
                Fleft_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q[:, i-1], transpose(Tvb_right_face_val) * Q[:, i])
                Fleft .= Tvb_right_face_val * Fleft_numeric
            end
            F_full .+= Fleft .- Fright
        end
        

        if claw.viscous
            
            if i == nelem
                F_full .+= Tvb_left_face_val*numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*transpose(Tvb_left_face_deriv)*Q[:,i],Qlr[nvar+1:end,2])
                F_full[end:end] .+= claw.ν*(1/Gjbnd[i,1])*Tvb_left_face_deriv[end]*((transpose(Tvb_left_face_val)*Q[:,i]) - ((transpose(Tvb_left_face_val)*Q[:,i]) + Qlr[1:nvar,2])*0.5)
                F_full .+= claw.pen*Tvb_left_face_val*(transpose(Tvb_left_face_val)*Q[:,i] - Qlr[1:nvar,2])
            else
                F_full .+= Tvb_left_face_val*numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*transpose(Tvb_left_face_deriv)*Q[:,i],(1/Gjbnd[i+1,2]).*transpose(Tvb_right_face_deriv) * Q[:, i+1])
                F_full[end:end] .+= claw.ν*(1/Gjbnd[i,1])*Tvb_left_face_deriv[end]*((transpose(Tvb_left_face_val)*Q[:,i]) - ((transpose(Tvb_left_face_val)*Q[:,i]) + (transpose(Tvb_right_face_val) * Q[:, i+1]))*0.5)
                F_full .+= claw.pen*Tvb_left_face_val*(transpose(Tvb_left_face_val)*Q[:,i] - transpose(Tvb_right_face_val) * Q[:, i+1])
            end
            if i == 1
                F_full .-= Tvb_right_face_val*numerical_flux_viscous(claw,nothing,nothing,Qlr[nvar+1:end,1],(1/Gjbnd[i,2]).*transpose(Tvb_right_face_deriv) * Q[:, i])
                F_full .-= claw.pen*Tvb_right_face_val*(transpose(Tvb_right_face_val)*Q[:,i] - Qlr[1:nvar,1])     # F_full .+= numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*Tvb_left_face_deriv,(1/Gjbnd[i+1,2]).*Tvb_right_face_deriv)*(transpose(Tvb_left_face_val)*Q_append[:,i] - transpose(Tvb_right_face_val) * Q_append[:, i+1])
                F_full[1:1] .-= claw.ν*(1/Gjbnd[i,2])*Tvb_right_face_deriv[1]*((transpose(Tvb_right_face_val)*Q[:,i]) - ((transpose(Tvb_right_face_val)*Q[:,i]) + Qlr[1:nvar,1])*0.5) 
            else
                F_full .-= Tvb_right_face_val*numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i-1,1]).*transpose(Tvb_left_face_deriv)*Q[:,i-1],(1/Gjbnd[i,2]).*transpose(Tvb_right_face_deriv) * Q[:, i])
                F_full .-= claw.pen*Tvb_right_face_val*(transpose(Tvb_right_face_val)*Q[:,i] - transpose(Tvb_left_face_val) * Q[:, i-1])     # F_full .+= numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*Tvb_left_face_deriv,(1/Gjbnd[i+1,2]).*Tvb_right_face_deriv)*(transpose(Tvb_left_face_val)*Q_append[:,i] - transpose(Tvb_right_face_val) * Q_append[:, i+1])
                F_full[1:1] .-= claw.ν*(1/Gjbnd[i,2])*Tvb_right_face_deriv[1]*((transpose(Tvb_right_face_val)*Q[:,i]) - ((transpose(Tvb_right_face_val)*Q[:,i]) + (transpose(Tvb_left_face_val) * Q[:, i-1]))*0.5) 
            end 
        end

        F[:, i] = Mi[:, :, i] * F_full
        if claw.src !== nothing
            F[:, i] .+= claw.src[:,i]
        end
    end
    return F
end

function eval_vel_dg(Q::AbstractMatrix,claw::ConservationLaw, T::AbstractArray, Tvb::AbstractArray, w::AbstractVector,Gj::AbstractMatrix, Gjbnd::AbstractMatrix,Qlr::Union{AbstractMatrix, Nothing})

    nv = size(Q, 1)
    nvar = size(T, 2)
    nelem = size(Q,2)
    nq = size(T, 4)
    F = zeros(Float64, size(Q))

    Qlr = if isnothing(Qlr) || isempty(Qlr)
        hcat(Q[end:end, end],Q[1:1, 1])
    else
        Qlr
    end
    Tvb_left_face_val = Tvb[:, :,1, 1]
    Tvb_right_face_val = Tvb[:, :,1, 2]
    Tvb_left_face_deriv = Tvb[:,:,2,1]
    Tvb_right_face_deriv = Tvb[:,:,2,2]
    if claw.viscous
        if isempty(Qlr[nvar+1:end,:])
            catted = hcat((1/Gjbnd[1,2]).*transpose(Tvb_right_face_deriv) * Q[:, 1:1],(1/Gjbnd[end,1]).*transpose(Tvb_left_face_deriv)*Q[:,end:end])
            Qlr = vcat(Qlr,catted)
        end
    end
    for i in 1:nelem
        volume_term = zeros(Float64, nv)
        Fright = zeros(Float64, nv)
        Fleft = zeros(Float64,nv)
        

        # vol_flux = zeros(Float64, nvar, nq)

        for j in 1:nq
            T_val_slice_j = T[:, :, 1, j]
            T_deriv_slice_j = T[:, :, 2, j]
            interpolated_state_j = transpose(T_val_slice_j) * Q[:, i]
            vol_flux = zeros(Float64,nvar)
            if claw.inviscid
                vol_flux .-= inviscid_flux(claw,interpolated_state_j)
            end
            if claw.viscous
                vol_flux .+= viscous_flux(claw,nothing,(1/Gj[i,j]).*transpose(T_deriv_slice_j) * Q[:, i])

            end
            volume_term .+= w[j] .* (T_deriv_slice_j * vol_flux)
        end
        F_full = volume_term 
        if claw.inviscid
            if i == nelem
                Fright_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q[:, i],Qlr[1:nvar, 2])
                Fright .= Tvb_left_face_val * Fright_numeric
            else
                Fright_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q[:, i], transpose(Tvb_right_face_val) * Q[:, i+1])
                Fright .= Tvb_left_face_val * Fright_numeric
            end
            if i ==1
                Fleft_numeric = numerical_flux_inviscid(claw,Qlr[1:nvar, 1], transpose(Tvb_right_face_val) * Q[:, i])
                Fleft .= Tvb_right_face_val * Fleft_numeric
            else
                Fleft_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q[:, i-1], transpose(Tvb_right_face_val) * Q[:, i])
                Fleft .= Tvb_right_face_val * Fleft_numeric
            end
            F_full .+= Fright .- Fleft
        end
        

        if claw.viscous
            
            if i == nelem
                F_full .-= Tvb_left_face_val*numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*transpose(Tvb_left_face_deriv)*Q[:,i],Qlr[nvar+1:end,2])
                F_full[end:end] .-= claw.ν*(1/Gjbnd[i,1])*Tvb_left_face_deriv[end]*((transpose(Tvb_left_face_val)*Q[:,i]) - ((transpose(Tvb_left_face_val)*Q[:,i]) + Qlr[1:nvar,2])*0.5)
                F_full .-= claw.pen*Tvb_left_face_val*(transpose(Tvb_left_face_val)*Q[:,i] - Qlr[1:nvar,2])
            else
                F_full .-= Tvb_left_face_val*numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*transpose(Tvb_left_face_deriv)*Q[:,i],(1/Gjbnd[i+1,2]).*transpose(Tvb_right_face_deriv) * Q[:, i+1])
                F_full[end:end] .-= claw.ν*(1/Gjbnd[i,1])*Tvb_left_face_deriv[end]*((transpose(Tvb_left_face_val)*Q[:,i]) - ((transpose(Tvb_left_face_val)*Q[:,i]) + (transpose(Tvb_right_face_val) * Q[:, i+1]))*0.5)
                F_full .-= claw.pen*Tvb_left_face_val*(transpose(Tvb_left_face_val)*Q[:,i] - transpose(Tvb_right_face_val) * Q[:, i+1])
            end
            if i == 1
                F_full .+= Tvb_right_face_val*numerical_flux_viscous(claw,nothing,nothing,Qlr[nvar+1:end,1],(1/Gjbnd[i,2]).*transpose(Tvb_right_face_deriv) * Q[:, i])
                F_full .+= claw.pen*Tvb_right_face_val*(transpose(Tvb_right_face_val)*Q[:,i] - Qlr[1:nvar,1])     # F_full .+= numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*Tvb_left_face_deriv,(1/Gjbnd[i+1,2]).*Tvb_right_face_deriv)*(transpose(Tvb_left_face_val)*Q_append[:,i] - transpose(Tvb_right_face_val) * Q_append[:, i+1])
                F_full[1:1] .+= claw.ν*(1/Gjbnd[i,2])*Tvb_right_face_deriv[1]*((transpose(Tvb_right_face_val)*Q[:,i]) - ((transpose(Tvb_right_face_val)*Q[:,i]) + Qlr[1:nvar,1])*0.5) 
            else
                F_full .+= Tvb_right_face_val*numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i-1,1]).*transpose(Tvb_left_face_deriv)*Q[:,i-1],(1/Gjbnd[i,2]).*transpose(Tvb_right_face_deriv) * Q[:, i])
                F_full .+= claw.pen*Tvb_right_face_val*(transpose(Tvb_right_face_val)*Q[:,i] - transpose(Tvb_left_face_val) * Q[:, i-1])     # F_full .+= numerical_flux_viscous(claw,nothing,nothing,(1/Gjbnd[i,1]).*Tvb_left_face_deriv,(1/Gjbnd[i+1,2]).*Tvb_right_face_deriv)*(transpose(Tvb_left_face_val)*Q_append[:,i] - transpose(Tvb_right_face_val) * Q_append[:, i+1])
                F_full[1:1] .+= claw.ν*(1/Gjbnd[i,2])*Tvb_right_face_deriv[1]*((transpose(Tvb_right_face_val)*Q[:,i]) - ((transpose(Tvb_right_face_val)*Q[:,i]) + (transpose(Tvb_left_face_val) * Q[:, i-1]))*0.5) 
            end 
        end

        F[:, i] = F_full
    end
    return F
end
