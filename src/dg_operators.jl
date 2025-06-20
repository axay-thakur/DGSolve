export eval_vel_dg
# --- 9. eval_vel_dg.m ---
function eval_vel_dg(t::Real, Q::AbstractMatrix,claw::ConservationLaw, T::AbstractArray, Tvb::AbstractArray, w::AbstractVector, Mi::AbstractArray, Qlr::Union{AbstractMatrix, Nothing})
    #nelem = size(Q, 2)
    nv = size(Q, 1)
    nvar = size(T, 2)
    nq = size(T, 4)
    F = zeros(Float64, size(Q))

    Q_append = if isnothing(Qlr) || isempty(Qlr)
        hcat(Q[:, end], Q, Q[:, 1])
    else
        hcat(Qlr[:, 1], Q, Qlr[:, 2])
    end
    nelem_plus_2 = size(Q_append, 2)

    for i in 2:(nelem_plus_2 - 1)
        volume_term = zeros(Float64, nv)
        # vol_flux = zeros(Float64, nvar, nq)

        for j in 1:nq
            T_val_slice_j = T[:, :, 1, j]
            T_deriv_slice_j = T[:, :, 2, j]
            interpolated_state_j = transpose(T_val_slice_j) * Q_append[:, i]
            vol_flux_j = inviscid_flux(claw,interpolated_state_j)
            if claw.viscous
                vol_flux_j -= viscous_flux(claw,nothing,transpose(T_deriv_slice_j) * Q_append[:, i])

            end

            volume_term .+= w[j] .* (T_deriv_slice_j * vol_flux_j)
        end

        Tvb_left_face_val = Tvb[:, :,1, 1]
        Tvb_right_face_val = Tvb[:, :,1, 2]
        
        Fright_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q_append[:, i], transpose(Tvb_right_face_val) * Q_append[:, i+1])
        Fright = Tvb_left_face_val * Fright_numeric

        Fleft_numeric = numerical_flux_inviscid(claw,transpose(Tvb_left_face_val) * Q_append[:, i-1], transpose(Tvb_right_face_val) * Q_append[:, i])
        Fleft = Tvb_right_face_val * Fleft_numeric
        F_full = volume_term .- Fright .+ Fleft
        if claw.viscous
            Tvb_left_face_deriv = Tvb[:,:,2,1]
            Tvb_right_face_deriv = Tvb[:,:,2,2]
            F_full .+= (Tvb_left_face_val-Tvb_right_face_val)* numerical_flux_viscous(claw,nothing,nothing,transpose(Tvb_left_face_deriv)*Q_append[:,i],transpose(Tvb_right_face_deriv) * Q_append[:, i+1])
            F_full .-= (Tvb_right_face_val-Tvb_left_face_val)*numerical_flux_viscous(claw,nothing,nothing,transpose(Tvb_left_face_deriv)*Q_append[:,i-1],transpose(Tvb_right_face_deriv) * Q_append[:, i])
            F_full .+= numerical_flux_viscous(claw,nothing,nothing,Tvb_left_face_deriv,Tvb_right_face_deriv)*(transpose(Tvb_left_face_val)*Q_append[:,i] - transpose(Tvb_right_face_val) * Q_append[:, i+1])
            F_full .-= numerical_flux_viscous(claw,nothing,nothing,Tvb_left_face_deriv,Tvb_right_face_deriv)*(transpose(Tvb_right_face_val)*Q_append[:,i] - transpose(Tvb_left_face_val) * Q_append[:, i-1])
            F_full .-= claw.pen*(Tvb_left_face_val-Tvb_right_face_val)*(transpose(Tvb_left_face_val)*Q_append[:,i] - transpose(Tvb_right_face_val) * Q_append[:, i+1])
            F_full .-= claw.pen*(Tvb_right_face_val-Tvb_left_face_val)*(transpose(Tvb_right_face_val)*Q_append[:,i] - transpose(Tvb_left_face_val) * Q_append[:, i-1])
        end

        F[:, i-1] = Mi[:, :, i-1] * F_full
    end
    return F
end