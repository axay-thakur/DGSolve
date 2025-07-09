# Helper function for run_sim_rk4_dg
export advance_ode_nonlin_rk4, run_sim_rk4_dg
function advance_ode_nonlin_rk4(y0, t0, dt, odefun)
    k1 = dt * odefun(t0, y0)
    k2 = dt * odefun(t0 + dt/2, y0 + k1/2)
    k3 = dt * odefun(t0 + dt/2, y0 + k2/2)
    k4 = dt * odefun(t0 + dt, y0 + k3)
    return y0 + (k1 + 2*k2 + 2*k3 + k4) / 6
end

# --- 10. run_sim_rk4_dg.m ---
function run_sim_rk4_dg(claw::ConservationLaw,Q0::AbstractMatrix, T::AbstractArray, Tvb::AbstractArray, wq::AbstractVector, Mi::AbstractArray, nstep::Int, dt::Real,Gj::AbstractMatrix,Gjbnd::AbstractMatrix, Qlr::Union{AbstractMatrix, Nothing})
    Q_store = zeros(Float64, nstep + 1, size(Q0, 1), size(Q0, 2))
    Q_store[1, :, :] = Q0
    Q = Q0

    for j in 1:nstep
        Q = advance_ode_nonlin_rk4(Q, 0.0, dt, (t_val, Q_val) -> eval_vel_dg(t_val, Q_val,claw,T, Tvb, wq, Mi,Gj,Gjbnd, Qlr))
        Q_store[j + 1, :, :] = Q
    end
    return Q_store
end