export solve_linear_advection_dg_nodal, solve_linear_advection_dg_modal
function solve_linear_advection_dg_nodal(;c=1.0,a = -1.0,b = 1.0,nelem = 50,porder = 2,nvar = 1,dt = 0.002)
    gr() # Set the plotting backend for Plots.jl

    

    msh = mesh_utility(a, b, nelem, porder)
    xe = msh.e2vcg
    Q0 = exp.(-(xe.^2) ./ 0.1)

    w, z = define_quad_onedim_gaussleg(porder + 4)

    zk = collect(range(-1.0, stop=1.0, length=porder + 1))
    Q_basis = eval_interp_onedim_lagrange(zk, z)

    xq, Gj = compute_transf_quant(xe, Q_basis)

    Tv = create_elem_basis(nvar, Q_basis)
    Tvb = Ilrbasis(Val(:nodal),nvar, porder + 1,zk)

    Minv = inv_mass_matrix(w, Gj, Tv)

    claw = ScalarAdvectionDiffusion(c,1)

    T_final = 2.0
    nstep1 = ceil(Int, T_final / dt)

    Q11 = run_sim_rk4_dg(claw, Q0, Tv, Tvb, w, Minv, nstep1, dt, nothing)

    # Plotting
    # Create the initial plot object
    plt = plot(title="Linear Advection - DG P-order = $(porder)",
               xlabel="x", ylabel="Q(x)", legend=:topright)

    # Plot initial condition (t=0)
    plot_dg1d_soln(Val(:nodal), plt, xe, Q11[1, :, :], ["r-"], label="time = 0 s")

    # Plot final condition (t=T_final)
    plot_dg1d_soln(Val(:nodal), plt, xe, Q11[end, :, :], ["k-"], label="time = $(T_final) s")

    # Save and display
    # savefig(plt, "linaddgpdim2.png")
    display(plt)
    return Q11, Q_basis
end

function solve_linear_advection_dg_modal(;c=1.0,a = -1.0,b = 1.0,nelem = 50,porder = 2,nvar = 1,dt = 0.002)
    gr() # Set the plotting backend for Plots.jl

    msh = mesh_utility(a, b, nelem, porder)
    xe = msh.e2vcg
    
    zk = collect(range(-1.0, stop=1.0, length=porder + 1))
    w, z = define_quad_onedim_gaussleg(porder + 4)

    # zk = collect(range(-1.0, stop=1.0, length=porder + 1))
    Q_basis = eval_modal_onedim_legendre(porder, z)
    Q_basis_lagrange = eval_interp_onedim_lagrange(zk, z)
    xq, Gj = compute_transf_quant(xe, Q_basis_lagrange)
    Q0 = exp.(-(transpose(xq).^2) ./ 0.1)
    Q0_modal = pinv(transpose(Q_basis[:,1,:]))*Q0

    Tv = create_elem_basis(nvar, Q_basis)
    Tvb = Ilrbasis(Val(:modal),nvar, porder)

    Minv = inv_mass_matrix(w, Gj, Tv)

    claw = ScalarAdvectionDiffusion(c,1)

    T_final = 2.0
    nstep1 = ceil(Int, T_final / dt)

    Q11 = run_sim_rk4_dg(claw,Q0_modal, Tv, Tvb, w, Minv, nstep1, dt, nothing)

    # Plotting
    # Create the initial plot object
    plt = plot(title="Linear Advection - DG Modal P-order = $(porder)",
               xlabel="x", ylabel="Q(x)", legend=:topright)

    # Plot initial condition (t=0)
    plot_dg1d_soln(Val(:modal), plt, xe, Q11[1, :, :], ["r-"], label="time = 0 s")

    # Plot final condition (t=T_final)
    plot_dg1d_soln(Val(:modal), plt, xe, Q11[2, :, :], ["k-"], label="time = $(T_final) s")

    # Save and display
    # savefig(plt, "linaddgmpdim2.png")
    display(plt)
    return Q11, Q_basis
end