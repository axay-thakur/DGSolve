export solve_linear_advection_dg_nodal, solve_linear_diffusion_dg_nodal_steady
function solve_linear_advection_dg_nodal(;c = 1.0,a = -1.0,b = 1.0,nelem = 50,porder = 2, nvar = 1, dt = 0.002)
    gr() # Set the plotting backend for Plots.jl

    msh = mesh_utility(a, b, nelem, porder)
    xe = msh.e2vcg
    Q0 = exp.(-(xe.^2) ./ 0.1)

    w, z = define_quad_onedim_gaussleg(porder + 4)

    zk = collect(range(-1.0, stop=1.0, length=porder + 1))
    Q_basis = eval_interp_onedim_lagrange(zk, z)
    Qbnd = eval_interp_onedim_lagrange(zk,zk[[end,1]])


    xq, Gj = compute_transf_quant(xe, Q_basis)
    xqnd, Gjbnd = compute_transf_quant(xe, Qbnd)
    Gjbnd = vcat(Gjbnd[1:1,:],Gjbnd, Gjbnd[end:end,:])

    Tv = create_elem_basis(nvar, Q_basis)
    Tvb = Ilrbasis(Val(:nodal),nvar, porder + 1,zk)

    
    Minv = inv_mass_matrix(w, Gj, Tv)

    claw = ScalarAdvectionDiffusion(c,1)

    T_final = 2.0
    nstep1 = ceil(Int, T_final / dt)

    Q = run_sim_rk4_dg(claw, Q0, Tv, Tvb, w, Minv, nstep1, dt,Gj,Gjbnd, nothing)

    # Plotting
    # Create the initial plot object
    plt = plot(title="Linear Advection - DG P-order = $(porder)",
               xlabel="x", ylabel="Q(x)", legend=:topright)

    # Plot initial condition (t=0)
    plot_dg1d_soln(Val(:nodal), plt, xe, Q[1, :, :], ["r-"], label="time = 0 s")

    # Plot final condition (t=T_final)
    plot_dg1d_soln(Val(:nodal), plt, xe, Q[end, :, :], ["k-"], label="time = $(T_final) s")

    # Save and display
    # savefig(plt, "linaddgpdim2.png")
    display(plt)
    return Q
end

function solve_linear_diffusion_dg_nodal_steady(;ν=1.0,a=-1.0,b=1.0,nelem=4,porder=16,nvar=1,Cip=1e4)
    zk,_ = legendre_gauss_lobatto(porder+1)
    msh = mesh_utility(a,b,nelem,porder,zk)
    xe = msh.e2vcg
    Q0 = sin.(2*pi*xe)
    Qsol = sin.(2*pi*xe)
    dQsol = 2*pi*cos.(2*pi*xe)
    src = (2*pi)^2*sin.(2*pi*xe)
    ndof = nelem*(porder+1) 

    z, w = legendre_gauss_lobatto(porder + 1)
    Gj = ones(Float64,nelem,porder+1)*(b-a)/nelem/2
    Gjbnd = ones(Float64,nelem,2)*(b-a)/nelem/2

    Q_basis = eval_interp_onedim_lagrange(zk, z)
    Tv = create_elem_basis(nvar, Q_basis)
    Tvb = Ilrbasis(Val(:nodal),nvar, porder + 1,zk)

    claw = ScalarAdvectionDiffusion(ν,1,Cip,porder,(b-a)/nelem)

    Qlr = vcat(hcat(Qsol[1:1,1],Qsol[end:end,end]),hcat(dQsol[1:1,1],dQsol[end:end,end]))
    Lmatrix = zeros(Float64,ndof,ndof)
    rhs = zeros(Float64,size(Q0))
    for i in 1:nelem
        for j in 1:(porder+1)
        Q = zeros(Float64, porder+1, nelem)
        Q[j,i] = 1.0
        Lmatrix[:,(i-1)*(porder+1)+j] = reshape(eval_vel_dg(Q,claw,Tv, Tvb, w,Gj,Gjbnd, Qlr),:)
        end
        rhs[:,i] .= element_mass_matrix(w, Gj[i,:], Tv)*src[:,i]
    end
    rhs = reshape(rhs,:)
    Q = reshape(Lmatrix\rhs,:,nelem);
    error = norm(Q - Qsol, 2)/norm(Qsol,2)
    println("L2 error in steady diffusion solution: ", error)
    plt = plot(title="Diffusion Steady - DG Nodal P-order = $(porder)",
                 xlabel="x", ylabel="Q(x)", legend=:topright)
    plot_dg1d_soln(Val(:nodal), plt, xe, Q[:, :], ["r-"], label="Numerical")
    plot_dg1d_soln(Val(:nodal), plt, xe, Qsol, ["k-"], label="Analytical")
    display(plt)
    return Q
end
