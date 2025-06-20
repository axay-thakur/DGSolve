export plot_dg1d_soln
# --- 12. plot_dg1d_soln.m ---
function plot_dg1d_soln(::Val{B},plot_obj, xdg::AbstractArray, Q::AbstractArray,pltsty::Vector{String}=nothing; label::String="") where {B}
    if ndims(Q) == 3
        Q_sol = Q
    elseif ndims(Q) == 2
        Q_sol = Q
    else
        error("Unsupported dimensions for Q. Expected 2D ((PORDER+1), NELEM) or 3D (time_idx, (PORDER+1), NELEM).")
    end

    

    line_color = :black
    line_style = :solid
    if !isnothing(pltsty) && !isempty(pltsty)
        style_str = pltsty[1]
        if occursin("k", style_str)
            line_color = :black
        elseif occursin("r", style_str)
            line_color = :red
        elseif occursin("b", style_str)
            line_color = :blue
        elseif occursin("g", style_str)
            line_color = :green
        end
        if occursin("--", style_str)
            line_style = :dash
        elseif occursin(":", style_str)
            line_style = :dot
        end
    end

    if isnothing(plot_obj)
        plot_obj = plot()
    end

    nelem = size(xdg, ndims(xdg))
    if B == :nodal
        nvpe = size(xdg, ndims(xdg) - 1)
        
        porder = nvpe - 1
        N_plot_points_per_elem = (porder > 1) ? 20 : 2
        # To avoid repeated labels for each segment, we add the label only once
        first_segment_plotted = false
        for e in 1:nelem
            local_xdg_nodes = (ndims(xdg) == 3) ? vec(xdg[1, :, e]) : vec(xdg[:, e])
            local_Q_nodes = vec(Q_sol[:, e])
            
            p_poly_obj = fit(local_xdg_nodes, local_Q_nodes, porder)

            x_plot_elem = collect(range(local_xdg_nodes[1], stop=local_xdg_nodes[end], length=N_plot_points_per_elem))
            y_plot_elem = p_poly_obj.(x_plot_elem)

            if !first_segment_plotted
                plot!(plot_obj, x_plot_elem, y_plot_elem, line=(line_color, line_style, 2), label=label)
                first_segment_plotted = true
            else
                plot!(plot_obj, x_plot_elem, y_plot_elem, line=(line_color, line_style, 2), label="")
            end
        end
    elseif B == :modal
        # For modal representation, we plot each polynomial degree separately
        porder = size(Q_sol, 1) - 1
        N_plot_points_per_elem = (porder > 1) ? 20 : 2
        first_segment_plotted = false
        for e in 1:nelem
            local_Q_nodes = vec(Q_sol[:, e])
            z_pts = collect(range(-1.0, stop=1.0, length=N_plot_points_per_elem))
            zk = collect(range(-1.0, stop=1.0, length=porder + 1))
            Q = eval_modal_onedim_legendre(porder, z_pts)
            Q_gmap = eval_interp_onedim_lagrange(zk, z_pts)
            x_plot_elem = vec(compute_transf_quant(xdg[:, e:e], Q_gmap)[1])
            y_plot_elem = vec(local_Q_nodes'*Q[:,1,:])
            if !first_segment_plotted
                plot!(plot_obj, x_plot_elem, y_plot_elem, line=(line_color, line_style, 2), label=label)
                first_segment_plotted = true
            else
                plot!(plot_obj, x_plot_elem, y_plot_elem, line=(line_color, line_style, 2), label="")
            end
        end
    else
        error("Unsupported basis type: $B")          
    end
    return plot_obj
end