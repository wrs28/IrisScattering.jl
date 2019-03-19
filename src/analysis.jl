#TODO: analyze_into_waveguides account for phase, for now it's not accounted for at all
#TODO: waveguide analysis only for perpendicular waveguides
################################################################################
### ANALYZE FIELD
################################################################################
#@code_warntype checked
"""
    c = analyze_output(sim, k, ψ, m)

c is the output coefficient in the mth channel

S is constructed from cs for unit input on each channel
"""
function analyze_field(sim::Simulation, k, ψ, m; direction=:out)
    if isempty(sim.sys.waveguides)
        c = analyze_into_angular_momentum(sim, k, ψ, m, direction)
    elseif any(ispc.(sim.sys.domains[get_waveguide_domains(sim,sim.sct.channels[m].waveguide)]))
        c = analyze_into_pc_waveguide(sim, k, ψ, m, direction)
    elseif (sim.sys.domains[ScalarFDFD.get_waveguide_   domains(sim, sim.sct.channels[m].waveguide)][end].domain_type) == :halfspace_waveguide
        c = analyze_into_halfspace_waveguide(sim, k, ψ, m, direction)
    elseif (sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, sim.sct.channels[m].waveguide)][end].domain_type) == :planar_waveguide
        c = analyze_into_planar_waveguide(sim, k, ψ, m, direction)
    else
        throw(ArgumentError("unrecognized waveguide type"))
    end

    return c
end

#
"""
    c = analyze_into_waveguides_pc(sim, k, ψ, m)
"""
function analyze_into_pc_waveguide(sim::Simulation, k, ψ, m, direction)

    xs = sim.dis.x[1][1] .+ (0:sim.dis.N[1]-1)*sim.dis.dx[1]
    ys = sim.dis.x[2][1] .+ (0:sim.dis.N[2]-1)*sim.dis.dx[2]
    ψ_itp = CubicSplineInterpolation((xs,ys),reshape(ψ[:,1],sim.dis.N[1],sim.dis.N[2]), bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))

    u, _ = incident_pc_mode(sim, k, m, direction)
    u_itp = CubicSplineInterpolation((xs,ys),reshape(u[:,1],sim.dis.N[1],sim.dis.N[2]), bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))

    wg_sim = ScalarFDFD.extract_waveguide_simulation(sim, sim.sct.channels[m].waveguide)

    if ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) == :left
        wg_sim.bnd.∂Ω[1,1] = sim.bnd.∂Ω_tr[1,1]
        wg_sim.bnd.∂Ω[2,1] = sim.bnd.∂Ω_tr[1,1] + wg_sim.lat.a
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) == :right
        wg_sim.bnd.∂Ω_tr[2,1] = sim.bnd.∂Ω_tr[2,1]
        wg_sim.bnd.∂Ω_tr[1,1] = sim.bnd.∂Ω_tr[2,1] - wg_sim.lat.a
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) == :bottom
        wg_sim.bnd.∂Ω_tr[1,2] = sim.bnd.∂Ω_tr[1,2]
        wg_sim.bnd.∂Ω_tr[2,2] = sim.bnd.∂Ω_tr[1,2] + wg_sim.lat.b
    elseif ScalarFDFD.get_asymptote(sim, sim.sct.channels[m].waveguide) == :top
        wg_sim.bnd.∂Ω_tr[2,2] = sim.bnd.∂Ω_tr[2,2]
        wg_sim.bnd.∂Ω_tr[1,2] = sim.bnd.∂Ω_tr[2,2] - wg_sim.lat.b
    else
        throw(ArgumentError("invalid asymptote $(get_asymptote(sim,sim.sct.channels[i].waveguide)) for channel $i"))
    end

    wg_sim2 = Simulation(wg_sim)
    xs = wg_sim2.dis.x[1] .+ wg_sim2.dis.dx[1]
    ys = wg_sim2.dis.x[2] .+ wg_sim2.dis.dx[2]

    c = quadrature(wg_sim2, (ψ_itp.(xs,ys).*conj(u_itp.(xs,ys)))[:]; weight=:ε)/quadrature(wg_sim2, abs2.(u_itp.(xs,ys))[:]; weight=:ε)

    return c::ComplexF64
end


"""
    c = analyze_into_waveguides_planar(sim, k, ψ, m, direction; η_init=-3)
"""
function analyze_into_planar_waveguide(sim::Simulation, k, ψ, m, direction; η_init=-3)

    side = sim.bnd.waveguides[sim.sct.channels[m].waveguide].side
    xs = sim.dis.x[1][1]:sim.dis.dx[1]:sim.dis.x[1][end]
    ys = sim.dis.x[2][1]:sim.dis.dx[2]:sim.dis.x[2][end]
    Ψ = CubicSplineInterpolation((xs,ys), reshape(ψ,sim.dis.N[1],sim.dis.N[2]))

    β, u_itp, collapsed_dim, full_dim = transverse_mode(sim, k, m)
    φ = u_itp(sim.dis.x[full_dim])
    if side == :l
        P = Ψ(sim.dis.x[1][1],sim.dis.x[2])
    elseif side == :r
        P = Ψ(sim.dis.x[1][end],sim.dis.x[2])
    elseif side == :b
        P = Ψ(sim.dis.x[1],sim.dis.x[2][1])
    elseif side == :t
        P = Ψ(sim.dis.x[1],sim.dis.x[2][end])
    else
        throw(ErrorException("Invalid side $(side)"))
    end

    c = sqrt(real(β))*quadrature((conj(φ).*P)[:], sim.dis.dx[full_dim])
    return c
end


#@code_warntype checked
"""
    c = analyze_into_waveguides_halfspace(sim,k,m)
"""
function analyze_into_halfspace_waveguide(sim, k, ψ, m)

    _,_,_,β = incident_mode_halfspace(sim, k, m)

    if get_asymptote(sim, sim.sct.channels[m].waveguide) == :left
        ind = sim.dis.X_idx[1]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) == :right
        ind = sim.dis.X_idx[end]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) == :bottom
        ind = sim.dis.X_idx[1]
    elseif get_asymptote(sim, sim.sct.channels[m].waveguide) == :top
        ind = sim.dis.X_idx[end]
    else
        throw(ArgumentError("invalid asymptote $(get_asymptote(sim,sim.sct.channels[i].waveguide)) for channel $i"))
    end
    c = ψ[ind]*β
end


"""
    c = analyze_into_angular_momentum(sim, k, ψ, m, direction)
"""
function analyze_into_angular_momentum(sim::Simulation, k, ψ, m, direction)

    nθ = ANALYSIS_QUADRATURE_POINTS
    θ = LinRange(0,2π,nθ)
    dθ = θ[2]-θ[1]

    # R is radius at which to interpolate
    R = (findmin(abs.(sim.bnd.∂Ω_tr)))[1] - findmax(sim.dis.dx)[1]
    X = R*cos.(θ)
    Y = R*sin.(θ)

    # interpolate wavefunction at r=R, result is P(θ)
    xs = sim.dis.x[1][1]:sim.dis.dx[1]:sim.dis.x[1][end]
    ys = sim.dis.x[2][1]:sim.dis.dx[2]:sim.dis.x[2][end]
    Ψ = CubicSplineInterpolation((xs,ys), reshape(ψ,sim.dis.N[1],sim.dis.N[2]))
    P = Ψ.(X,Y)

    q = sim.sct.channels[m].quantum_number

    if direction == :in
        c = quadrature(exp.(-1im*q*θ).*P, dθ)/(π*hankelh2(q,k*R))
    elseif direction == :out
        c = quadrature(exp.(-1im*q*θ).*P, dθ)/(π*hankelh1(q,k*R))
    else
        throw(ErrorException("Invalid direction $(direction). Must be one of :in or :out"))
    end

    return c
end


################################################################################
### ABSORPTION AND FLUX
################################################################################
#@code_warntype checked
"""
    flux, (left, right, bottom, top) = surface_flux(sim, ψ)
"""
function surface_flux(sim::Simulation,Ψ)

    flux = Array{Float64}(undef,size(Ψ,2))
    top = Array{Float64}(undef,size(flux))
    bottom = Array{Float64}(undef,size(flux))
    left = Array{Float64}(undef,size(flux))
    right = Array{Float64}(undef,size(flux))

    for i in 1:size(Ψ,2)
        ψ = reshape(Ψ[sim.dis.X_idx,i],sim.dis.N_tr[1],sim.dis.N_tr[2])
        kx = imag(conj(ψ[1:end-1,:]).*ψ[2:end,:]/sim.dis.dx[1])
        ky = imag(conj(ψ[:,1:end-1]).*ψ[:,2:end]/sim.dis.dx[2])
        if sim.dis.N[1]==1
            bottom[i] = ky[1,1]
            top[i] = ky[1,end]
            left[i]  = 0
            right[i]  = 0
        elseif sim.dis.N[2]==1
            bottom[i] = 0
            top[i] = 0
            left[i]  = kx[1,1]
            right[i]  = kx[end,1]
        else
            slice = 2
            bottom[i] = sum(ky[slice+1:end-slice+1, slice]              )*sim.dis.dx[1]
            top[i] =    sum(ky[slice+1:end-slice+1, end-slice+1]        )*sim.dis.dx[1]
            left[i]  =  sum(kx[slice              , slice+1:end-slice+1])*sim.dis.dx[2]
            right[i]  = sum(kx[end-slice+1        , slice+1:end-slice+1])*sim.dis.dx[2]
        end
        flux[i] = top[i] - bottom[i] + right[i] - left[i]
    end

    return flux, (left,right,bottom,top)
end


#@code_warntype checked
"""
    absorption = compute_loss(sim, k, ψ)
"""
function bulk_absorption(sim::Simulation, k, ψ)

    absorption = Array{Float64}(undef, length(k))

    for i ∈ eachindex(k)
        A = imag(sim.sys.ε[:]*k[i]^2).*abs2.(ψ[:,i])
        absorption[i] = quadrature(sim, A)
    end

    return absorption
end
