#TODO: higher order gaps

################################################################################
### DISPERSION MANAGERS
################################################################################
#@code_warntype checked
"""
    build_dispersions!(sim; num_wg_bands=15, num_bloch=17, num_free_bands=2)
"""
function build_dispersions!(sim::Simulation;
    interpolation=:cubic, β_start=0, β_stop=20, nβ=100,
    num_wg_modes=maximum([sim.sct.channels[i].quantum_number for i ∈ eachindex(sim.sct.channels)])+1,
    free_band_zone=:half, parallel=nprocs()>1, num_bloch=17, num_free_bands=2, num_wg_bands_multiplier=num_free_bands)

    num_wg = length(sim.sys.waveguides)
    if isempty(sim.sct.channels)
        throw(ErrorException("no channels specified"))
    end
    for i ∈ eachindex(sim.sys.waveguides)
        ch_inds = findall([sim.sct.channels[j].waveguide==sim.sys.waveguides[i] for j ∈ eachindex(sim.sct.channels)])
        if any(isempty.([sim.sct.channels[j].dispersion for j ∈ ch_inds]))
            println("Building dispersion for waveguide $i")

            wg_bands, ks, bands, gaps = waveguide_dispersion(sim, sim.sys.waveguides[i];
            interpolation=interpolation, β_start=β_start, β_stop=β_stop, nβ=nβ,
            num_wg_modes=num_wg_modes, free_band_zone=free_band_zone, parallel=parallel,
            num_free_bands=num_free_bands, num_bloch=num_bloch, num_wg_bands_multiplier=num_wg_bands_multiplier)

            if isempty(wg_bands)
                throw(ErrorException("no guided modes found. Either there are none in the first bandgap, or increase num_wg_modes from $num_wg_modes or num_wg_bands_multiplier from $num_wg_bands_multiplier"))
            end
            for j ∈ ch_inds
                if isempty(sim.sct.channels[j].dispersion)
                    if length(wg_bands) < sim.sct.channels[j].quantum_number
                        throw(ErrorException("channel $j (waveguide $i) has higher mode number than computed. ensure the higher mode exists in the regime you are working."))
                    end
                    push!(sim.sct.channels[j].dispersion, wg_bands[sim.sct.channels[j].quantum_number])
                    if !isempty(gaps[1])
                        for l ∈ eachindex(gaps[1])
                            push!(sim.sct.channels[j].gaps, [gaps[1][l], gaps[2][l]])
                        end
                    else
                        nothing
                    end
                end
            end
        end
    end
    return nothing
end

#@code_warntype checked
"""
    remove_dispersions!(sim)

removes dispersion curves from simulation
"""
function remove_dispersions!(sim::Simulation)
    for i ∈ 1:length(sim.sct.channels)
        for j ∈ 1:length(sim.sct.channels[i].dispersion)
            pop!(sim.sct.channels[i].dispersion)
        end
    end
    return nothing
end


"""
    waveguide_dispersion(sim, waveguide)
"""
function waveguide_dispersion(sim::Simulation, waveguide::Int;
    interpolation=:cubic, β_start=0, β_stop=20, nβ=31,
    num_wg_modes=maximum([sim.sct.channels[i].quantum_number for i ∈ eachindex(sim.sct.channels)])+1,
    free_band_zone=:half, parallel=nprocs()>1, num_bloch=17, num_free_bands=2,
    num_wg_bands_multiplier=2, disp_opt=true)

    if isempty(sim.sys.waveguides)
        nothing
    elseif ishalfspace(sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, waveguide)][1])
        wg_dispersion, ks, bands, gaps = ScalarFDFD.halfspace_waveguide_dispersion(sim, waveguide)
    elseif any(ispc.(sim.sys.domains[get_waveguide_domains(sim,waveguide)]))
        wg_dispersion, ks, bands, gaps = ScalarFDFD.pc_waveguide_dispersion(sim, waveguide;
            num_wg_bands_multiplier=num_wg_bands_multiplier, num_bloch=num_bloch,
            num_free_bands=num_free_bands, parallel=parallel, free_band_zone=free_band_zone,
            interpolation=interpolation, disp_opt=disp_opt)
    elseif ScalarFDFD.isplanar(sim.sys.domains[ScalarFDFD.get_waveguide_domains(sim, waveguide)[1]])
        wg_sim = extract_waveguide_simulation(sim, waveguide)
        wg_dispersion, ks, bands, gaps = ScalarFDFD.planar_waveguide_dispersion(wg_sim;
            nβ=nβ, num_wg_modes=num_wg_modes, β_start=β_start, β_stop=β_stop, interpolation=interpolation)
    else
        throw(ArgumentError("unrecognized waveguide type"))
    end

    return wg_dispersion, ks, bands, gaps
end


################################################################################
### WAVEGUIDE DISPERSIONS
################################################################################
"""
    wg_dispersion, bands, gaps, ks = trivial_waveguide_dispersion(sim, waveguide, num_bloch)
"""
function halfspace_waveguide_dispersion(sim::Simulation, waveguide)

    wg_dispersion = Array{AbstractInterpolation}(undef,1)

    ks = 0 .+ 0:1
    wg_inds = ScalarFDFD.get_waveguide_domains(sim, waveguide)

    bands = Array{Float64}(undef, 2, 1)
    bands[:,1] = sim.sys.domains[wg_inds[1]].n₁*ks

    wg_dispersion[1] = extrapolate(interpolate(bands[:,1], BSpline(Linear())), Line())
    wg_dispersion[1] = scale(wg_dispersion[1], ks)

    return wg_dispersion, ks, Array{Float64}(undef,0,0), (Float64[], Float64[])
end


"""
    wg_bands = pc_waveguide_dispersion(sim::Simulation, waveguide, gaps; num_free_bands=2, num_wg_bands=15, num_bloch=17, interpolation=:cubic)
"""
function pc_waveguide_dispersion(sim::Simulation, waveguide::Int; num_wg_bands_multiplier,
    num_bloch, num_free_bands, interpolation, parallel, free_band_zone, disp_opt=true)

    wg_sim = extract_waveguide_simulation(sim, waveguide)
    num_wg_bands=round(Int,num_wg_bands_multiplier*max(wg_sim.lat.a/wg_sim.lat.b,wg_sim.lat.b/wg_sim.lat.a))
    pc_sim = Simulation(wg_sim.sys.domains[2]; bnd=Boundary(bc=:p), dis=wg_sim.dis)

    if wg_sim.lat.a < wg_sim.lat.b
        ka_bloch = -π/wg_sim.lat.a:2π/wg_sim.lat.a/(num_bloch-1):π/wg_sim.lat.a
        kb_bloch = fill(0,length(ka_bloch))
        k_bloch = ka_bloch
    else
        kb_bloch = -π/wg_sim.lat.b:2π/wg_sim.lat.b/(num_bloch-1):π/wg_sim.lat.b
        ka_bloch = fill(0,length(kb_bloch))
        k_bloch = kb_bloch
    end

    _, gaps = band_structure(pc_sim, num_bloch; num_bands=num_free_bands, parallel=parallel, interpolation=interpolation, zone=free_band_zone, calc_type="waveguide band structure ", disp_opt=disp_opt)
    free_bands, _, _ = band_structure(pc_sim, (ka_bloch, kb_bloch); num_bands=num_free_bands, parallel=parallel, interpolation=interpolation, pg=Progress(Int(1e5), 1e5), disp_opt=disp_opt)
    _, _, bands = band_structure(wg_sim, (ka_bloch, kb_bloch); num_bands=num_wg_bands, parallel=parallel, calc_type="waveguide dispersion ", disp_opt=disp_opt)

    mode_bool = fill(false, num_wg_bands)
    for i ∈ eachindex(bands)
        # this is a crude way to determine guided modes TODO
        gap_size = gaps[2][1]-gaps[1][1]
        mode_bool[i] = any(gaps[1][1]+.1*gap_size .< bands[i].(1:.01:num_bloch) .< gaps[2][1]-.1*gap_size)
    end
    mode_inds = findall(mode_bool)

    wg_dispersion = Array{AbstractInterpolation}(undef,length(mode_inds))
    for j in 1:length(mode_inds)
        wg_dispersion[j] = scale(bands[mode_inds[j]], k_bloch)
    end
    return wg_dispersion, k_bloch, free_bands, gaps
end


#@code_warntype check
"""
    wg_dispersion, ks, bands, gaps = planar_waveguide_dispersion(sim, nβ, num_wg_modes, k_start, k_stop, interpolation)
"""
function planar_waveguide_dispersion(sim::Simulation; nβ, num_wg_modes, β_start, β_stop, interpolation)

    K = Array{ComplexF64}(undef, nβ, num_wg_modes)
    ks = β_start .+ (0:nβ-1)*(β_stop-β_start)/(nβ-1)
    for i ∈ eachindex(ks)
        K[i,:] = ScalarFDFD.planar_kl(sim, ks[i], .001, num_wg_modes)[1]
    end

    wg_dispersion = Array{AbstractInterpolation}(undef, num_wg_modes)
    if interpolation == :linear
        interpolator = BSpline(Linear())
    elseif interpolation == :quadratic
        interpolator = BSpline(Quadratic(Line(OnGrid())))
    elseif interpolation == :cubic
        interpolator = BSpline(Cubic(Line(OnGrid())))
    else
        throw(ArgumentError("invalid order $(interpolation)"))
    end
    for i ∈ eachindex(wg_dispersion)
        wg_dispersion[i] = extrapolate(interpolate(real(K[:,i]), interpolator), Line())
        wg_dispersion[i] = scale(wg_dispersion[i], ks)
    end

    return wg_dispersion, ks, Array{Float64}(undef,0,0), (Float64[], Float64[])
end
