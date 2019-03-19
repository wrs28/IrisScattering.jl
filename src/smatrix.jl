#TODO: nonlinear smatrix routines, possibly using Optim?
"""
    S, flux_sct, flux_tot, absorption = smatrix(sim, k; channels, is_linear=true, F=[1], disp_opt=true, fileName="")

scattering matrix `S` at frequences `k` on `channels`
"""
function smatrix(sim::Simulation, k; is_linear=true, disp_opt=true, file_name="", is_parallel=nprocs()>1, num_channel_blocks=1)

    if is_linear && !is_parallel
        S, scattered_flux, total_flux, absorption = smatrix_l(sim, k, 1:length(sim.sct.channels), disp_opt, file_name)
    elseif is_linear && is_parallel
        S, scattered_flux, total_flux, absorption = smatrix_lp(sim, k, disp_opt, file_name, num_channel_blocks)
    else
        throw("no nonlinear method has yet been defined.")
    end

    return S, scattered_flux, total_flux, absorption
end


"""
    smatrix_l(sim, k, channels, F, disp_opt, file_name)
"""
function smatrix_l(sim::Simulation, k, channels, disp_opt, file_name)

    build_dispersions!(sim)

    nc = length(channels)
    nk = length(k)
    M = length(sim.sct.channels)

    S = Array{ComplexF64}(undef, nk, nc, M)
    flux_sct = Array{Float64}(undef, nk, nc)
    flux_tot = Array{Float64}(undef, nk, nc)
    absorption = Array{Float64}(undef, nk, nc)
    a = zeros(ComplexF64,M)

    if disp_opt
        pg = Progress(length(k)*M*nc, PROGRESS_UPDATE_TIME::Float64, "smatrix ")
    end

    for i ∈ eachindex(k)

        ζ = factorize_scattering_operator(sim, k[i])

        for m in 1:nc
            a = 0a
            a[channels[m]] = 1.
            ψ_sct, ψ_tot, ψ_out, ζ = scattering(sim, k[i], a; H=ζ, is_linear=true)
            for m′ ∈ 1:M
                if isempty(sim.sys.waveguides)
                    c = analyze_field(sim, k[i], ψ_out, m′)
                else
                    if sim.sct.channels[m′].waveguide == sim.sct.channels[m].waveguide
                        c = analyze_field(sim, k[i], ψ_sct, m′)
                    else
                        c = analyze_field(sim, k[i], ψ_tot, m′)
                    end
                end
                S[i, m, m′] = c
                if disp_opt
                    next!(pg)
                end
            end
            flux_sct[i, m] = surface_flux(sim, ψ_sct)[1][1]
            flux_tot[i, m] = surface_flux(sim, ψ_tot)[1][1]
            absorption[i, m] = bulk_absorption(sim, k[i], ψ_tot)[1]
        end

        if !isempty(file_name)
            fid = open(file_name,"w")
            serialize(fid, (S, flux_sct, flux_tot, absorption, sim, i))
            close(fid)
        end
    end

    return S, flux_sct, flux_tot, absorption
end


















################################################################################
### S-MATRIX LINEAR ROUTINES
################################################################################


# end of fuction smatrix_l


#
# ################################################################################
# ### S-MATRIX NONLINEAR ROUTINES
# ################################################################################
# """
# S =  smatrix_nl(input::InputStruct; N=10, N_type="D", isNonLinear=false, F=1.,
#     dispOpt = true, ψ_init = [], fileName = "")
#
#     N is the number of steps to go from D0 = 0 to given D0
# """
# function smatrix_nl(input1::InputStruct, k::Array{ComplexF64,1};
#     N::Int=1, N_type::String="D", isNonLinear::Bool=false, F::Array{Float64,1}=[1.],
#     dispOpt::Bool=true, ψ_init::Array{ComplexF64,1}=ComplexF64[],
#     fileName::String = "")::Array{ComplexF64,4}
#
#     if !isempty(input1.bc ∩ ["o", "open", "pml_in"])
#         input = deepcopy(input1)
#         for i in 1:4
#             if input.bc[i] in ["o", "open", "pml_in"]
#                 input.bc[i] = "pml_out"
#                 updateInputs!(input, :bc, input.bc)
#             end
#         end
#     else
#         input = input1
#     end
#
#     M = length(input.channels)
#     D₀ = input.D₀
#     A = input.a
#
#     ψ₋ = Array{ComplexF64}(prod(input.N_ext))
#     ψ  = Array{ComplexF64}(prod(input.N_ext))
#
#     # if isempty(ψ_init) && (N>1)
#     S = NaN*ones(ComplexF64,length(input.k),N,M,M)
#     # else
#         # S = NaN*ones(ComplexF64,length(input.k),1,2,2)
#     # end
#
#     if N > 1
#         D = linspace(0,D₀,N)
#         A_vec = linspace(0.001,1,N)
#     else
#         D = [input.D₀]
#         A_vec = [1]
#     end
#
#     for ii in 1:length(input.k)
#
#         k = input.k[ii]
#
#         if (ii/1 == round(ii/1)) & dispOpt
#             if typeof(k)<:Real
#                 printfmtln("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}.",ii,length(input.k),k)
#             else
#                 printfmtln("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}{4:+2.3f}i.",ii,length(input.k),real(k),imag(k))
#             end
#         end
#
#         if isempty(ψ_init)
#
#             for j in 1:N
#
#                 if N_type == "D"
#                     updateInputs!(input, :D₀, D[j])
#                 elseif N_type == "A"
#                     updateInputs!(input, :a, A_vec[j]*A)
#                 end
#
#                 if isNonLinear
#                     if isnan(ψ[1]) | j==1
#                         ψ, W = computePsi(input,k,isNonLinear = true, F = F)
#                     else
#                         ψ, W = computePsi(input,k,isNonLinear = true, F = F, ψ_init = ψ)
#                     end
#                 else
#                     ψ = 0.
#                 end
#
#                 # compute linearized S-matrix
#                 ζ = lufact(speye(1,1)) # initialize lufact variable (speeds up repeated inhomogeneous solves)
#                 for m1 in 1:M
#                     # set flux-normalized input amplitudes to unity
#                     at = zeros(ComplexF64,M)
#                     at[m1] = 1.
#                     updateInputs!(input,:a,at)
#                     # solve for scattered field
#                     ψ₋, dummy1, ζ, inputs_s = computePsi(input, k; isNonLinear=false, F=F./(1+abs2.(γ(k,input)*ψ)), A=ζ)
#                     # analyze into channels
#                     println("here 1")
#                     for m2 in m1:M
#                         dt=0
#                         if input.channelBoundaries[m2] in [1,2]
#                             t = inputs_s.x₂_ext
#                             u = inputs_s.x₁_inds
#                             dt = input.dx̄[2]
#                         else
#                             t = inputs_s.x₁_ext
#                             u = inputs_s.x₂_inds
#                             dt = input.dx̄[1]
#                         end
#                         println("here 2")
#                         ϕ = zeros(ComplexF64,length(t))
#                         for q in 1:length(t)
#                              ϕ[q] = conj(input.incidentWave(k, m2, input.∂R[input.channelBoundaries[m2]], t[ii], input.∂R, input.geometry, input.geoParams)[1])
#                         end
#                         println("here 3")
#                         P = 0
#                         if input.channelBoundaries[m2] in [1,3]
#                             println("here 4")
#                             println(size(ψ₋))
#                             println(size(t))
#                             println(size(u))
#                             P = reshape(ψ₋,:,length(t))[u[1],:]
#                             println("here 5")
#                         else
#                             println("here 6")
#                             P = reshape(ψ₋,:,length(t))[u[end],:]
#                         end
#                         println("here 7")
#                         println(size(P))
#                         println(size(S))
#                         println(j)
#                         println(ii)
#                         println(size(ϕ))
#                         S[ii,j,m1,m2] = sum(P.*ϕ)*dt
#                         S[ii,j,m2,m1] = S[ii,j,m1,m2]
#                     end
#                 end
#     println("here 8")
#                 updateInputs!(input,:D₀, D₀)
#                 updateInputs!(input,:a , A )
#
#                 if !isempty(fileName)
#                     fid = open(fileName,"w")
#                     serialize(fid,(input,D,S,ii,j))
#                     close(fid)
#                 end
#
#             end
#
#         else
#             if isNonLinear
#                 ψ = scattering(input,k,isNonLinear = true, F = F, ψ_init = ψ_init)
#             else
#                 ψ = 0.
#             end
#
#             input.a = [1.0,0.0]
#             ψ₊, W = scattering(input,k,isNonLinear = false, F = F./(1+abs2.(γ(input,k)*ψ)))
#
#             input.a = [0.0,1.0]
#             ψ₋, W = scattering(input,k,isNonLinear = false, F = F./(1+abs2.(γ(input,k)*ψ)), A = W)
#
#             S[1,1,ii,1] = ψ₊[x_inds[1]]*exp(+1im*input.dx*k)
#             S[2,1,ii,1] = ψ₊[x_inds[end]]
#             S[1,2,ii,1] = ψ₋[x_inds[1]]
#             S[2,2,ii,1] = ψ₋[x_inds[end]]*exp(-1im*input.dx*k)
#
#             input.D₀ = D₀
#             input.a = A
#
#             if !isempty(fileName)
#                 foo2(fid) = serialize(fid,(input,D,S,ii,1))
#                 open(foo2,fileName,"w")
#             end
#
#         end
#
#     end
#     println("here 9")
#     return S
# end # end of fuction smatrix_nl
