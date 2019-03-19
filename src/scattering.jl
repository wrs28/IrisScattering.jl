# TODO: nonlinear scattering
################################################################################
####### SCATTERING
################################################################################
"""
    ψ_sct, ψ_tot, H = scattering(sim, k, a; H=lu(sparse(I,1,1)), F=[1], file_name="",
        is_linear=true, disp_opt=false, ψ_init=[], ftol=2e-8, iter=150, num_wg_bands=15)

Solves linear inhomogeneous problem.

`k` is the frequency.

`a` is a vector of amplitudes for each channel given in input.

`H` is a factorized wave operator (useful for repeated calculation)

`F` is a multiplicative modifier to the pump profile (can implement saturation this way)

if `file_name` is not empty, output is saved to file

`num_bloch` is the number of bloch frequencies which are used to determine gaps

`num_wg_bands` is the number of bands in the presence of the waveguide to use t
look for guided modes

Note: use build_dispersions! for more control over the dispersion-building process.
"""
function scattering(sim::Simulation, k, a;
        H=lu(sparse(complex(1.,0)*I,1,1)),
        file_name="", is_linear=true, ψ_init=[], ftol=2e-8, iter=150)

    if length(a) < length(sim.sct.channels)
        throw(ArgumentError("number of input amplitudes $(length(a)) less than the number of input channels $(length(sim.sct.channels))"))
    elseif length(a) > length(sim.sct.channels)
        @warn "number of input amplitudes $(length(a)) greater than the number of input channels $(length(sim.sct.channels)), will ignore last $(length(a)-length(sim.sct.channels)) amplitudes"
    end

    build_dispersions!(sim)

    # if is_linear
        ψ, φ₊₋, φ₊, φ₋, H = scattering_l(sim, k, a, H)
    # else
        # ψ, (φ₊₋, φ₊, φ₋, H) = scattering_nl(sim, k, a, H, F, disp_opt=disp_opt, ψ_init=ψ_init, ftol=ftol, iter=iter)
    # end

    if !isempty(file_name)
        fid = open(file_name, "w")
        serialize(fid, (sim, k, ψ) )
        close(fid)
    end

    ψ_sct = ψ
    ψ_tot = ψ + φ₊₋ + φ₊ + φ₋
    ψ_out = ψ + φ₋

    return ψ_sct, ψ_tot, ψ_out, H
end


function scattering_l(sim::Simulation, k, a, H)
    j, φ₊₋, φ₊, φ₋ = ScalarFDFD.synthesize_source(sim, k, a)
    if (H.m == H.n == 1)
        H = ScalarFDFD.factorize_scattering_operator(sim,k)
    end
    ψ = H\j
    return ψ, φ₊₋, φ₊, φ₋, H
end


function factorize_scattering_operator(sim::Simulation, k)
    k²= k^2
    bl_original = ScalarFDFD.set_bl!(sim, :pole)
    try
        ∇² = laplacian(sim, k)
        N = prod(sim.dis.N)
        εt = ε_bl(sim; k=k)

        ɛk² = sparse(1:N, 1:N, ɛt*k², N, N)
        # χk² = sparse(1:N, 1:N, float(sim.tls.D₀)*γ(sim,k)*F.*Ft[:]*k², N, N)

        return lu(∇²+ɛk²)
    finally
        reset_bl!(sim, bl_original)
    end
end











# ################################################################################
# ### NONLINEAR SOLVERS
# ################################################################################
# """
# ψ, ϕ, A, input = scattering_nl(input, k; dispOpt=false, ψ_init=[],
#     F=[1.], A=[], fileName="", ftol=2e-8, iter=750)
#
#     Solves inhomogeneous problem with source defined in incident wave file.
#
#     k is the injection frequency.
#
#     A is a factorized wave operator.
# """
# function scattering_nl(input1::InputStruct, k::Complex128;
#     dispOpt::Bool = false, ψ_init::Array{Complex128,1}=Complex128[], F::Array{Float64,1}=[1.],
#     A::Base.SparseArrays.UMFPACK.UmfpackLU=lufact(speye(1,1)),
#     ftol::Float64=2e-8, iter::Int=150)::Tuple{Array{Complex128,1},Array{Complex128,1},Base.SparseArrays.UMFPACK.UmfpackLU,InputStruct}
#
#         input = open_to_pml(input1)
#
#         j, ∇², φ₊, φ₋ = createJ(input, k, m)
#
#         N = prod(input.N_ext); ε_sm = input.ε_sm; k²= k^2;
#         D₀ = input.D₀; F_sm = input.F_sm
#         ɛk² = sparse(1:N, 1:N, ɛ_sm[:]*k², N, N, +)
#         χk² = sparse(1:N, 1:N, D₀*γ(k,input)*F.*F_sm[:]*k², N, N, +)
#
#         f!(Ψ, fvec) = scattering_residual(Ψ, fvec, j, ∇², εk², k, input)
#         jac!(Ψ, jacarray) = scattering_jacobian(Ψ, jacarray, j, ∇², εk², k, input)
#         df = DifferentiableSparseMultivariateFunction(f!, jac!)
#
#         Ψ_init = Array{Float64}(2*length(j))
#         Ψ_init[1:length(j)]     = real(ψ_init)
#         Ψ_init[length(j)+1:2*length(j)] = imag(ψ_init)
#
#         z = nlsolve(df, Ψ_init, show_trace=dispOpt, ftol=ftol, iterations=iter)
#
#         if converged(z)
#             ψ = z.zero[1:length(j)] + 1im*z.zero[length(j)+1:2*length(j)]
#         else
#             ψ = NaN*ψ
#             println("Warning, solve_scattered did not converge. Returning NaN.")
#         end
#
#         if !isempty(fileName)
#             if truncate
#                 fid = open(fileName,"w")
#                 serialize(fid, ((ψ-φ₊)[input.x̄_inds], φ₊[input.x̄_inds], input, k) )
#                 close(fid)
#             else
#                 fid = open(fileName,"w")
#                 serialize(fid, ((ψ-φ₊), φ₊, input, k) )
#                 close(fid)
#             end
#         end
#
#         if truncate
#             return (ψ-φ₊)[input.x̄_inds], φ₊[input.x̄_inds], A, input
#         else
#             return ψ-φ₊, φ₊, A, input
#         end
# end # end of function computePsi_nonlinear
#
# ################################################################################
# ###  NONLINEAR AUXILLIARIES
# ################################################################################
# """
# scattering_residual(Ψ, fvec, j, ∇², εk², k, input)
# """
# function scattering_residual(Ψ::Array{Float64,1}, fvec::Array{Float64,1}, j,
#                                     ∇², εk², k::Complex128, input::InputStruct)
#
#     ψ = similar(j,Complex128)
#     ind_r = 1:length(j)
#     ind_i = length(j)+1:2*length(j)
#     ψ = Ψ[ind_r] + 1im*Ψ[ind_i]
#     temp = (∇²+ɛk²+χ(ψ,k,input)*k^2)*ψ - j
#     fvec[ind_r] = real(temp)
#     fvec[ind_i] = imag(temp)
# end
#
# """
# scattering_jacobian(Ψ, jacarray, j, ∇², εk², k, input)
# """
# function scattering_jacobian(Ψ::Array{Float64,1}, jacarray, j, ∇², εk², k::Complex128, input::InputStruct)
#     ψ = similar(j,Complex128)
#     ind_r = 1:length(j)
#     ind_i = length(j)+1:2*length(j)
#     ψ = Ψ[ind_r] + 1im*Ψ[ind_i]
#     temp = ∇²+ɛk²+χ(ψ,k,input)*k^2
#     tempr = similar(temp,Float64)
#     tempi = similar(temp,Float64)
#     tr = nonzeros(tempr)
#     ti = nonzeros(tempi)
#     tr[:] = real((nonzeros(temp)))
#     ti[:] = imag((nonzeros(temp)))
#     tempj = [tempr+χʳʳ′(ψ,k,input) -tempi+χʳⁱ′(ψ,k,input); tempi+χⁱʳ′(ψ,k,input) tempr+χⁱⁱ′(ψ,k,input)]
#     jacarray[:,:] = tempj[:,:]
# end
#
# """
# χ(ψ, k, input)
# """
# function χ(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Complex{Float64},Int64}
#     N = prod(input.N_ext)
#     h = hb(ψ,k,input)
#     V = input.F_sm[:].*γ(k,input)*input.D₀./(1+h)
#     return sparse(1:N, 1:N, V, N, N, +)
# end
#
# """
# """
# function hb(ψ::Array{Complex128,1},k::Complex128,input::InputStruct)::Array{Float64,1}
#     N = prod(input.N_ext)
#     h = abs2.(γ.(k, input)*ψ)
#     return h
# end
#
# """
# χʳʳ′(Ψ, k, input)
# """
# function χʳʳ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
#     N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
#     h = hb(ψ,k,input)
#     V = -2D₀.*F_sm[:].*abs2(γt).*real.(γt.*ψ).*real.(ψ)./(1+h).^2
#     return sparse(1:N,1:N,V,N,N,+)
# end
#
# """
# χⁱʳ′(Ψ, k, input)
# """
# function χⁱʳ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
#     N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
#     h = hb(ψ,k,input)
#     V = -2D₀.*F_sm[:].*abs2(γt).*imag.(γt.*ψ).*real.(ψ)./(1+h).^2
#     return sparse(1:N,1:N,V,N,N,+)
# end
#
# """
# χʳⁱ′(Ψ, k, input)
# """
# function χʳⁱ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
#     N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
#     h = hb(ψ,k,input)
#     V = -2D₀.*F_sm[:].*abs2(γt).*real.(γt.*ψ).*imag.(ψ)./(1+h).^2
#     return sparse(1:N,1:N,V,N,N,+)
# end
#
# """
# χⁱⁱ′(Ψ, k, input)
# """
# function χⁱⁱ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
#     N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
#     h = hb(ψ,k,input)
#     V = -2D₀.*F_sm[:].*abs2(γt).*imag.(γt.*ψ).*imag.(ψ)./(1+h).^2
#     return sparse(1:N,1:N,V,N,N,+)
# end
