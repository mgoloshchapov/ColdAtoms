module ColdAtoms

using PyPlot
using Statistics, Distributions, Random
using PhysicalConstants.CODATA2018: c_0, k_B, m_u
using Unitful
using LinearAlgebra
using QuantumOptics
using DifferentialEquations
using SplitApplyCombine
using Interpolations

export 
    w0_to_z0, trap_frequencies,
    release_recapture,
    samples_generate, samples_visualise, R, V, get_trap_params,
    Sϕ, ϕ_amplitudes, ϕ,
    simulation, Ω_twophoton, T_twophoton, δ_twophoton, Ωr_required, 
    g, p, r, gt
    
        
include("utilities.jl")
include("basic_experiments.jl")
include("lasernoise_sampler.jl")
include("atom_sampler.jl")
include("rydberg_model.jl")
include("plot_rus.jl")

end