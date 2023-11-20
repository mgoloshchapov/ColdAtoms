using QuantumOptics
using ProgressMeter
using BenchmarkTools
using Profile, ProfileView

include("utilities.jl")
include("atom_sampler.jl")
include("rydberg_model.jl")


m = 86.9091835;    
T = 50.0;
U0 = 1000.0;
w0 = 1.0;
λ = 0.852;
z0 = π*w0^2/λ;

atom_params = [m, T];
trap_params = [U0, w0, z0];
detuning_params = [2*π * 740.0, 0.0];
red_laser_params = [2π * 14.0, 5.0, 5.0*3.68];
blue_laser_params = [2π * 14.0, 5.0, 5.0*3.68];

f = [0.01:0.01:10.0;];
red_laser_phase_amplitudes = [];
blue_laser_phase_amplitudes = [];

N = 10;
samples, acc_rate = samples_generate(trap_params, atom_params, N; freq=10, skip=1000, harmonic=false)
tspan = [0.0:1.0ew:10.0;];
ψ0 = g;


ProfileView.@profview simulation_schroedinger(
    tspan, ψ0, 
    trap_params,
    atom_params,

    red_laser_params,
    blue_laser_params,
    detuning_params,
    
    samples,
    
    f,
    red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes;

    laser_noise = false
)