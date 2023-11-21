include("../src/ColdAtoms.jl")
using .ColdAtoms
using QuantumOptics
using PyPlot
using BenchmarkTools


#Default simulation parameters
include("../params/default.jl")

N = 100;
samples, acc_rate = samples_generate(trap_params, atom_params, N; skip=5000, freq=1000);

Ωr = 2π * 35.0;
red_laser_params = [Ωr, wr, zr];
detuning_params = [Δ0, δ_twophoton(Ωr, Ωb, Δ0)];

T0 = T_twophoton(Ωr, Ωb, Δ0);
tspan = [0.0:T0/20:2.5*T0;];
ψ0 = g;

@time simulation(
    tspan, ψ0, 
    
    atom_params,
    trap_params,
    samples,
    
    f,
    red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes,
    
    red_laser_params,
    blue_laser_params,
    
    detuning_params,
    decay_params
    );


@time simulation(
    tspan, ψ0, 
    
    atom_params,
    trap_params,
    samples,
    
    f,
    red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes,
    
    red_laser_params,
    blue_laser_params,
    
    detuning_params,
    decay_params
    );