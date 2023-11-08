include("utilities.jl")
using QuantumOptics
using ProgressMeter
using BenchmarkTools


#Basis states
basis = NLevelBasis(4);
g = nlevelstate(basis, 1);
p = nlevelstate(basis, 2);
r = nlevelstate(basis, 3);
gt = nlevelstate(basis, 4);

#Operators
σgp = g ⊗ dagger(p);
σpg = p ⊗ dagger(g);

σpr = p ⊗ dagger(r);
σrp = r ⊗ dagger(p);

np = p ⊗ dagger(p);
nr = r ⊗ dagger(r);

σgtp = gt ⊗ dagger(p);
σpgt = p ⊗ dagger(gt);


#Due to atom dynamics
function Ω(x, y, z, laser_params)
    Ω0, w0, z0 = laser_params;
    return Ω0 .* A(x, y, z, w0, z0) .* A_phase(x, y, z, w0, z0);
end;


#Due to Doppler shift for red laser
function Δ(vz, laser_params)
    Ω0, w0, z0 = laser_params
    k = 2 * z0/w0^2;
    return k * vz
end;



#Due to Doppler shifts for red and blue lasers
function δ(vz, red_laser_params, blue_laser_params)
    Ωr0, wr0, zr0 = red_laser_params;
    Ωb0, wb0, zb0 = blue_laser_params;
    
    kr = 2 * zr0/wr0^2;
    kb = 2 * zb0/wb0^2;
    
    return (kr - kb) * vz
end;


#Two-photon Rydberg hamiltonian for 1 atom
function Hamiltonian(t, Ωr, Ωb, Δ, δ)
    return TimeDependentSum(
        [
        t -> -Δ(t),
        t -> -δ(t),
        t -> Ωr(t) ./2.0,
        t -> conj.(Ωr(t)) ./2.0,
        t -> Ωb(t)/2.0,
        t -> conj.(Ωb(t)) ./2.0,
        ],
        
        [
        np,
        nr,
        σgp,
        σpg,
        σpr,
        σrp  
        ]
    )
end;



#Jump operators for master equation 
function JumpOperators(decay_params)
    Γg, Γgt = decay_params;
    
    return [sqrt(Γg)*σgp, sqrt(Γgt)*σgtp], [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt]
end;



function simulation_stable(
    tspan, ρ0, 
    
    samples,
    
    f,
    red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes,
    
    red_laser_params,
    blue_laser_params,
    
    detuning_params,
    decay_params;
    
    laser_noise=true,
    atom_dynamics=true
)

N = length(samples);
ωr, ωz = trap_frequencies(atom_params, trap_params);

Δ0, δ0 = detuning_params;
J, Jdagger = JumpOperators(decay_params);

ρ_mean = [zero(ρ0) for _ ∈ 1:length(tspan)];

@showprogress for i ∈ 1:N
    #Atom initial conditions
    xi, yi, zi, vxi, vyi, vzi = samples[i];
    

    #Atom trajectories
    X(t) = R(t, xi, vxi, ωr);
    Y(t) = R(t, yi, vyi, ωr);
    Z(t) = R(t, zi, vzi, ωz);
    Vz(t) = V(t, zi, vzi, ωz);


    #Generate phase noise traces for red and blue lasers
    ϕ_red_res = ϕ(tspan, f, red_laser_phase_amplitudes);
    ϕ_blue_res = ϕ(tspan, f, blue_laser_phase_amplitudes);

    #Interpolate phase noise traces to pass to hamiltonian
    nodes = (tspan, );
    ϕ_red = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
    ϕ_blue = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));
    
    #Hamiltonian params trajectories
    δ_temp(t) = δ(Vz(t), red_laser_params, blue_laser_params) .+ δ0;
    Δ_temp(t) = Δ(Vz(t), red_laser_params) .+ Δ0;
    Ωr_temp(t) = exp.(1.0im .* ϕ_red(t)) .* Ω(X(t), Y(t), Z(t), red_laser_params);
    Ωb_temp(t) = exp.(1.0im .* ϕ_blue(t)) .* Ω(X(t), Y(t), Z(t), blue_laser_params);
    
    #Hamiltonian
    H_temp(t) = Hamiltonian(t, Ωr_temp, Ωb_temp, Δ_temp, δ_temp);
    
    #Returns hamiltonian and jump operators in a form required by timeevolution.master_dynamic
    super_operator(t, rho) = [H_temp(t), J, Jdagger];

    #Time evolution
    tout, ρ = timeevolution.master_dynamic(tspan, ρ0, super_operator;);
    
    ρ_mean = ρ_mean .+ ρ;
end;

return ρ_mean/N
end;