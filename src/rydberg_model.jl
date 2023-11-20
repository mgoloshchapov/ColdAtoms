#Basis states
const basis = NLevelBasis(4);
const g = nlevelstate(basis, 1);
const p = nlevelstate(basis, 2);
const r = nlevelstate(basis, 3);
const gt = nlevelstate(basis, 4);

#Operators
const σgp = g ⊗ dagger(p);
const σpg = p ⊗ dagger(g);

const σpr = p ⊗ dagger(r);
const σrp = r ⊗ dagger(p);

const np = p ⊗ dagger(p);
const nr = r ⊗ dagger(r);

const σgtp = gt ⊗ dagger(p);
const σpgt = p ⊗ dagger(gt);


#Due to atom dynamics
"""
laser_params: [Ω₀, w₀, z₀]
"""
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
function δ(vz, red_laser_params, blue_laser_params; parallel=true)
    Ωr0, wr0, zr0 = red_laser_params;
    Ωb0, wb0, zb0 = blue_laser_params;
    
    kr = 2 * zr0/wr0^2;
    kb = 2 * zb0/wb0^2;

    if parallel
        return (kr + kb) * vz
    else
        return (kr - kb) * vz
    end;
end;



#Two-photon Rydberg hamiltonian for 1 atom
function Hamiltonian(Ωr, Ωb, Δ, δ)
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
"""
decay_params: [Γg, Γgt]
"""
function JumpOperators(decay_params)
    Γg, Γgt = decay_params;
    return [sqrt(Γg)*σgp, sqrt(Γgt)*σgtp], [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt]
end;


function simulation(
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
    decay_params;

    atom_motion=true,
    laser_noise=true,
    spontaneous_decay=true,
    parallel=false
    )

    N = length(samples);

    ωr, ωz = trap_frequencies(atom_params, trap_params);
    Δ0, δ0 = detuning_params;

    if spontaneous_decay
        decay_params_temp = decay_params;
    else
        decay_params_temp = [0.0, 0.0];
    end;

    J, Jdagger = JumpOperators(decay_params_temp);

    ρ0 = ψ0 ⊗ dagger(ψ0);


    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    #Second moment for error estimation of level populations. 
    #Not sure if I can use it for error estimation of arbitrary operators.
    ρ2_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    for i ∈ 1:N
        if atom_motion
            #Atom initial conditions
            xi, yi, zi, vxi, vyi, vzi = samples[i];
        else
            xi, yi, zi, vxi, vyi, vzi = zeros(6);
        end;
        
        #Atom trajectories
        X(t) = R(t, xi, vxi, ωr);
        Y(t) = R(t, yi, vyi, ωr);
        Z(t) = R(t, zi, vzi, ωz);
        Vz(t) = V(t, zi, vzi, ωz);

        
        if laser_noise
            red_laser_phase_amplitudes_temp = red_laser_phase_amplitudes;
            blue_laser_phase_amplitudes_temp = blue_laser_phase_amplitudes;
        else
            red_laser_phase_amplitudes_temp = zero(red_laser_phase_amplitudes);
            blue_laser_phase_amplitudes_temp = zero(blue_laser_phase_amplitudes);
        end;

        #Generate phase noise traces for red and blue lasers
        tspan_noise = [0.0:tspan[end]/1000:tspan[end];];
        ϕ_red_res = ϕ(tspan_noise, f, red_laser_phase_amplitudes_temp);
        ϕ_blue_res = ϕ(tspan_noise, f, blue_laser_phase_amplitudes_temp);

        #Interpolate phase noise traces to pass to hamiltonian
        nodes = (tspan_noise, );
        ϕ_red = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
        ϕ_blue = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));

        
        #Hamiltonian params trajectories
        δ_temp(t) = δ(Vz(t), red_laser_params, blue_laser_params; parallel=parallel) + δ0;
        Δ_temp(t) = Δ(Vz(t), red_laser_params) + Δ0;
        Ωr_temp(t) = exp(1.0im * ϕ_red(t)) * Ω(X(t), Y(t), Z(t), red_laser_params);
        Ωb_temp(t) = exp(1.0im * ϕ_blue(t)) * Ω(X(t), Y(t), Z(t), blue_laser_params);
        
        #Hamiltonian
        H_temp = TimeDependentSum(
        [
            t -> -Δ_temp(t),
            t -> -δ_temp(t),
            t -> Ωr_temp(t)/2.0,
            t -> conj(Ωr_temp(t))/2.0,
            t -> Ωb_temp(t)/2.0,
            t -> conj(Ωb_temp(t))/2.0
        ],
        
        [
            np,
            nr,
            σgp,
            σpg,
            σpr,
            σrp  
        ]
        );
        
        #Returns hamiltonian and jump operators in a form required by timeevolution.master_dynamic
        function super_operator(t, rho)
            return H_temp, J, Jdagger;
        end;
        

        tout, ρ = timeevolution.master_dynamic(tspan, ρ0, super_operator);

        ρ_mean = ρ_mean + ρ;
        ρ2_mean = ρ2_mean + ρ .^ 2;
    end;

    return ρ_mean/N, ρ2_mean/N
end;


function Ω_twophoton(Ωr, Ωb, Δ)
    return Ωb * Ωr / (2.0 * Δ)
end;


function T_twophoton(Ωr, Ωb, Δ)
    return 2.0*π / Ω_twophoton(Ωr, Ωb, Δ);
end;


function δ_twophoton(Ωr, Ωb, Δ)
    return (Ωb^2 - Ωr^2)/(4.0 * Δ)
end;


function Ωr_required(Ω, Ωb, Δ)
    return 2.0 * Δ * Ω / Ωb
end;



using Distributed
function simulation_parallel(
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
    decay_params;

    atom_motion=true,
    laser_noise=true,
    spontaneous_decay=true,
    parallel=false
    )

    N = length(samples);

    ωr, ωz = trap_frequencies(atom_params, trap_params);
    Δ0, δ0 = detuning_params;

    if spontaneous_decay
        decay_params_temp = decay_params;
    else
        decay_params_temp = [0.0, 0.0];
    end;

    J, Jdagger = JumpOperators(decay_params_temp);

    ρ0 = ψ0 ⊗ dagger(ψ0);


    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    ρ_temp = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    #Second moment for error estimation of level populations. 
    #Not sure if I can use it for error estimation of arbitrary operators.
    ρ2_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    @distributed for i ∈ 1:N
        if atom_motion
            #Atom initial conditions
            xi, yi, zi, vxi, vyi, vzi = samples[i];
        else
            xi, yi, zi, vxi, vyi, vzi = zeros(6);
        end;
        
        #Atom trajectories
        X(t) = R(t, xi, vxi, ωr);
        Y(t) = R(t, yi, vyi, ωr);
        Z(t) = R(t, zi, vzi, ωz);
        Vz(t) = V(t, zi, vzi, ωz);

        
        if laser_noise
            red_laser_phase_amplitudes_temp = red_laser_phase_amplitudes;
            blue_laser_phase_amplitudes_temp = blue_laser_phase_amplitudes;
        else
            red_laser_phase_amplitudes_temp = zero(red_laser_phase_amplitudes);
            blue_laser_phase_amplitudes_temp = zero(blue_laser_phase_amplitudes);
        end;

        #Generate phase noise traces for red and blue lasers
        tspan_noise = [0.0:tspan[end]/1000:tspan[end];];
        ϕ_red_res = ϕ(tspan_noise, f, red_laser_phase_amplitudes_temp);
        ϕ_blue_res = ϕ(tspan_noise, f, blue_laser_phase_amplitudes_temp);

        #Interpolate phase noise traces to pass to hamiltonian
        nodes = (tspan_noise, );
        ϕ_red = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
        ϕ_blue = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));

        
        #Hamiltonian params trajectories
        δ_temp(t) = δ(Vz(t), red_laser_params, blue_laser_params; parallel=parallel) + δ0;
        Δ_temp(t) = Δ(Vz(t), red_laser_params) + Δ0;
        Ωr_temp(t) = exp(1.0im * ϕ_red(t)) * Ω(X(t), Y(t), Z(t), red_laser_params);
        Ωb_temp(t) = exp(1.0im * ϕ_blue(t)) * Ω(X(t), Y(t), Z(t), blue_laser_params);
        
        #Hamiltonian
        H_temp = TimeDependentSum(
        [
            t -> -Δ_temp(t),
            t -> -δ_temp(t),
            t -> Ωr_temp(t)/2.0,
            t -> conj(Ωr_temp(t))/2.0,
            t -> Ωb_temp(t)/2.0,
            t -> conj(Ωb_temp(t))/2.0
        ],
        
        [
            np,
            nr,
            σgp,
            σpg,
            σpr,
            σrp  
        ]
        );
        
        #Returns hamiltonian and jump operators in a form required by timeevolution.master_dynamic
        function super_operator(t, rho)
            return H_temp, J, Jdagger;
        end;
        

        tout, ρ_temp = timeevolution.master_dynamic(tspan, ρ0, super_operator);

        ρ_mean = ρ_mean + ρ_temp;
        ρ2_mean = ρ2_mean + ρ_temp .^ 2;
    end;

    return ρ_mean/N, ρ2_mean/N
end;


function Ω_twophoton(Ωr, Ωb, Δ)
    return Ωb * Ωr / (2.0 * Δ)
end;


function T_twophoton(Ωr, Ωb, Δ)
    return 2.0*π / Ω_twophoton(Ωr, Ωb, Δ);
end;


function δ_twophoton(Ωr, Ωb, Δ)
    return (Ωb^2 - Ωr^2)/(4.0 * Δ)
end;


function Ωr_required(Ω, Ωb, Δ)
    return 2.0 * Δ * Ω / Ωb
end;