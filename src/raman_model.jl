
#Basis states
const basis = NLevelBasis(2);
const g = nlevelstate(basis, 1);
const e = nlevelstate(basis, 2);

#Operators
const σge = g ⊗ dagger(e);
const σeg = e ⊗ dagger(g);

const operators = [σge, σeg];


#Due to atom dynamics
"""
laser_params: [Ω₀, w₀, z₀]
"""
function Ω(x, y, z, laser_params)
    Ω0, w0, z0 = laser_params;
    return Ω0 .* I(x, y, z, w0, z0);
end;


function simulation_raman(
    tspan, ψ0, 
        
    atom_params,
    laser_params,
    trap_params,

    samples,
    phase;
    distance=0.0
    )
    """
    atom_params:        [m(a.u.), T(μK)]
    laser_params:       [Ωr(MHz), wr(μm), zr(μm)]
    trap_params:        [U0(μK), w0(μm), z0(μm)]
    samples:            samples of atom initial coordinates and velocities
    phase:              function of time phase(t)
    """
    N = length(samples);

    ωr, ωz = trap_frequencies(atom_params, trap_params);

    # ρ0 = ψ0 ⊗ dagger(ψ0);

    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    ρ_temp = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    #Second moment for error estimation of level populations. 
    ρ2_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    
    for i ∈ 1:N

        #Atom initial conditions
        xi, yi, zi, vxi, vyi, vzi = samples[i];
        
        #Atom trajectories
        X = t -> R(t, xi, vxi, ωr; free=false) - distance;
        Y = t -> R(t, yi, vyi, ωr; free=false);
        Z = t -> R(t, zi, vzi, ωz; free=false);


        #Hamiltonian params trajectories
        Ht = TimeDependentSum(
        [
            t -> exp(1.0im * phase(t)) * Ω(X(t), Y(t), Z(t), laser_params) / 2.0;
            t -> conj(exp(1.0im * phase(t)) * Ω(X(t), Y(t), Z(t), laser_params) / 2.0);
        ],
        operators
        );

        _, ψ_temp = timeevolution.schroedinger_dynamic(tspan, ψ0, Ht);
        ρ_temp = ψ_temp .⊗ dagger.(ψ_temp);

        ρ_mean = ρ_mean + ρ_temp;
        ρ2_mean = ρ2_mean + ρ_temp .^ 2;
    end;

    return ρ_mean/N, ρ2_mean/N
end;
