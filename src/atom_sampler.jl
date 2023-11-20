#Constants and scales
#------------------------------------
const c = ustrip(u"m/s", c_0);  #Speed of light
const kB = ustrip(u"J/K", k_B)  #Boltzmann constant
const mu = ustrip(u"kg", m_u);  #Unit of atomic mass

const E0 = kB * 1e-6;       #Characteristic energy in μK
const g0 = 9.81 * 1e-6;     #Gravity free fall acceleration
const vconst = sqrt(E0/mu); #Useful constant for kinetic energy
const r0 = 1e-6;            #Characteristic distance in m
#------------------------------------


"""
Functions
"""
#Potential energy in gaussian beam
function Π(cord, trap_params)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    return U0 .* (1.0 .- A(x, y, z, w0, z0) .^2);
end;


#Potential energy in gaussian beam in harmonic approximation
function Π_Harmonic(cord, trap_params)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    
    r2 = x .^2 .+ y .^2;
    return U0 .* (2.0*r2 ./w0^2 + (z ./z0).^2);
end;



#Kinetic energy
function K(cord, trap_params, m)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    return m/vconst^2 *(vx .^2 + vy .^2 + vz .^2) / 2.0
end;



#Total energy
function H(cord, trap_params, m; harmonic=false)
    if harmonic
        return Π_Harmonic(cord, trap_params) .+ K(cord, trap_params, m)
    else
        return Π(cord, trap_params) .+ K(cord, trap_params, m)
    end;
end;


#Target distribution for Monte-Carlo
function prob_boltzmann(cord, trap_params, atom_params; harmonic=false)    
    m, T = atom_params;
    return exp.(- H(cord, trap_params, m; harmonic) ./ T)
end;


"""
trap_params: [U0(μK), w0(μm), z0(μm)]
atom_params: [m(a.u.), T(μK)]
N:            number of samples
freq:         make step equal to freq between samples to reduce correlation between them in MCMC
skip:         skip first samples, so Markov Chain can converge to desired distribution
harmonic:     make harmonic approximation, false by default
"""
function samples_generate(trap_params, atom_params, N; freq=10, skip=1000, harmonic=false)
    U0, w0, z0 = trap_params;
    m, T = atom_params;

    mean = zeros(6);
    vstep = vconst*sqrt(T/m);
    rstep = sqrt(T/U0)/2;
    cov = Diagonal(([w0*rstep, w0*rstep, z0*sqrt(2)*rstep, vstep, vstep, vstep]) .^ 2);
    d = MvNormal(mean, cov);
    
    samples = [[0.0, 0.0, 0.0, vstep/sqrt(3), vstep/sqrt(3), vstep/sqrt(3)]];
    u_acc = rand(Uniform(0.0, 1.0), N*freq + skip);
    acc_rate = 0;
    
    for i ∈ 1:N*freq + skip - 1
        cord_last = samples[end];
        cord_new = cord_last + rand(d);
        p_acc = prob_boltzmann(cord_new, trap_params, atom_params; harmonic)/prob_boltzmann(cord_last, trap_params, atom_params; harmonic);
        
        if p_acc > u_acc[i] && H(cord_new, trap_params, m; harmonic) < U0
            push!(samples, cord_new);
            acc_rate += 1; 
        else
            push!(samples, cord_last);
        end;
    end;
        
    return samples[1+skip:freq:end], acc_rate/(N*freq + skip)
end;


#Generate coordinate trajectory from Monte-Carlo initial conditions
function R(t, ri, vi, ω)
    return ri * cos.(ω .* t) + vi/ω * sin.(ω .* t);
end;    


#Generate velocity trajectory from Monte-Carlo initial conditions
function V(t, ri, vi, ω)
    return vi * cos.(ω .* t) - ri * ω * sin.(ω .* t);
end;       