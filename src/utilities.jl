#Beam waist radius
function w(z, w0, z0)
    return w0 .* sqrt.(1.0 .+ (z ./z0) .^2)
end;



#Converter
function w0_to_z0(w0, λ)
    return π*w0^2/λ
end;



#Amplitude of gaussian beam with |E0|=1
function A(x, y, z, w0, z0)
    return (w0 ./ w(z, w0, z0)) .* exp.(- (x .^2 .+ y .^2) ./ (w(z, w0, z0) .^2))
end;



#Phase of gaussian beam
function A_phase(x, y, z, w0, z0)
    k = 2.0 * z0 / w0^2;
    return exp.(-1.0im * k * z .* (0.5*(x .^2 + y .^ 2) ./ (z .^2  + z0 .^2)) + 1.0im * atan.(z ./ z0));
end;



#Complex amplitude of gaussian beam with |E0|=1
function E(x, y, z, w0, z0)
    return A(x,y,z,w0,z0) .* Phase(x,y,z,w0,z0)
end;


#ωr, ωz trap frequencies
function trap_frequencies(atom_params, trap_params)
    m, T = atom_params;
    U0, w0, z0 = trap_params;
    ω = vconst/w0 * sqrt(U0/m);
    
    return 2*ω, sqrt(2)*ω
end;