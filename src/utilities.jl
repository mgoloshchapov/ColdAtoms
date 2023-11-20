#Beam waist radius
function w(z, w0, z0)
    return w0 .* sqrt.(1.0 .+ (z ./z0) .^2);
end;



#Converter
function w0_to_z0(w0, λ)
    return π*w0^2/λ;
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


"""
atom_params: [m(a.u.), T(μK)]
trap_params: [U0(μK), w0(μm), z0(μm)]
"""
function trap_frequencies(atom_params, trap_params)
    m, T = atom_params;
    U0, w0, z0 = trap_params;
    ω = vconst/w0 * sqrt(U0/m);
    
    return 2*ω, sqrt(2)*ω
end;


#Function for visualisation of samples
function samples_visualise(samples)
    x, y, z, vx, vy, vz = invert(samples);
    figure(figsize=(6,6))
    subplot(221)
    hist2D(x, z, bins=50, range=[[-1.0, 1.0], [-2.5, 2.5]],cmap="plasma", rasterized=true);
    xlabel("x, μm")
    ylabel("z, μm")
    title("Coordinate distribution");

    subplot(222)
    hist(x, bins=[-1.5:0.02:1.5;], density=true, alpha=0.5, label="x")
    hist(z, bins=[-2.0:0.1:2.0;], density=true, alpha=0.5, label="z")
    xlabel("μm")
    ylabel("pdf")
    title("Coordinate distribution");
    legend()

    subplot(223)
    hist2D(vx, vz, bins=50, range=[[-0.3, 0.3], [-0.3, 0.3]],cmap="plasma", rasterized=true);
    xlabel("\$ v_x \$, \$ \\mu m/ \\mu s \$")
    ylabel("\$ v_z \$, \$ \\mu m/ \\mu s \$")
    title("Velocity distribution");


    subplot(224)
    hist(vx, bins=[-0.3:0.01:0.3;], density=true, alpha=0.5, label="\$ v_x \$")
    hist(vz, bins=[-0.3:0.01:0.3;], density=true, alpha=0.5, label="\$ v_z \$")
    xlabel("\$ \\mu m/ \\mu s \$")
    ylabel("pdf")
    title("Velocity distribution");
    legend()

    tight_layout()
end;