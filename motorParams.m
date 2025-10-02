function params = motorParams()
    % Physical and electrical parameters for PMSM
    params.p    = 2;              % Number of pole pairs
    params.J    = 2.35*10e-4;        % Rotor inertia [kg·m²]
    params.Bv   = 1.1*10e-4;         % Viscous damping [N·m·s]
    params.Lq   = 7*10e-3;           % q-axis inductance [H]
    params.Ld   = params.Lq;      % d-axis inductance [H] (assuming equal)
    params.R    = 2.98;           % Stator resistance [Ω]
    params.phi  = 0.125;          % Permanent magnet flux linkage [Wb]
end
