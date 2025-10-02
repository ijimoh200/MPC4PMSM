function [Ac,Bc] = LinearisedPmsmDynamics()
% Load parameters
params = motorParams();

% Define symbolic variables
syms id iq we vd vq Tl

% Electromagnetic torque
Te = (3/2) * (params.p * params.phi * iq);

% State equations
X1 = (1 / params.Ld) * (vd - params.R * id + we * params.Lq * iq);
X2 = (1 / params.Lq) * (vq - params.R * iq - we * params.Ld * id - we * params.phi);
X3 = (params.p / params.J) * (Te - params.Bv / params.p * we - Tl);

% Compute Jacobian matrices
Ac = jacobian([X1 X2 X3]', [id iq we]);   % w.r.t. states
Bc = jacobian([X1 X2 X3]', [vd vq]);      % w.r.t. inputs

% Define operating point (example values)
id0 = 0;     % d-axis current [A]
iq0 = 0;     % q-axis current [A]
we0 = 15.7;  % Electrical speed [rad/s]

% Evaluate Ac and Bc at x0 = [id0, iq0, we0]
Ac_eval = subs(Ac, {'id', 'iq', 'we'}, {id0, iq0, we0});
Ac = double(Ac_eval);  % Convert symbolic to numeric
% Evaluate the system matrices
Bc = eval(Bc);
end

  

   
