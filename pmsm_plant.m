function X = pmsm_plant(t,x,u,dist)
%where u is the applied control input 
% and x represent the initial states
% and Tl is the torque disturbance

% Load model parameters
params = motorParams();

% Define state and input variables
id = x(1);
iq = x(2);
we = x(3);
vd = u(1);
vq = u(2);

% Electromagnetic torque
Te = (3/2) * (params.p * params.phi * iq);
%Extract disturbances
Tl = dist(1,:);
alpha = dist(2,:);
% State equations
X(1,1) = (1 / params.Ld) * ((1-alpha)*vd - params.R * id + we * params.Lq * iq);
X(2,1) = (1 / params.Lq) * ((1-alpha)*vq - params.R * iq - we * params.Ld * id - we * params.phi);
X(3,1) = (params.p / params.J) * (Te - params.Bv / params.p * we - Tl);

end