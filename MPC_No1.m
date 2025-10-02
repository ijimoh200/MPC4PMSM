function [xi, zi, xtildei, ui, dui, Tli] = MPC_No1(Ac,Bc,Cc,Dc,Ts,ref,dist, NoS,umin,umax)
% Increment form without model augmentation but with state and input estimation
H = [1 0 0; 0 0 1];

sys = ss(Ac,Bc,Cc,Dc);

sys_d = c2d(sys, Ts);

n = size(Ac, 1);

m = size(Bc, 2);

p = size(Cc, 1);

nz = 2;


N = 6;

Nu = N;

[Ad, Bd, Cd, Dd] = ssdata(sys_d);

Atilde = [Ad Bd; zeros(m,n) eye(m)];

Btilde = [Bd; eye(m)];


% weighting matrices S, Q and R
Q = blkdiag(10, .001, .1)*10;
R = blkdiag(1, 1)*1;
[X,~,~,~] = dare(Ad,Bd, Q,R,[],[]);
S = X;

% Contruction of prediction model matrices
A_p = predA(Cd,Ad,N);
B_p = predB(Cd,Ad, Bd, N, Nu);
Q_p = predQ(Q,S,N);
R_p =  kron(eye(Nu), R);
Hess = B_p'*Q_p*B_p + R_p;
Hessian = (Hess+Hess')/2;

% Constraints implementation
Umax = kron(ones(Nu,1),umax);
Umin = kron(ones(Nu,1),umin);
G = buildOmega1(m, N);
Aqprog = [-G; G];

% Memory Locations
xtilde       = zeros(n+m,NoS);
xhat       = zeros(n,NoS);
x       = zeros(n,NoS);
y       = zeros(p,NoS);
z       = zeros(nz,NoS);    % Outputs
u       = zeros(m,NoS);
du      = zeros(m,NoS);
u_      = zeros(m, 1);     % ui(k-1) i.e., previous computed velocity
uhat       = zeros(m,NoS);
Ctilde_obs = [H zeros(nz,m)];

%Kalman Filtering
Qw = blkdiag(1^2,1^2,1^2,1^2,1^2)*1; % Proces Noise Covariance Matrix w
Rv = blkdiag(1^2,1^2)*1; % Measurement Noise Covariance Matrix v
[~,K,~] = idare(Atilde',Ctilde_obs',Qw,Rv,zeros(5,2),eye(5));

L = K';

% Simulation Loop

for k = 1:NoS-1

    z(:,k) = H*x(:,k);

    xref = kron(ones(N,1), [0; 0;ref(:,k)]);

    uhat_p = kron(ones(Nu,1), uhat(:,k));

    f = B_p'*Q_p*(A_p*xhat(:,k) - xref + B_p*uhat_p);

    bin = [Umax;-Umin];

    Bin = [kron(ones(Nu,1), eye(m)); -kron(ones(Nu,1), eye(m))];

    bqprog = bin-Bin*u_;

    options = optimset('Display', 'off');

    V = quadprog(Hessian,f, Aqprog, bqprog,[],[],[],[],[],options);

    u(:,k) = V(1:m);

    du(:,k) = V(1:m);

    u(:,k) = u_ + du(:,k);

    x_next = runge_kutta_pmsm(x(:,k), Ts, u(:,k), dist(:,k));

    x(:,k+1) = awgn(x_next,30);

    xtilde(:,k+1) = Atilde*xtilde(:,k) + Btilde*du(:,k) + L*(z(:,k)-Ctilde_obs*xtilde(:,k));

    xhat(:,k+1) = xtilde(1:3,k+1);

    uhat(:,k+1) = xtilde(4:5,k+1);

    u_ = u(:,k) ;

    sprintf('Iteration number = %d of %d',k, NoS)

end

xi = x;

zi = H*xi;

ui = u;

dui = du;

xtildei = xtilde;

Tli = dist;

end





