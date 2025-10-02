function [xi, zi, ui, dui, Tli]=MPC_No2(Ac,Bc,Cc,Dc,Ts,ref,dist, NoS,umin,umax)
% Delta input form with no input estimation

H = [1 0 0; 0 0 1];

sys = ss(Ac,Bc,Cc,Dc);

sys_d = c2d(sys, Ts);

n = size(Ac, 1);

m = size(Bc, 2);

p = size(Cc, 1);

N = 6;

Nu = N;

[Ad, Bd, Cd, Dd] = ssdata(sys_d);

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
x       = zeros(n,NoS);
u       = zeros(m,NoS);
du      = zeros(m,NoS);
u_      = zeros(m, 1);     % ui(k-1) i.e., previous computed velocity

% Simulation Loop

for k = 1:NoS-1

    yref = kron(ones(N,1), [0; 0;ref(:,k)]);

    u_p = kron(ones(Nu,1), u_);

    f = B_p'*Q_p*(A_p*x(:,k) - yref + B_p*u_p);

    bin = [Umax;-Umin];

    Bin = [kron(ones(Nu,1), eye(m)); -kron(ones(Nu,1), eye(m))];

    bqprog = bin-Bin*u_;

    options = optimset('Display', 'off');

    V = quadprog(Hessian,f, Aqprog, bqprog,[],[],[],[],[],options);

    du(:,k) = V(1:m);

    u(:,k) = u_ + du(:,k);

    x_next = runge_kutta_pmsm(x(:,k), Ts, u(:,k), dist(:,k));

    x(:,k+1) = awgn(x_next,30);

    u_ = u(:,k);

    sprintf('Iteration number = %d of %d',k, NoS)

end

xi = x;

zi = H*xi;

ui = u;

dui = du;

Tli = dist;
