function [xi, zi, ui, dui, Tli]=MPC_No3(Ac,Bc,Cc,Dc,Ts,ref,Tl, NoS,umin,umax)
% Velocity form MPC
H = [1 0 0; 0 0 1];

sys = ss(Ac,Bc,Cc,Dc);

sys_d = c2d(sys, Ts);

nx = size(Ac, 1);

nu = size(Bc, 2);

ny = size(Cc, 1);


nz = 2;


N = 6;

Nu = N;

[Ad, Bd, Cd, Dd] = ssdata(sys_d);

% weighting matrices S, Q and R
Q = blkdiag(10, .001, 2)*10;
R = blkdiag(1, 1)*1;
[X,~,~,~] = dare(Ad,Bd, Q,R,[],[]);
S = X;

A_tilde = [Ad zeros(nx,ny);
    Cd eye(ny)];

B_tilde = [Bd;
    zeros(ny, nu)];

C_tilde = [Cd eye(ny)];

%construction of bbarQ, bbarT, bbarR, bbarH and bbarF
A_p = predA(C_tilde,A_tilde,N);
B_p = predB(C_tilde,A_tilde, B_tilde, N, Nu);
Q_p = predQ(Q,S,N);
R_p = kron(eye(Nu), R);
Hess = B_p'*Q_p*B_p + R_p;
Hessian = (Hess+Hess')/2;

% Memory Location definition
tildeX = zeros(nx+ny, NoS);
x = zeros(nx,NoS); %States initialised
dX = zeros(nx,NoS); %Change in states
u = zeros(nu,NoS); %Control input signal initialised
u_ = zeros(nu,1); %
dU = zeros(nu,NoS); %Change in control input signal which carries out integral function
y = zeros(ny, NoS); %Outputs initialisation

% Constraints implementation
Umax = kron(ones(Nu,1),umax);
Umin = kron(ones(Nu,1),umin);
Aqprog = [eye(Nu*nu); eye(Nu*nu)];

for k = 1:NoS-1

    y(:,k) = Cd*x(:,k);

    yref = kron(ones(N,1), [0; 0;ref(:,k)]);

    u_p = kron(ones(Nu,1), u_);

    bqprog = [Umax;-Umin]-[-u_p; u_p];

    f = B_p'*Q_p*(A_p*tildeX(:,k) - yref);

    V = quadprog(Hessian, f, Aqprog, bqprog);

    dU(:,k) = V(1:nu);

    u(:,k) = u_ + dU(:,k);

    x_next =  runge_kutta_pmsm(x(:,k), Ts, u(:,k), Tl(:,k));

    x(:,k+1) = awgn(x_next,30);

    dX(:,k+1) = x(:,k+1) - x(:,k);

    tildeX(:,k+1)=[dX(:,k+1) ; y(:,k)];

    u_ = u(:,k); %Updates the the previous control input

    sprintf('Iteration number = %d of %d',k, NoS)

end

xi = x;

zi = H*xi;

ui = u;

dui = dU;

Tli = Tl;

end