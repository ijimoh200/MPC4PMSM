% Linearised system model
[Ac,Bc] = LinearisedPmsmDynamics();
Cc =eye(3);
Dc = 0;
% Specify sampling period (s)
Ts = 2e-3;
% Simulation duraction (s)
Tf = 2;
NoS = round(Tf/Ts);
% Define Reference signal
ref = zeros(1,NoS);
% Define unmeasured external torque disturbance
Tl  = zeros(1, NoS);
alpha  = zeros(1, NoS); % assume alpha1 = apha2
% Compute reference and disturbance signals for simulation
for k = 1:NoS
    if k <= round(NoS*.3)
        ref(:,k) = 15.7*2;
    elseif k <= round(NoS*.67)
        ref(:,k) = 15.7*3;
    else
        ref(:,k) = 15.7*2.5;
    end
    if k <= NoS*.5
        Tl(:,k) = 0;
        alpha(:,k) = 0;
    else
        Tl(:,k) = 0;
        alpha(:,k) = 0.6;
    end
end
dist = [Tl;
    alpha];

t = 0:Ts:Tf-Ts; % Simulation time

umax = [25.17; 51.96]; umin = -umax;

% Simulations

[x1, z1, xtildehat, u1, du1, Tl1] = MPC_No1(Ac,Bc,Cc,Dc,Ts,ref,dist, NoS,umin,umax);

[x2, z2, u2, du2, Tl2] = MPC_No2(Ac,Bc,Cc,Dc,Ts,ref,dist, NoS,umin,umax);

[x3, z3, u3, du3, Tl3] = MPC_No3(Ac,Bc,Cc,Dc,Ts,ref,dist, NoS,umin,umax);

yref1 = zeros(1,NoS);

yref2 = ref;

uhat = [0 0 0 1 0; 0 0 0 0 1]*xtildehat;
xhat = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0]*xtildehat;

% plot results
u1(:,end) = u1(:,end-1);
u2(:,end) = u2(:,end-1);


figure
subplot(311)
stairs(t, yref1,'-.k', 'linewidth', 1.5)
hold on
stairs(t, x1(1, 1:size(t,2)),'-r', 'linewidth', 1.2)
stairs(t, x2(1, 1:size(t,2)),'--b', 'linewidth', 1.2)
%stairs(t, x3(1, 1:size(t,2)),'-.g', 'linewidth', 1.2)
ylabel('$i_d$ [A]', 'interpreter','latex')
legend('Reference','MPC1','MPC2','MPC3', 'interpreter','latex')
xlabel('Time [s]')
ylim([-2 2.5])
grid on
grid minor

subplot(312)
stairs(t, x1(2, 1:size(t,2)),'-r', 'linewidth', 1.2)
hold on
stairs(t, x2(2, 1:size(t,2)),'--b', 'linewidth', 1.2)
%stairs(t, x3(2, 1:size(t,2)),'-.g', 'linewidth', 1.2)
ylabel('$i_q$ [A]', 'interpreter','latex')
xlabel('Time [s]')
ylim([-2 8])
grid on
grid minor

subplot(313)
stairs(t, 100/15.7*yref2,'-.k', 'linewidth', 1.2)
hold on
stairs(t, 100/15.7*x1(3, 1:size(t,2)),'r', 'linewidth', 1.2)
stairs(t, 100/15.7*x2(3, 1:size(t,2)),'--b', 'linewidth', 1.2)
%stairs(t, 100/15.7*x3(3, 1:size(t,2)),'-.g', 'linewidth', 1.2)
ylabel('$\omega_e$ [rpm]', 'interpreter','latex')
xlabel('Time [s]')
ylim([0 350])
grid on
grid minor


figure
subplot(211)
stairs(t, u1(1, 1:size(t,2)),'r', 'linewidth', 1.2)
hold on
stairs(t, u2(1, 1:size(t,2)),'--b', 'linewidth', 1.2)
%stairs(t, u3(1, 1:size(t,2)),'-.g', 'linewidth', 1.2)
stairs(t, umax(1)*ones(1,NoS),'-.k', 'linewidth', 1.2)
stairs(t, umin(1)*ones(1,NoS),'-.k', 'linewidth', 1.2)
ylabel('$V_d$ [V]', 'interpreter','latex')
legend('MPC1','MPC2','Limits', 'interpreter','latex')
xlabel('Time [s]')
ylim([-30 30])
grid on
grid minor

subplot(212)
stairs(t, u1(2, 1:size(t,2)),'r', 'linewidth', 1.2)
hold on
stairs(t, u2(2, 1:size(t,2)),'--b', 'linewidth', 1.2)
%stairs(t, u3(2, 1:size(t,2)),'-.g', 'linewidth', 1.2)
stairs(t, umax(2)*ones(1,NoS),'-.k', 'linewidth', 1.2)
stairs(t, umin(2)*ones(1,NoS),'-.k', 'linewidth', 1.2)
ylabel('$V_q$ [V]', 'interpreter','latex')
xlabel('Time [s]')
ylim([-60 60])
grid on
grid minor

figure
subplot(211)
plot(t,u1(1,:),'r', 'linewidth', 1.2)
hold on
plot(t,uhat(1,:),'--b', 'linewidth', 1.2)
legend('$V_d$','$\hat{V}_d$','interpreter','latex')
xlabel('Time [s]')
ylabel('$V_d$ [V]', 'interpreter','latex')
ylim([-30 30])
grid on
grid minor

subplot(212)
plot(t,u1(2,:),'r', 'linewidth', 1.2)
hold on
plot(t,uhat(2,:),'--b', 'linewidth', 1.2)
legend('$V_q$','$\hat{V}_q$','interpreter','latex')
ylabel('$V_q$ [V]', 'interpreter','latex')
xlabel('Time [s]')
ylim([-60 60])
grid on
grid minor



% Error vectors
error1 = z1 - [0*ref; ref];
error2 = z2 - [0*ref; ref];
error3 = z3 - [0*ref; ref];

% Preallocate result vectors
rmse1 = zeros(2,1); 
rmse2 = zeros(2,1); 
rmse3 = zeros(2,1); 
chdu1 = zeros(2,1); 
chdu2 = zeros(2,1); 
chdu3 = zeros(2,1); 

% Loop through each row (state variable)
for i = 1:2
    % Controller 1
    err1 = error1(i,50:end); % 0.1/Ts = 50
    rmse1(i,1) = sqrt(mean(err1.^2));

    % Controller 2
    err2 = error2(i,50:end);
    rmse2(i,1) = sqrt(mean(err2.^2));

    % Controller 3
    err3 = error3(i,50:end);
    rmse3(i,1) = sqrt(mean(err3.^2));
end

for i = 1:2
    % Controller 1
    err1 = du1(i,50:end); %
    chdu1(i,1) = sqrt(mean(err1.^2));
    
    % Controller 2
    err2 = du2(i,50:end);
    chdu2(i,1) = sqrt(mean(err2.^2));

    % Controller 3
    err3 = du3(i,50:end);
    chdu3(i,1) = sqrt(mean(err3.^2));
end

% Optional: State labels
labels = {'i_d','w_e'};
% Display results
fprintf('\n--- RMSE Comparison ---\n');
for i = 1:2
    fprintf('%s: Ctrl1 = %.4f, Ctrl2 = %.4f, Ctrl3 = %.4f\n', ...
        labels{i}, rmse1(i), rmse2(i), rmse3(i));
end

labels = {'V_d','V_q'};
fprintf('\n--- RMS du Comparison ---\n');
for i = 1:2
    fprintf('%s: Ctrl1 = %.4f, Ctrl2 = %.4f, Ctrl3 = %.4f\n', ...
        labels{i}, chdu1(i), chdu2(i), chdu3(i));
end