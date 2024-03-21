%% Research code by Agus Hasan
%% This code is to generate Figure 9, Figure 10, and Figure 11
clear;
clc;

%% Time horizon
tf  = 10;
dt  = 0.001;
t   = dt:dt:tf;

%% Define the number of variables
n = 3;      % number of measured variables
r = 10;     % number of parameters

%% True parameters
% parameters for the inertia matrix
m    = 23.8;
Iz   = 1.76;
xg   = 0.046;
Xud  = -2;
Yvd  = -10;
Yrd  = 0;
Nvd  = 0;
Nrd  = -1;
% true hydrodynamic coefficients
Xu   = -0.7225;
Xuu  = -1.3274;
Yv   = -0.8612;
Yvv  = -36.2823;
Yr   = 0.1079;
Nv   = 0.1052;
Nr  = -0.5;
Nrr = -1;

%% System description
M = [m-Xud 0 0;0 m-Yvd m*xg-Yrd;0 m*xg-Nvd Iz-Nrd];
B = [zeros(n);inv(M)];
C = [0 0 0 1 0 0;0 0 0 0 1 0; 0 0 0 0 0 1];

%% Noise
RF = 0.3;

%% State initialization
x        = [3;3;0;0;0;0];
y        = [0;0;0];
vbar     = [0;0;0];
thetabar = zeros(r,1);
 
%% Known paramaters
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);
m33 = M(3,3);
mt  = m22*m33-m23*m32;

%% Initial control inputs
u     = [50 10 0]';

%% For plotting
uArray          = [];
xArray          = [];
yArray          = [];
vbarArray       = [];
thetabarArray   = [];

%% Initialization for estimator

lambdav = 0.99;
lambdat = 0.9999;
Rv      = 0.001*eye(n);
Rt      = 0.001*eye(n);
Pv      = 0.001*eye(n);
Pt      = 0.001*eye(r);
Gamma   = zeros(n,r);

%% simulation
for i=1:(tf/dt)
    %u     = [50 10 0]';                                     % Figure 9 and Figure 10 
    u     = [100*sin(i*dt) 10*cos(i*dt) 1*cos(i*dt)]';     % Figure 11

    uArray         = [uArray u];
    xArray         = [xArray x];
    yArray         = [yArray y];
    vbarArray      = [vbarArray vbar];
    thetabarArray  = [thetabarArray thetabar]; 

    c13 = -(m-Yvd)*x(5)-(m*xg-Yrd)*x(6);
    c23 = (m-Xud)*x(4);
    Cv  = [0 0 c13; 0 0 c23;-c13 -c23 0];
    Dv  = [-Xu-Xuu*abs(x(4)) 0 0;0 -Yv-Yvv*abs(x(5)) -Yr;0 -Nv -Nr-Nrr*abs(x(6))];

    x = x+dt*[cos(x(3))*x(4)-sin(x(3))*x(5);sin(x(3))*x(4)+cos(x(3))*x(5);x(6);-inv(M)*(Cv+Dv)*[x(4);x(5);x(6)]]+B*dt*u;
    y = C*x+RF^2*randn(n,1);

   Phi = [y(1) abs(y(1))*y(1) 0 0 0 0 0 0 0 0;
          0 0 y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3) 0 0 0 0;
          0 0 0 0 0 0 y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3)];

    bk  = dt*inv(M)*(u-Cv*y);
    
    % Estimation using adaptive observer
    Kv = Pv*inv(Pv+Rv);
    Kt = Pt*Gamma'*inv(Gamma*Pt*Gamma'+Rt);
    Gamma = (eye(n)-Kv)*Gamma;

    vbar = vbar+(Kv+Gamma*Kt)*(y-vbar);
    thetabar = thetabar-Kt*(y-vbar);

    vbar = eye(n)*vbar+bk+Phi*thetabar;
    thetabar = thetabar;
    Pv = (1/lambdav)*eye(n)*(eye(n)-Kv)*Pv*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*Gamma)*Pt;
    Gamma = eye(n)*Gamma-Phi;

end

Temp1 = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetabarArray(3,:);thetabarArray(7,:)];
Temp2 = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetabarArray(4,:);thetabarArray(8,:)];

figure(1)
subplot(3,1,1)
plot(t,yArray(1,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vbarArray(1,:), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('u [m/s]','FontSize',48)
subplot(3,1,2)
plot(t,yArray(2,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vbarArray(2,:), ':m', 'LineWidth', 6)
grid on;
grid minor;
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
ylabel('v [m/s]','FontSize',48)
legend('measured','estimated','FontSize',48)
subplot(3,1,3)
plot(t,yArray(3,:), '-k', 'LineWidth', 6)
hold on;
plot(t,vbarArray(3,:), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('r [rad/s]','FontSize',48)
xlabel('time (s)','FontSize',48)

figure(2)
subplot(4,2,1)
plot(t,Xu*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,m11*thetabarArray(1,:)/dt, ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Xu-4 Xu+4]);
ylabel('X_u','FontSize',48)
subplot(4,2,2)
plot(t,Xuu*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,m11*thetabarArray(2,:)/dt, ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('X_{uu}','FontSize',48)
legend('true parameter','estimated parameter','FontSize',36)
ylim([Xuu-4 Xuu+4]);
subplot(4,2,3)
plot(t,Nv*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,Temp1(1,:), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Nv-4 Nv+4]);
ylabel('N_v','FontSize',48)
subplot(4,2,4)
plot(t,Yv*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,Temp1(2,:), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('Y_v','FontSize',48)
ylim([Yv-4 Yv+4]);
subplot(4,2,5)
plot(t,Nr*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,Temp2(1,:), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Nr-4 Nr+4]);
ylabel('N_r','FontSize',48)
subplot(4,2,6)
plot(t,Yr*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,Temp2(2,:), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('Y_r','FontSize',48)
ylim([Yr-4 Yr+4]);
subplot(4,2,7)
plot(t,Yvv*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,mt*thetabarArray(5,:)/(dt*m33), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlabel('time (s)','FontSize',48)
ylim([Yvv-4 Yvv+4]);
ylabel('Y_{vv}','FontSize',48)
subplot(4,2,8)
plot(t,Nrr*ones(length(t),1), '-k', 'LineWidth', 6)
hold on;
plot(t,-mt*thetabarArray(6,:)/(dt*m23), ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('N_{rr}','FontSize',48)
xlabel('time (s)','FontSize',48)
ylim([Nrr-4 Nrr+4]);