%% Research code by Agus Hasan
%% This code is to generate Figure 6
clear;
clc;

%% Time horizon
tf  = 16;
dt  = 0.001;
t   = dt:dt:tf;

%% Define the number of variables
n = 3;      % number of measured variables
r = 27;     % number of parameters

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
B = [zeros(3);inv(M)];
C = [0 0 0 1 0 0;0 0 0 0 1 0; 0 0 0 0 0 1];

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

%% Simulation
for i=1:(tf/dt)
    
%    u     = [50*sin(i*dt) 10*cos(i*dt) 1*cos(i*dt)]';  

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
    y = C*x;

    Phi = [y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)^3 y(2)^3 y(3)^3 zeros(1,18);
          zeros(1,9) y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)^3 y(2)^3 y(3)^3 zeros(1,9);
          zeros(1,18) y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)^3 y(2)^3 y(3)^3];

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
subplot(9,3,1)
plot(t,thetabarArray(1,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_1','FontSize',24)
subplot(9,3,2)
plot(t,thetabarArray(2,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_2','FontSize',24)
subplot(9,3,3)
plot(t,thetabarArray(3,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_3','FontSize',24)
subplot(9,3,4)
plot(t,thetabarArray(4,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_4','FontSize',24)
subplot(9,3,5)
plot(t,thetabarArray(5,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_5','FontSize',24)
subplot(9,3,6)
plot(t,thetabarArray(6,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_6','FontSize',24)
subplot(9,3,7)
plot(t,thetabarArray(7,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_7','FontSize',24)
subplot(9,3,8)
plot(t,thetabarArray(8,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_8','FontSize',24)
subplot(9,3,9)
plot(t,thetabarArray(9,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_9','FontSize',24)
subplot(9,3,10)
plot(t,thetabarArray(10,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{10}','FontSize',24)
subplot(9,3,11)
plot(t,thetabarArray(11,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{11}','FontSize',24)
subplot(9,3,12)
plot(t,thetabarArray(12,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{12}','FontSize',24)
subplot(9,3,13)
plot(t,thetabarArray(13,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{13}','FontSize',24)
subplot(9,3,14)
plot(t,thetabarArray(14,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{14}','FontSize',24)
subplot(9,3,15)
plot(t,thetabarArray(15,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{15}','FontSize',24)
subplot(9,3,16)
plot(t,thetabarArray(16,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{16}','FontSize',24)
subplot(9,3,17)
plot(t,thetabarArray(17,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{17}','FontSize',24)
subplot(9,3,18)
plot(t,thetabarArray(18,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{18}','FontSize',24)
subplot(9,3,19)
plot(t,thetabarArray(19,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{19}','FontSize',24)
subplot(9,3,20)
plot(t,thetabarArray(20,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{20}','FontSize',24)
subplot(9,3,21)
plot(t,thetabarArray(21,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{21}','FontSize',24)
subplot(9,3,22)
plot(t,thetabarArray(22,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{22}','FontSize',24)
subplot(9,3,23)
plot(t,thetabarArray(23,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{23}','FontSize',24)
subplot(9,3,24)
plot(t,thetabarArray(24,:)/dt, ':m', 'LineWidth', 6)
set(gca,'xticklabel',[],'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{24}','FontSize',24)
subplot(9,3,25)
plot(t,thetabarArray(25,:)/dt, ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{25}','FontSize',24)
xlabel('time (s)','FontSize',48)
subplot(9,3,26)
plot(t,thetabarArray(26,:)/dt, ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{26}','FontSize',24)
xlabel('time (s)','FontSize',48)
subplot(9,3,27)
plot(t,thetabarArray(27,:)/dt, ':m', 'LineWidth', 6)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\theta_{27}','FontSize',24)
xlabel('time (s)','FontSize',48)

Coeff = [round(thetabarArray(1:r/n,end)/dt,4) round(thetabarArray((r/n)+1:2*r/n,end)/dt,4) round(thetabarArray((2*r/n)+1:r,end)/dt,4)];
round(thetabarArray(:,end)/dt,4)