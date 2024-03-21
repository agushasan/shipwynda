%% Research code by Agus Hasan
%% This code is to generate Figure 7 and Figure 8 (change the input below)
clear;
clc;

%% Time horizon
tf  = 20;
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

%% State initialization
x        = [3;3;0;0;0;0];
xD       = [3;3;0;0;0;0];

%% Initial control inputs
u     = [50 10 0]';        % Figure 7
%u     = [25 10 1]';        % Figure 8(a)
%u     = [5 10 0]';         % Figure 8(b)
%u     = [10 12 1]';        % Figure 8(c)
%u     = [10 10 0]';        % Figure 8(d)

%% For plotting
uArray          = [];
xArray          = [];
xDArray         = [];

%% Parameters from data-driven discovery
THETA = [-0.0229 -0.0365 -0.0304 -0.0540 0.0046 -0.0055 -0.0004 -0.0303 -0.0051 -0.0007 -0.0226 0.0128 0.0005 -1.0819 -0.0116 0.0001 -0.0040 0.0002 0.0265 -0.1236 -0.3624 -0.0157 0.3132 0.345 -0.0039 0.0213 -0.0039]';

%% simulation
for i=1:(tf/dt)

    uArray         = [uArray u];
    xArray         = [xArray x];
    xDArray        = [xDArray xD];    

    c13 = -(m-Yvd)*x(5)-(m*xg-Yrd)*x(6);
    c23 = (m-Xud)*x(4);
    Cv  = [0 0 c13; 0 0 c23;-c13 -c23 0];
    Dv  = [-Xu-Xuu*abs(x(4)) 0 0;0 -Yv-Yvv*abs(x(5)) -Yr;0 -Nv -Nr-Nrr*abs(x(6))];
    x = x+dt*[cos(x(3))*x(4)-sin(x(3))*x(5);sin(x(3))*x(4)+cos(x(3))*x(5);x(6);-inv(M)*(Cv+Dv)*[x(4);x(5);x(6)]]+B*dt*u;    

    c13D = -(m-Yvd)*xD(5)-(m*xg-Yrd)*xD(6);
    c23D = (m-Xud)*xD(4);
    CvD  = [0 0 c13D; 0 0 c23D;-c13D -c23D 0];
    DvD = [xD(4) xD(5) xD(6) xD(4)^2 xD(5)^2 xD(6)^2 xD(4)^3 xD(5)^3 xD(6)^3 zeros(1,18);
          zeros(1,9) xD(4) xD(5) xD(6) xD(4)^2 xD(5)^2 xD(6)^2 xD(4)^3 xD(5)^3 xD(6)^3 zeros(1,9);
          zeros(1,18) xD(4) xD(5) xD(6) xD(4)^2 xD(5)^2 xD(6)^2 xD(4)^3 xD(5)^3 xD(6)^3]*THETA;
    xD = xD+dt*[cos(xD(3))*xD(4)-sin(xD(3))*xD(5);sin(xD(3))*xD(4)+cos(xD(3))*xD(5);xD(6);-inv(M)*(CvD)*[xD(4);xD(5);xD(6)]+DvD]+B*dt*u;

end

figure(1)
plot(xArray(1,:),xArray(2,:), ':b', 'LineWidth', 6)
hold on;
plot(xDArray(1,:),xDArray(2,:), ':k', 'LineWidth', 6)
hold on;
plot(xArray(1,1),xArray(2,1), 'ok', 'LineWidth', 16)
hold on;
plot(xArray(1,end),xArray(2,end), 'or', 'LineWidth', 16)
hold on;
plot(xDArray(1,1),xDArray(2,1), 'ok', 'LineWidth', 16)
hold on;
plot(xDArray(1,end),xDArray(2,end), 'or', 'LineWidth', 16)
set(gca,'color','#BDEDFF','LineWidth',3,'FontSize',36)
grid on;
grid minor;
%legend('trajectory - first principle','trajectory - data-driven','start','end','','','FontSize',48)
xlabel('x (m)','FontSize',72)
ylabel('y (m)','FontSize',72)

rmse(xArray(1,:),xDArray(1,:))+rmse(xArray(2,:),xDArray(2,:))
