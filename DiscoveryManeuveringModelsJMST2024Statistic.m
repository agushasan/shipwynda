%% Research code by Agus Hasan
%% This code is to generate Figure 12
clear;
clc;

%% Time horizon
tf  = 5;
dt  = 0.001;
t   = dt:dt:tf;

%% Define the number of variables
n = 3;      % number of measured variables
r = 10;     % number of parameters
l = 10;
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
Cb = eye(3);

%% State initialization
x        = [0;0;0;0;0;0];
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
u     = [100 0 1]';

%% For plotting
uArray              = [];
xArray              = [];
yArray              = [];
vbarArray           = [];
thetabarArray       = [];
ABSErrorArray       = [];
ABSErrorArrayTotal  = [];

%% Initialization for estimator

lambdav = 0.999;
lambdat = 0.999;
Rv      = 0.001*eye(n);
Rt      = 0.001*eye(n);
Pv      = 1*eye(n);
Pt      = 1*eye(r);
Gamma   = zeros(n,r);

for k=1:3

    if k==1
        Ra = -1+2*rand;
        Rb = 0;
        Rc = 0;
    end
    if k==2
        Ra = 0;
        Rb = -1+2*rand;
        Rc = 0;
    end
    if k==3
        Ra = 0;
        Rb = 0;
        Rc = -1+2*rand;
    end

    ABSErrorArray = [];

for j = 1:l

uArray          = [];
xArray          = [];
yArray          = [];
vbarArray       = [];
thetabarArray   = [];

%% Simulation
for i=1:(tf/dt)
    
    u     = [(10+10*Ra)*sin(i*dt) (10+10*Rb)*cos(i*dt) (10+10*Rc)*cos(i*dt)]';

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

    Phi = [y(1) abs(y(1))*y(1) 0 0 0 0 0 0 0 0;
          0 0 y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3) 0 0 0 0;
          0 0 0 0 0 0 y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3)];

    bk  = dt*inv(M)*(u-Cv*y);
    
    % Estimation using adaptive observer
    Kv = Pv*Cb'*inv(Cb*Pv*Cb'+Rv);
    Kt = Pt*Gamma'*Cb'*inv(Cb*Gamma*Pt*Gamma'*Cb'+Rt);
    Gamma = (eye(n)-Kv*Cb)*Gamma;

    vbar = vbar+(Kv+Gamma*Kt)*(y-Cb*vbar);
    thetabar = thetabar-Kt*(y-Cb*vbar);

    vbar = eye(n)*vbar+bk+Phi*thetabar;
    thetabar = thetabar;
    Pv = (1/lambdav)*eye(n)*(eye(n)-Kv*Cb)*Pv*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*Cb*Gamma)*Pt;
    Gamma = eye(rank(Cb))*Gamma-Phi;

end

Temp1 = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetabarArray(3,:);thetabarArray(7,:)];
Temp2 = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetabarArray(4,:);thetabarArray(8,:)];

ABSError = sqrt(abs(Xu-m11*thetabarArray(1,end)/dt)^2+abs(Xuu-m11*thetabarArray(2,end)/dt)^2+abs(Nv-Temp1(1,end))^2+abs(Yv-Temp1(2,end))^2+abs(Nr-Temp2(1,end))^2+abs(Yr-Temp2(2,end))^2+abs(Yvv-mt*thetabarArray(5,end)/(dt*m33))^2+abs(Nrr+mt*thetabarArray(6,end)/(dt*m23)))^2;

ABSErrorArray = [ABSErrorArray ABSError];

end

ABSErrorArrayTotal = [ABSErrorArrayTotal ABSErrorArray];

end

figure(1)
hAx=gca;
boxplot([ABSErrorArrayTotal(1:l)',ABSErrorArrayTotal(l+1:2*l)',ABSErrorArrayTotal(2*l+1:3*l)'],'Labels',{'\tau_u','\tau_v','\tau_r'})
hAx.YAxis.Scale ="log";
hAx.XAxis.TickLabelInterpreter='tex';
set(gca,'color','#BDEDFF','FontSize',48)
set(gca,'linew',3)
bx = findobj('Tag','boxplot');
set(bx.Children,'LineWidth',3)
ylabel('log(RMSE)')
grid on;
grid minor;