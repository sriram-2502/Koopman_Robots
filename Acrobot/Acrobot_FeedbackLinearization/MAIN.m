%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acrobot simulation with swing-up controls. The acrobot is a
% double pendulum with only a motor at the elbow (i.e. the second joint).
%
%   This is an implementation of "Partial Feedback Linearization of Underactuated
%   Mechanical Systems" by Mark Spong. I am not affiliated with him or this
%   paper.
%
% Files:
% MAIN - Execute this file; parameters here.
%
% Matthew Sheen, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear all;

%deriverRelativeAngles; % Should not need to run unless you want to redo
%all the equation derivations.

% Pick controller
p.controller = 'noncollocated'; % Choose: noncollocated, collocated, none
% 1) noncollocated controller is really crazy and can stabilize to any target
% angle! The downside is that it requires basically boundless torque.
% 2) collocated controller does a reasonable "pumping" motion for swing-up.
% The downside is that it can't stabilize without some linear controller to
% keep it at the top.
% 3) none - no controller mode. Play around with the free pendulum.
%
% NOTE: This code does not transition to LQR at the top.

% Initial conditions:
init = [0    0    0   0]';

% Simulation duration
duration = 20;
animationSpeed = 1;

% For link 1 linearization (noncollocated):
p.kd1 = 6.4;
p.kp1 = 15;
p.target = pi/2; % Target stabilization angle

% For link 2 linearization (collocated):
p.alpha = pi/6; % "pumping" angle
p.kd2 = 200;
p.kp2 = 2000;

%%%%%%%% System Parameters %%%%%%%%

p.g = 9.81; % Gravity
p.m1 = 1; % Mass of link 1.
p.m2 = 1; % Mass of link 2.
p.l1 = 1; % Total length of link 1.
p.l2 = 1; % Total length of link 2.
p.d1 = p.l1/2; % Center of mass distance along link 1 from the fixed joint.
p.d2 = p.l2/2; % Center of mass distance along link 2 from the fixed joint.
p.I1 = 0.083;% 1/12*p.m1*p.l1^2; %Moment of inertia of link 1 about COM
p.I2 = 1/3;% 1/12*p.m2*p.l2^2; %Moment of inertia of link 2 about COM

p.T1 = 0; % Motor torques default to zero. Not really relevant unless you run this without control (or maybe disturbances at first joint?)
p.T2 = 0; % Joint torque 2 only affects the system when controls aren't used.
p.sat = 10000; %Actuator saturation threshold (Nm). Obviously unrealistic for this kind of controller


%%%%%%%% Integrate %%%%%%%%

options1 = odeset('AbsTol', 1e-6,'RelTol',1e-6, 'Events', @LQR_switch); %Transition from swing up to linear balance controller when conditions met.
[tarray, zarray, te, ze, ie] = ode15s(@Dynamics, [0 duration], init, options1, p);

init_LQR = zarray(end,:);
[t_LQR, z_LQR] =  AcrobatLqr(p, init_LQR, te);

%Torques are post-calculated due to the difficulty of pulling numbers out
%of ODE45. Therefore, we also have to post-cast the values within the
%actuator limits.
if strcmp(p.controller, 'noncollocated')
    % Linearized link 1
    th1d = p.target;
    Tarray = ControlTorque1(p.I1,p.I2,p.d1,p.d2,p.g,p.kd1,p.kp1,p.l1,p.m1,p.m2,zarray(:,1),zarray(:,3),pi/2,zarray(:,2),zarray(:,4));
elseif strcmp(p.controller, 'collocated')
    % Linearized link 2
    th2d = p.alpha*atan(zarray(:,2));
    Tarray = ControlTorque2(p.I1,p.I2,p.d1,p.d2,p.g,8,16,p.l1,p.m1,p.m2,zarray(:,1),zarray(:,3),pi/2,zarray(:,2),zarray(:,4));
else % No controller
   Tarray = zeros(size(zarray,1),1); 
end

Tarray(Tarray>p.sat) = p.sat;
Tarray(Tarray<-p.sat) = -p.sat;

% Energy calculations
energy = TotEnergy(p.I1,p.I2,p.d1,p.d2,p.g,p.l1,p.m1,p.m2,zarray(:,1),zarray(:,3),zarray(:,2),zarray(:,4));

% Do animation
Plotter

% event function to stop swing up control and switch to LQR
function [value, isterminal, direction] = LQR_switch(tarray, zarray, p)
    value      = (zarray(3) < -19); %sweet spot for stopping when link 2 is close enough to eqb
    isterminal = 1;   % Stop the integration
    direction  = 0;
end


%%%%%%%% LQR stuff %%%%%%%%

%LQR Controller
function [t_LQR, z_LQR] = AcrobatLqr(p, init_LQR, te)
global m1 m2 l1 l2 lc1 lc2 I1 I2 g  A B Q R xd
m1 = p.m1 ; m2 = p.m2;
l1 = p.l1 ; l2 = p.l2; g = p.g;
lc1 = p.d1; lc2 = p.d2;
I1 = p.I1 ;I2 = p.I2; 
xd = [pi 0 0 0];

%x0 = [pi+0.1 0.003 0 0 ];

dtdq = [g*(m1*lc1+m2*l1+m2*lc2) m2*g*lc2;
       m2*g*lc2 m2*g*lc2];
   b = [0 ;1];
   M = [I1+I2+m2*l1^2+2*m2*l1*lc2  I2+m2*l1*lc2;
        I2+m2*l1*lc2  I2];
A    = [zeros(2) eye(2);
        M\dtdq zeros(2)];
B = [zeros(2,1);M\b];
Q = diag([10 10 1 1]);
R = 1;

options2 = odeset('AbsTol', 1e-6,'RelTol',1e-6);
[t_LQR,z_LQR] = ode45(@dynamics_LQR,[te duration], init_LQR, options2);
end

function xdot = dynamics_LQR(~,x)
global m1 m2 l1  lc1 lc2 I1 I2 g A B Q R xd
q1 = x(1); q2 = x(2); dq1 = x(3); dq2 = x(4);
m11 = I1+ I2 + m2*l1^2 + 2*m2*l1*lc2*cos(q2);
m12 = I2 + m2*l1*lc2*cos(q2);
m22 = I2;
M = [m11 m12 ; m12 m22];
C = [-2*m2*l1*lc2*sin(q2)*dq2 -m2*l1*lc2*sin(q2)*dq2;
     m2*l1*lc2*sin(q2)*dq1 0];
G = [m1*g*lc1*sin(q1) + m2*g*(l1*sin(q1)+lc2*sin(q1+q2));
     m2*g*lc2*sin(q1+q2)];
dq = [dq1;dq2]; B1 = [0;1];

K = lqr(A,B,Q,R);
u = -K*(x-xd');

ddq = M\(B1*u-C*dq -G);
xdot = [dq;ddq];
disp(xdot)
end


