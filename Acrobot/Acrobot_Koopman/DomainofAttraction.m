%% Control based eigenfunctions
% Eigenfunctions based on Least-Square
clear; close all; clc;
set(0,'DefaultLineLineWidth',2.3) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% Acrobot dynamics
% Automated Control Design for a Piecewise- Affine Approximation of a Class of Nonlinear Systems
% states are [q1 q2 q1dot q2dot]
x0 = [0; 0; 0; 0] ; %initial condition for new model with zero vertical
f = @(t,x)([ ...
    x(3,:); ...
    x(4,:); ...
    (2.*(cos(x(2,:)) + 2).*((981.*sin(pi-x(1,:) + x(2,:)))./200 - (-x(3,:).^2.*sin(x(2,:)))./2))./(cos(x(2,:)).^2 - 8) - (4.*((sin(x(2,:)).*-x(4,:).^2)./2 - x(3,:).*sin(x(2,:)).*x(4,:) + (981.*sin(pi- x(1,:) + x(2,:)))./200 + (2943.*sin(pi-x(1,:)))./200))./(cos(x(2,:)).^2 - 8);...
    (2.*(cos(x(2,:)) + 2).*((sin(x(2,:)).*x(4,:).^2)./2 -x(3,:).*sin(x(2,:)).*x(4,:) + (981.*sin(pi- x(1,:) + x(2,:)))./200 + (2943.*sin(pi-x(1,:)))./200))./(cos(x(2,:)).^2 - 8) - (4.*(cos(x(2,:)) + 3).*((981.*sin(pi - x(1,:) + x(2,:)))./200 - (-x(3,:).^2.*sin(x(2,:)))./2))./(cos(x(2,:)).^2 - 8);...
    ]);
g = @(t,x)[0; 0; -4./(cos(x(2,:)).^2 - 8); (2.*(cos(x(2,:)) + 2))./(cos(x(2,:)).^2 - 8)];
ft = @(t,x,u) f(t, x) + g(t, x)*u;
%% linearize the system
n=4; % dimension
x = sym('x',[n;1]);
A = double(subs(jacobian(f(0,x),x),x,x0)); 
u = sym('u',[1;1]);
B = double(subs(jacobian(ft(0,x,u),u),[x;u],[x0;0]));
C = eye(4,4);
D = 0;
sys = ss(A,B,C,D);
%% LQR Controller
Q = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %cost on states
R = 1; %cost on inputs
[K,S,E] = lqr(sys,Q,R);
uLQR = @(x) -K*(x);
f_LQR = @(t,x,u)([A*x+B*u]);
%% Simulation Parameters
% theta_initial = 45;
% dtheta_initial = 0;
% 
% %Run Simulation
% sim('part4sim.slx');

%% Basin of Attraction Analysis
numq1 = 40;
numq2 = 40;

q1_initials = linspace(0,180,numq1+1);
q2_initials = linspace(0,180,numq2+1);

counter = 0;
dataMatrix = [0 0 1];
for i=1:numq1
    for j=1:numq2
        counter = counter + 1;
        q1_initial = q1_initials(i); %[deg]
        q2_initial = q2_initials(j); %[deg/s]
        tic;
        %sim('part4sim.slx');
        options = odeset('RelTol',1e-3,'AbsTol',1e-10);
        tspan = 0:0.01:10;
        x0 = [q1_initial; q2_initial; 0; 0]; %for car pole
        [tspan2, xLQR] = ode45(@(t,x) f_LQR(t, x, uLQR(x)), tspan, x0, options);
        if (xLQR(end,1)<=1e-3) && (xLQR(end,2)<=1e-3)
            success = 1;
        else
            success = 0;
        end
        %dataMatrix = [dataMatrix; theta_initial dtheta_initial success];
        dataMatrix = [dataMatrix; q1_initial q2_initial success];
        disp('Iteration: ');
        disp(counter);
        disp('Percent Done:');
        disp(100*(counter/(numq1*numq2 + 2)));
        disp('Success: ');
        disp(success);
        disp('Sim Time');
        disp(toc);
    end
end

%dataMatrix = load('part4BasinOfAttraction.mat');
%dataMatrix = dataMatrix.dataMatrix;


%%Fill out results
% dataMatrixFull = [dataMatrix; ...
%     [-dataMatrix(:,1) dataMatrix(:,2) dataMatrix(:,3)];...
%     [dataMatrix(:,1) -dataMatrix(:,2) dataMatrix(:,3)];...
%     [-dataMatrix(:,1) -dataMatrix(:,2) dataMatrix(:,3)]];

%% Plot Results

dataMatrixFull = [dataMatrix; ...
    [-dataMatrix(:,1) dataMatrix(:,2) dataMatrix(:,3)];...
    [dataMatrix(:,1) -dataMatrix(:,2) dataMatrix(:,3)];...
    [-dataMatrix(:,1) -dataMatrix(:,2) dataMatrix(:,3)]];


figure(1);
hold on;
for i=1:length(dataMatrixFull)
    if (dataMatrixFull(i,3)==1)
        plot(dataMatrixFull(i,1), dataMatrixFull(i,2), 'g*');
    else
        plot(dataMatrixFull(i,1), dataMatrixFull(i,2), 'r*');
    end
end
hold off;
xlabel('q1 [deg]');
ylabel('q2 [deg]');
title('Basin of Attraction');
