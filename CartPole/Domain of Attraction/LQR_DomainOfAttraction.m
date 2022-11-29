%% Domain of Attraction for cart pole
% Eigenfunctions based on Least-Square
clear; close all; clc;
set(0,'DefaultLineLineWidth',2.3) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% Cart Pole dynamics
% Automated Control Design for a Piecewise- Affine Approximation of a Class of Nonlinear Systems
% states are [x xdot theta thetadot]
m = 2; M = 8;
l = 0.5; gr = 9.8; 

A = [0 1 0 0; 0 0 m*gr/M 0;0 0 0 1;0 0 (m+M)*gr/(M*l) 0];
B = [0; 1/M; 0; 1/(M*l)];
C = eye(4);
D = 0;
sys = ss(A,B,C,D);

f = @PendCartModel_F;
%% LQR Controller
Q = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %cost on states
R = 1; %cost on inputs
[K,S,E] = lqr(sys,Q,R);
uLQR = @(x) -K*(x);
%f_LQR = @(t,x,u)([A*x+B*u]);
%% Simulation Parameters
theta_initial = 45;
dtheta_initial = 0;

%Run Simulation
%sim('part4sim.slx');
%% Basin of Attraction Analysis 
numThetas = 25;
numDThetas = 25;
numx=25;
theta_initials = linspace(0,(pi/2)-0.1,numThetas+1);
dtheta_initials = linspace(0,deg2rad(100),numDThetas+1);
x_initials = linspace(0,5,numx+1);

%% x vs theta
counter = 0;
dataMatrix2 = [0 0 1];
for i=1:numx
    for j=1:numThetas
        counter = counter + 1;
        theta_initial = theta_initials(i); %[deg]
        x_initial = x_initials(j); 
        tic;
        %sim('part4sim.slx');
        options = odeset('RelTol',1e-10,'AbsTol',1e-20);
        tspan = 0:0.01:10;
        x02 = [x_initial;0;theta_initial;0]; % for theta vs dtheta
        [tspan2, xLQR2] = ode45(@(t,x) f(t, x, uLQR(x)), tspan, x02, options);
        %if (xLQR2(end,1)<=1e-5) && (xLQR2(end,3)<=1e-5) %check for x and theta
        if(abs(xLQR2(end,3))<=1e-5)
            success = 1;
        else
            success = 0;
        end
        %x vs theta
        dataMatrix2 = [dataMatrix2; x_initial theta_initial success];
        
        disp('Iteration: ');
        disp(counter);
        disp('Percent Done:');
        disp(100*(counter/(numThetas*numx + 2)));
        disp('Success: ');
        disp(success);
        disp('Sim Time');
        disp(toc);
    end
end

%% Plot Results

dataMatrixFull2 = [dataMatrix2; ...
    [-dataMatrix2(:,1) dataMatrix2(:,2) dataMatrix2(:,3)];...
    [dataMatrix2(:,1) -dataMatrix2(:,2) dataMatrix2(:,3)];...
    [-dataMatrix2(:,1) -dataMatrix2(:,2) dataMatrix2(:,3)]];

figure(1);
% for x vs theta
hold on;
for i=1:length(dataMatrixFull2)
    if (dataMatrixFull2(i,3)==1)
        plot(dataMatrixFull2(i,1), dataMatrixFull2(i,2), 'g*');
    else
        plot(dataMatrixFull2(i,1), dataMatrixFull2(i,2), 'r*');
    end
end
hold off;
xlabel('x [m]');
ylabel('Theta [rad]');
title('Basin of Attraction');


%% theta vs dtheta
counter = 0;
dataMatrix1 = [0 0 1];
for i=1:numThetas
    for j=1:numDThetas
        counter = counter + 1;
        theta_initial = theta_initials(i); %[deg]
        dtheta_initial = dtheta_initials(j); %[deg/s]
        tic;
        %sim('part4sim.slx');
        options = odeset('RelTol',1e-10,'AbsTol',1e-20);
        tspan = 0:0.01:10;
        x01 = [0;0;theta_initial;dtheta_initial]; % for theta vs dtheta
        [tspan1, xLQR1] = ode45(@(t,x) f(t, x, uLQR(x)), tspan, x01, options);
        %if (xLQR2(end,1)<=1e-5) && (xLQR2(end,3)<=1e-5) %check for x and theta
        %if (xLQR1(end,3)<=1e-3) && (xLQR1(end,4)<=1e-3) %check for theta and dtheta
        if(abs(xLQR1(end,3))<=1e-5)
            success = 1;
        else
            success = 0;
        end
        %theta vs dtheta
        dataMatrix1 = [dataMatrix1; theta_initial dtheta_initial success];
               
        disp('Iteration: ');
        disp(counter);
        disp('Percent Done:');
        disp(100*(counter/(numThetas*numx + 2)));
        disp('Success: ');
        disp(success);
        disp('Sim Time');
        disp(toc);
    end
end
%% Plot Results

dataMatrixFull1 = [dataMatrix1; ...
    [-dataMatrix1(:,1) dataMatrix1(:,2) dataMatrix1(:,3)];...
    [dataMatrix1(:,1) -dataMatrix1(:,2) dataMatrix1(:,3)];...
    [-dataMatrix1(:,1) -dataMatrix1(:,2) dataMatrix1(:,3)]];

figure(2);
%for theta vs dtheta
hold on;
for i=1:length(dataMatrixFull1)
    if (dataMatrixFull1(i,3)==1)
        plot(dataMatrixFull1(i,1), dataMatrixFull1(i,2), 'g*');
    else
        plot(dataMatrixFull1(i,1), dataMatrixFull1(i,2), 'r*');
    end
end
hold off;
xlabel('Theta [rad]');
ylabel('dTheta [rad/s]');
title('Basin of Attraction');

