function AcrobatLqr()
clc;clear all; close all;
global m1 m2 l1 l2 lc1 lc2 I1 I2 g  A B Q R xd
m1 = 1 ; m2 = 1;
l1 = 1 ; l2 = 1;g = 9.81;
lc1 = l1/2; lc2 = l2/2;
I1 = 1 ;I2 = 1; 
xd = [pi 0 0 0];
x0 = [pi+0.1 0.003 0 0 ];
tspan = 0:0.01:5;
dtdq = [g*(m1*lc1+m2*l1+m2*lc2) m2*g*lc2;
       m2*g*lc2 m2*g*lc2];
   b = [0 ;1];
   M = [I1+I2+m2*l1^2+2*m2*l1*lc2  I2+m2*l1*lc2;
        I2+m2*l1*lc2  I2];
A    = [zeros(2) eye(2);
        M\dtdq zeros(2)];
eig(A)
B = [zeros(2,1);M\b];
Q = diag([10 10 1 1]);
R = 1;
[t,x] = ode45(@dynamics,tspan,x0);
Animation(t,x)
end

function xdot = dynamics(~,x)
global m1 m2 l1  lc1 lc2 I1 I2 g 
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
u = control(x);
ddq = M\(B1*u-C*dq -G);
xdot = [dq;ddq];
end

function u = control(x)
global A B Q R xd
K = lqr(A,B,Q,R);
u = -K*(x-xd');
end

function Animation(t,x)
global  l1 l2 
N = length(t);
th1 = x(1,1);th2 = x(1,2);
x1 = l1*sin(th1);y1 = -l1*cos(th1);
x2 = x1 + l2*sin(th1+th2);y2 = y1 - l2*cos(th1+th2);
link1 = line([0 x1],[0 y1],'color','k','LineWidth',4);
link2 = line([x1 x2],[y1 y2],'color','r','LineWidth',4);
T = t;
str = strcat(num2str(T)+"/",num2str(T(end)));
timestamp = text(1.6,-0.1,str,'HorizontalAlignment','right');
axis([-3 3 -1 3])
for i = 1:N
th1 = x(i,1);th2 = x(i,2);
x1 = l1*sin(th1);y1 = -l1*cos(th1);
x2 = x1 + l2*sin(th1+th2);y2 = y1 - l2*cos(th1+th2);
set(link1,'xdata',[0 x1],'ydata',[0 y1]);
set(link2,'xdata',[x1 x2],'ydata',[y1 y2]);
T = t(i);
str = strcat(num2str(T)+"/",num2str(t(end)));
set(timestamp,'String',str);
pause(0.02);
drawnow
end
end