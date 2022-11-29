function g = PendCartModel_g(t, x)
% states are [x xdot theta thetadot]
mp = 2; mc = 8;
l = 0.5; g = 9.8; 
Cx3 = cos(x(3,:));
Sx3 = sin(x(3,:));

M = mc+mp.*(Sx3.^2);
%% g(x) dynamics for xdot = f(t,x)+g(t,x)*u
g = [0; 1/M; 0; Cx3./(M.*l)];

end