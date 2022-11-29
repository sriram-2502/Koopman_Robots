function Fx = PendCartModel_F(t, x, u)
% states are [x xdot theta thetadot]
mp = 2; mc = 8;
l = 0.5; g = 9.8; 
Sx3 = sin(x(3,:));
Cx3 = cos(x(3,:));

M = mc+mp.*(Sx3.^2);
%% Dynamics wtih control xdot = f(t,x,u)

f1 = (-mp*l.*Sx3.*x(4,:).^2 + u + mp*g.*Cx3.*Sx3) ./ M;
f2 = (-mp*l.*Cx3.*Sx3.*x(4,:).^2 + u.*Cx3 + (mp+mc)*g.*Sx3) ./ (M.*l);

Fx = [x(2,:); f1; x(4,:); f2];
end