function g = Acrobot_g(t, x)
m1 = 1 ; m2 = 1;
l1 = 1 ; l2 = 1;g = 9.81;
lc1 = l1/2; lc2 = l2/2;
I1 = 1 ;I2 = 1; 
q1 = x(1,:); q2 = x(2,:); dq1 = x(3,:); dq2 = x(4,:);

m11 = I1+ I2 + m2*l1^2 + 2*m2*l1*lc2*cos(q2);
m12 = I2 + m2*l1*lc2*cos(q2);
m22 = I2;
M = [m11 m12 ; m12 m22];
B1 = [0;1];

g = [0;0; M\(B1)];
end