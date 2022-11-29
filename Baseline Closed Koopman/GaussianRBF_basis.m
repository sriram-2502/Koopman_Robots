function [Psi, DPsi] = GaussianRBF_basis(C, dim, gamma)
%  [Psi, DPsi] = GaussianRBF_basis(deg, dim) returns a gaussian RBF function and its derivative 
%   There is not a 1 included in Psi, only the linear and higher order terms
%   Centers = matrix of center points (size = (dim of x, K))
%   dim = number of states

x=sym('x',[dim,1]);
assume(x,'real')

Psi = [];
G = [];

for i=1:size(C,2)
    for j=1:dim
        g = exp(-gamma*(norm(x(j)-C(j,i)).^2));
        G = [G; g];
    end
end

Psi = G
DPsi = jacobian(Psi,x);

end