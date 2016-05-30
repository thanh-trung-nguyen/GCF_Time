function f = FEM_basis_expansion(phi,CoefValue)
% evaluate the value of 
% f(x) = 1 + sum_{k=1}^K phi_k(x)*beta(k)
% This is used to parametrize the coefficient.

f = ones(size(phi,1),1);
K = length(CoefValue);

for k = 1:K
    f = f + phi(:,k)*CoefValue(k);
end






