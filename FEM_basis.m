function phi = FEM_basis(X,x)
% calculate the FEM basis functions phi(x,k). Each column is the value of
% one function phi(:,k). 

if size(X,1) == 1
    X = X';
end
if size(x,1) == 1
    x = x';
end

K = length(X);

IdxBegin = find(x <= X(1),1,'last');
IdxEnd = find(x >= X(K),1,'first');

X2 = [x(IdxBegin); X; x(IdxEnd)]; 
phi = zeros(length(x),K);
for k = 1:K
    phi(:,k) = FEM_basis_function(X2(k),X2(k+1),X2(k+2),x);
end


