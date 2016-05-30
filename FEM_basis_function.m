function phi = FEM_basis_function(xleft,xcenter,xright, x)
% calculate the piecewise linear basis function phi(x) of the FEM method
% xleft: the left end
% xright: the right end
% xcenter: the peak
% x: value at which we evaluate the value of phi. 

EPS = 1e-10;

phi = 0*x;

if (xcenter-xleft > EPS)
    IdxLeft = find(x <= xcenter & x > xleft);
    phi(IdxLeft) = (x(IdxLeft) - xleft)/(xcenter - xleft);
end
if (xright - xcenter > EPS)
    IdxRight = find(x < xright & x >= xcenter);
    phi(IdxRight) = (x(IdxRight) - xright)/(xcenter - xright);
end

