function q = compute_q_from_laguerre_coef(qn,s,s0)
% compute v from the Laguerre's coefficient

[Nq,Nx] = size(qn);

q = zeros(length(s),Nx);

for k = 0:Nq-1    
    fn = laguerre(k,s,s0);
    
    q = q + fn*qn(k+1,:);
end
