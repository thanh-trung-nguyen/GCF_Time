function data_new = noise_truncate(data,noiselevel)
% truncate data points before the first peak and smaller than the noise
% level. This is helpful in the Laplace transform.
% data is a vector.


Nt = length(data);
n = 1;
while n < Nt && abs(data(n)) < noiselevel
    data(n) = 0;
    n = n+1;
end

data_new = data;
