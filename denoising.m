function denoised_data = denoising(data,t, max_freq)
% denoising a noisy data using Fourier transform


[fourier_dat,freq] = fourier_fft(data,t);

fourier_dat(abs(freq) > max_freq) = 0;

denoised_data = real(fourier_ifft(fourier_dat,freq,t));

if size(data,1) > 1 && size(denoised_data,2) > 1 
    denoised_data = denoised_data.';
end

if size(data,2) > 1 && size(denoised_data,1) > 1 
    denoised_data = denoised_data.';
end
