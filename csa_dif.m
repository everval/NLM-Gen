function [csax] = csa_dif(x,p,q)
%Generate long memory by using the MA representation of the cross-sectional
%aggregated process using the fast Fourier algorithm, Vera-Valdes(2019)

iT = size(x,1);
n = 2.^nextpow2(2*iT-1);
coefs = ( beta(p+(0:iT-1),q) ./ beta(p,q) ).^(1/2);
csax = ifft(fft(x, n).*fft(coefs', n)); 
csax = csax(1:T, :);
end