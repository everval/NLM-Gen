csa_dif = function(x, p, q){
  iT = length(x)
  n = nextn(2*iT - 1, 2)
  k = 0:(iT-1)
  b = (beta(p+k,q)/beta(p,q))^(1/2)
  dx = fft( fft(c(x, rep(0, n - iT))) * fft(c(b, rep(0, n - iT))), inverse = T) / n
  return(Re(dx[1:iT]))
}