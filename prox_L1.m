 function p = prox_L1(x, gamma)


p = sign(x) .* max(0, abs(x) - gamma);