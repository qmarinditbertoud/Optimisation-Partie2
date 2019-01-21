function y = op_direct(x,H,ind)

[n,m] = size(H);
xh = real(ifft2(H.*fft2(reshape(x,n,m))));
y = xh(ind)';