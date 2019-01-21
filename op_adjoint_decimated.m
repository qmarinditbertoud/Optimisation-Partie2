function y = op_adjoint_decimated(x,H,ind)

[n,m] = size(H);
yd = zeros(1,n*m);
yd(ind) = x;

y = real(ifft2(conj(H).*fft2(reshape(yd,n,m))));
y = y(:)';
