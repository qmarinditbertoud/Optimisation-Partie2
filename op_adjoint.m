function y = op_adjoint(x,H)

y = real(ifft2(conj(H).*fft2(x)));

