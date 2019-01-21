function y = op_direct(x,H)

y = real(ifft2(H.*fft2(x)));

