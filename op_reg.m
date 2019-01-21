function [y1,y2] = op_reg(x)

H = psf2otf([1/2,-1/2],size(x));
V = psf2otf([1/2,-1/2]',size(x));
y1 = real(ifft2(H.*fft2(x)));
y2 = real(ifft2(V.*fft2(x)));

