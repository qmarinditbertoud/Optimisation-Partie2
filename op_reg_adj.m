function x = op_reg_adj(y1,y2)

H = psf2otf([1/2,-1/2],size(y1));
V = psf2otf([1/2,-1/2]',size(y2));
x = real(ifft2(conj(H).*fft2(y1))) + real(ifft2(conj(V).*fft2(y2)));


