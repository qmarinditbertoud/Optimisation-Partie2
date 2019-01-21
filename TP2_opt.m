clear all;
close all;

%% question 1
Im_degrade = load('cameraman_bruit.mat');
I = Im_degrade.y;

figure();
imagesc(I);
colormap('gray');
title('image bruit√©e');

%% question 0
[z1,z2] = op_reg(I);
%affiche le gradient et permet d'isoler les contours du bruit et fait
%ressortir les grandes variations sur z1 et les petites variations sur z2

figure();
imagesc(z1);
colormap('gray');
figure();
imagesc(z2);
colormap('gray');

%% 


%fonction diff de gradient lipschitz : 1/2|| y + L*xu||^2
%gradient = L(L*( u) +y)
%gamma = 1 ;

N=500;
lambda=45;

gamma=1.1/lambda;

prox_F = @(x) x - gamma.*prox_L1(x/gamma,lambda/gamma);
gammaX=1;
lambdaX=1;
X1=I;
X2=I;

K=[];

for k=1:N
    
    [grad1,grad2] = op_reg(op_reg_adj(X1,X2)+I);
    
    X1 = (1-lambdaX) * X1 + lambdaX*prox_F(X1-gammaX*grad1);
    X2 = (1-lambdaX) * X2 + lambdaX*prox_F(X2-gammaX*grad2);
    K=[K 0.5*norm(op_reg_adj(X1,X2)).^2+lambda.*norm(op_reg(I+op_reg_adj(X1,X2)),1).^2]; %Crit√®re K(xk) de l'algo
end
  
figure()
plot([1:N],K),title('CritËre K(X^[K]) en fonction du nombre d''itÈrations')

X = op_reg_adj(X1,X2);
u = I + X;

figure();
subplot(1,2,1)
imagesc(u);
colormap('gray');
title('image reconstruite');


subplot(1,2,2)
imagesc(I);
colormap('gray');
title('image bruit√©e');
