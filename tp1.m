close all
clear all
clc

load cameraman.mat
load cameraman_cs.mat
load flou.mat
load laplacien.mat
load decimation.mat
Y=load('cameraman_bruit.mat');Y=Y.y;

lambda=40;
nu=lambda;
gamma=1.9/nu;

N=50;
lambdak=0.2;
gammak=1;

u1=Y;
u2=Y;
K=[];

prox_F=@(z) z-gamma.*prox_L1(z/gamma,lambda/gamma); %fonction proximal de l'indicatrice de norme infinie


%Algo Forward-Backward
for i=1:N
    [gradu1,gradu2]=op_reg(op_reg_adj(u1,u2)+Y); %grad de g appliqué à u
        
    u1=u1*(1-lambdak)+lambdak*(prox_F(u1-gammak.*gradu1)); %algo F-B appliqué pour u1
    u2=u2*(1-lambdak)+lambdak*(prox_F(u2-gammak.*gradu2)); %algo F-B appliqué pour u2
    
    U=op_reg_adj(u1,u2); %L*(u)
    
    K=[K 1/2*norm(U).^2+lambda.*norm(op_reg(Y+U),1).^2]; %Critère K(xk) de l'algo
end

figure(2)
subplot 121,imagesc(Y),colormap 'gray',title('Image bruitée')
subplot 122,imagesc(Y+U),colormap 'gray',title('Image après algo FB') %Y+U=X puisque U=L*(u)

figure(3)
plot([1:N],K),title('Critère K(X^{[K]}) en fonction du nombre d''itérations')
xlabel('Nombre d''itération')
ylabel('Critère')