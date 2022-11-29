clear all; close all; clc;
 
x=linspace(-10,10,100); 
t=linspace(0,4*pi,80); dt=t(2)-t(1);
[X,T]=meshgrid(x,t); 

%f1=sech(X+3).*(1-0.5*cos(2*T));
%f2=(sech(X).*tanh(X)).*(2.5*sin(3*T));
f1=sech(X+3).*(1*exp(i*2.3*T));
f2=(sech(X).*tanh(X)).*(2*exp(i*2.8*T));
f=f1+f2;
 
subplot(2,2,1), surfl(real(f1)), shading interp, colormap(gray), view(-20,60)
subplot(2,2,2), surfl(real(f2)), shading interp,colormap(gray), view(-20,60)
subplot(2,2,3), surfl(real(f)), shading interp,colormap(gray), view(-20,60) 

% PCA decomposition
 
[u,s,v]=svd(f');
 
figure(2)
subplot(4,3,1), plot(diag(s)/sum(diag(s)),'ko','Linewidth',[2])
subplot(4,1,2), plot(t,v(:,1)/max(v(:,1)),t,v(:,2)/max(v(:,2)),'Linewidth',[2])
subplot(4,1,3), plot(x,u(:,1)/max(u(:,1)),'Linewidth',[2])
subplot(4,1,4), plot(x,u(:,2)/max(u(:,2)),'Linewidth',[2])

X = f.'; X1 = X(:,1:end-1); X2 = X(:,2:end);
r=2;
[U2,Sigma2,V2] = svd(X1,'econ'); U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
 
% DMD J-Tu decomposition:  Use this one
    
Atilde = U'*X2*V/Sigma;    
[W,D] = eig(Atilde);    
Phi = X2*V/Sigma*W;
    
mu = diag(D);
omega = log(mu)/dt;

u0=f(1,:).';
y0 = Phi\u0;  % pseudo-inverse initial conditions
u_modes = zeros(r,length(t));
for iter = 1:length(t)
     u_modes(:,iter) =(y0.*exp(omega*t(iter)));
end
u_dmd = Phi*u_modes;

figure(1)
subplot(2,2,4)
surfl(real(u_dmd.')), shading interp,colormap(gray), view(-20,60)


figure(2), 
subplot(4,3,2), plot(diag(Sigma2)/sum(diag(Sigma2)),'ko','Linewidth',[2])
subplot(4,1,3), hold on, plot(x,-Phi(:,1)/max(Phi(:,1)),'Linewidth',[2])
subplot(4,1,4), hold on, plot(x,Phi(:,2)/max(Phi(:,2)),'Linewidth',[2])


subplot(4,1,3), plot(x,-sech(x).*tanh(x)/max(sech(x).*tanh(x)),':','Linewidth',[2])
subplot(4,1,4), plot(x,sech(x+3)/max(sech(x+3)),':','Linewidth',[2])

E1(1)=norm(u(:,1)/max(u(:,1))-(-sech(x).*tanh(x)/max(sech(x).*tanh(x))).');
E2(1)=norm(u(:,2)/max(u(:,2))-(sech(x+3)/max(sech(x+3))).');

E1(2)=norm(-Phi(:,1)/max(Phi(:,1))-(-sech(x).*tanh(x)/max(sech(x).*tanh(x))).');
E2(2)=norm(Phi(:,2)/max(Phi(:,2))-(sech(x+3)/max(sech(x+3))).');

subplot(4,3,3), bar([E1; E2].'), 

%%
figure(1)
for j=1:4
subplot(2,2,j)
set(gca,'Ylim',[0 80],'Ytick',[0 20 40 60 80],'Xlim',[0 100],'Xtick',[0 50 100],'Xticklabel',{'-10','0','10'}, ...
    'Zlim',[-1.5 1.5],'Ztick',[-1.5 0 1.5],'Fontsize',[13])
end

figure(2)
subplot(4,3,1)
set(gca,'Xlim',[0 20],'Xtick',[0 10 20],'Ylim',[0 0.7],'Ytick',[0 0.2 0.4 0.6],'Yticklabel',{'0%','20%','40%','60%'},'Fontsize',[13])

subplot(4,3,2)
set(gca,'Xlim',[0 20],'Xtick',[0 10 20],'Ylim',[0 0.7],'Ytick',[0 0.2 0.4 0.6],'Yticklabel',{'0%','20%','40%','60%'},'Fontsize',[13])

subplot(4,3,3)
set(gca,'Xlim',[0 4],'Xtick',[1 2 3],'Ylim',[0 3.5],'Ytick',[0 1 2 3],'Fontsize',[13])
legend('E_1','E_2')

subplot(4,1,2)
set(gca,'Xlim',[0 4*pi],'Xtick',[0 pi 2*pi 3*pi 4*pi],'Xticklabel',{'0','\pi','2\pi','3\pi','4\pi'},'Fontsize',[13])
legend('PCA 1','PCA 2','ICA 1','ICA 2')


subplot(4,1,3)
set(gca,'Xlim',[-10 10],'Xtick',[-10 -5 0 5 10],'Ylim',[-1.5 1],'Ytick',[-1 0 1],'Fontsize',[13])
legend('PCA','DMD','Exact')


subplot(4,1,4)
set(gca,'Xlim',[-10 10],'Xtick',[-10 -5 0 5 10],'Ylim',[-0.8 1],'Ytick',[-0.5 0 0.5 1],'Fontsize',[13])


