
omega = 2*pi/12.4/3600.
f = -1e-4
nu=1e-4
Nm=10
Nx=50
Nz=floor(5000/10)
rho0=1000.


% get N2 
load ../TasmaniaRays.mat

z = linspace(0,5000.,Nz)';
N2 = interp1(ray.z(~isnan(ray.N2)),ray.N2(~isnan(ray.N2)),z,'linear','extrap')';

N2 = repmat(N2',1,Nx+1)';

H0 = [linspace(5000,500,Nx) 0]';
X0 = [linspace(-10.e3,0.,Nx-1) 30.e3]';
H0 = [linspace(5000,300,Nx+1)]';
X0 = [linspace(-10.e3,0.,Nx)]';
A0 = zeros(Nm,1);
A0(2)=1.;
x = linspace(-100,30.,1000)'*1e3;
y = 0.
size(N2)

%%
ang=0.
[dat2.u dat2.p dat2.h dat2.A dat2.B dat2.c dat2.k dat2.K]=CELTang(H0,X0,A0,N2,omega,f,nu,x,y,z,ang);
dat1=dat2
dat1.Finc=sum(1/2*rho0*H0(1)*dat1.c(:,1).*abs(A0).^2);
dat2.Finc=sum(1/2*rho0*H0(1)*real(dat2.k(:,1)./dat2.K(:,1)).*dat2.c(:,1).*abs(A0).^2); % The incident flux is less by k/K because the wave is oblique

% Right going energy flux
dat1.FA=1/2*rho0*repmat(H0',[Nm 1]).*dat1.c.*abs(dat1.A).^2;
dat2.FA=1/2*rho0*repmat(H0',[Nm 1]).*real(dat2.k./dat2.K).*dat2.c.*abs(dat2.A).^2;

% Left going energy flux
dat1.FB=1/2*rho0*repmat(H0',[Nm 1]).*dat1.c.*abs(dat1.B).^2;
dat2.FB=1/2*rho0*repmat(H0',[Nm 1]).*real(dat2.k./dat2.K).*dat2.c.*abs(dat2.B).^2;

%%

%dat1=dat2


ang=0.
[dat2.u dat2.p dat2.h dat2.A dat2.B dat2.c dat2.k dat2.K]=CELTangJ(H0,X0,A0,N2,omega,f,nu,x,y,z,ang);
dat2.Finc=sum(1/2*rho0*H0(1)*real(dat2.k(:,1)./dat2.K(:,1)).*dat2.c(:,1).*abs(A0).^2); % The incident flux is less by k/K because the wave is oblique

% Right going energy flux
dat2.FA=1/2*rho0*repmat(H0',[Nm 1]).*real(dat2.k./dat2.K).*dat2.c.*abs(dat2.A).^2;

% Left going energy flux
dat2.FB=1/2*rho0*repmat(H0',[Nm 1]).*real(dat2.k./dat2.K).*dat2.c.*abs(dat2.B).^2;
%%
figure(12)
clf
subplot(2,1,1)
imagesc(x,z,real(dat1.u)')
caxis([-1.,1.]*2)
subplot(2,1,2)
imagesc(x,z,real(dat2.u)')
caxis([-1.,1.]*2)

%%
% Incident energy flux
save(sprintf('datJmkM%dNx%d',Nm,Nx),'dat2')

%%
figure(22)
clf
hold on
Nm2 = size(dat2.FB,1)
Nm1 = size(dat1.FB,1)
plot((1:Nm2)-1,dat2.FB(:,1)/dat2.Finc,'b')
plot((1:Nm2)-1,dat2.FA(:,end)/dat2.Finc,'m')
plot((1:Nm1)-1,dat1.FB(:,1)/dat1.Finc,'b--')
plot((1:Nm1)-1,dat1.FA(:,end)/dat1.Finc,'m--')
set(gca,'yscale','log')

sum(dat1.FB(:,1)/dat1.Finc)
sum(dat2.FB(:,1)/dat2.Finc)
%%
figure(23)
close all
%plot((1:Nm)-1,cumsum(dat1.FB(:,1)/dat1.Finc),'r')
hold on
ylim([0.1,1])

plot((1:Nm)-1,cumsum(dat2.FB(:,1)/dat2.Finc))
%set(gca,'yscale','log')

sum(dat1.FB(:,1)/dat1.Finc)

sum(dat2.FB(:,1)/dat2.Finc)
dat2.FB(:,1)/dat2.Finc


