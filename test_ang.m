% test convergence of solution vs dx for different slopes
clear

omega=1;
f=0;
nu=0;%10^(-5); % or set to 0 for no viscosity
rho0=1000;
Nm=32;
Nx=64;
ang=57;

% Slope and Nx
s=1; % Slope criticality
H=2/3;
x2=H/s;    

% Output grid
x=linspace(-x2*2,x2,300)';
y=linspace(-1,1,200)';
z=linspace(0,1,100)';

% Stratification
Nz=max([Nm/(1-H)+50 100]);
N2=2*ones([Nx+1 Nz]);
               
% Topography
H0=linspace(1,1-H,Nx+1)';
X0=linspace(-x2,0,Nx+1)';X0=X0(1:end-1);
                     
% Forcing
A0=zeros([Nm 1]);A0(2)=1;
B0=zeros([Nm 1]);
        
% Call CELT model
[dat1.u dat1.p dat1.h dat1.A dat1.B dat1.c]=CELT(H0,X0,A0,B0,N2,omega,f,nu,x,z);
[dat2.u dat2.p dat2.h dat2.A dat2.B dat2.c dat2.k dat2.K]=CELTang(H0,X0,A0,N2,omega,f,nu,x,y,z,ang);

% Incident energy flux
dat1.Finc=nansum(1/2*rho0*H0(1)*dat1.c(:,1).*abs(A0).^2);
dat2.Finc=nansum(1/2*rho0*H0(1)*real(dat2.k(:,1)./dat2.K(:,1)).*dat2.c(:,1).*abs(A0).^2); % The incident flux is less by k/K because the wave is oblique

% Right going energy flux
dat1.FA=1/2*rho0*repmat(H0',[Nm 1]).*dat1.c.*abs(dat1.A).^2;
dat2.FA=1/2*rho0*repmat(H0',[Nm 1]).*real(dat2.k./dat2.K).*dat2.c.*abs(dat2.A).^2;

% Left going energy flux
dat1.FB=1/2*rho0*repmat(H0',[Nm 1]).*dat1.c.*abs(dat1.B).^2;
dat2.FB=1/2*rho0*repmat(H0',[Nm 1]).*real(dat2.k./dat2.K).*dat2.c.*abs(dat2.B).^2;

% Conservation of energy 
% residual = energy flux in - energy flux out
dat1.res=nansum(dat1.FA(:,1)+dat1.FB(:,end))-nansum(dat1.FB(:,1)+dat1.FA(:,end));
dat2.res=nansum(dat2.FA(:,1)+dat2.FB(:,end))-nansum(dat2.FB(:,1)+dat2.FA(:,end));

disp(sprintf(['Normal energy loss: ',num2str(dat1.res),' W \nOblique energy loss: ',num2str(dat2.res),' W']));

% % Total generation
% gen=nansum(FB(2:end,1)+FA(2:end,end));
% 
% % Mode-1 generation
% gen=FB(2,1)+FA(2,end);

%%
ulims=[-1 1]*1.5;
xlims=[min(x) max(x)];

close(figure(1));figure(1);clf;colormap(gray)

subplot(2,1,1);
pcolor(x,z,real(dat1.u)');axis ij;shading interp;caxis(ulims);hold on;
set(gca,'tickdir','out','xlim',xlims,'ytick',0:.25:5)
xlabel('Distance')
ylabel('Depth')
text(.0,1.08,'Velocity at t=0')
colorbar

subplot(2,1,2);
pcolor(x,z,real(dat2.u)');axis ij;shading interp;caxis(ulims);hold on;
set(gca,'tickdir','out','xlim',xlims,'ytick',0:.25:5)
xlabel('Distance')
ylabel('Depth')
text(.0,1.08,'Velocity at t=0')
colorbar

