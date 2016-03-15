function [u p h A B c k K uA uB pA pB u0 u0A u0B p0 p0A p0B]=CELTang(H,X,A0,N2,omega,f,nu,x,y,z,theta)
% USAGE:  [u p h A B c k]=CELTang(H,X,A0,B0,N2,omega,f,nu,x,z)
% Solve the Coupling Equation for Linear Tides (CELT)
%
% The solution assumes linear, Boussinesq, hydrostate, f-plane mechanics.
% The barotropic mode propagates as a shallow water wave.
%
% INPUTS:
% H    [Nx+1 x 1] Height of flats (positive)
% X    [Nx  x 1]  Location of steps
% A0   [Nm  x 1]  Intenral-tide forcing from left
% N2   [Nx+1 x Nz]  Stratification (from shallow to deep)
% omega [1  x 1]  Frequency of waves
% f    [1   x 1]   Intertial frequency
% nu   [1   x 1]   Vertical viscosity
% x    [nx  x 1]  Horizontal coordinates of output (across slope)
% y    [ny  x 1]  Horizontal coordinates of output (along slope)
% z    [nz  x 1]  Vertical coordinates of output (positive)
% theta [1  x 1]  Angle of obliquity, 0 = normal (degrees)
%
% OUTPUTS:
% u    [nx x nz]  Complex amplitude of internal-tide velocity
% p    [nx x nz]  Complex amplitude of internal-tide pressure
% h    [nx x  1]  Topography mapped to output coordinates
% A    [Nm x Nx] Amplitudes of right-going waves
% B    [Nm x Nx] Amplitudes of left-going waves
% c    [Nm x Nx]  group speed
% k    [Nm x Nx]  Along-slope wavenumber
% K    [Nm x Nx]  Total wavenumber
% uA   [nx x nz]  Complex amplitude of right-going internal-tide velocity
% uB   [nx x nz]  Complex amplitude of left-going internal-tide velocity
% pA   [nx x nz]  Complex amplitude of right-going internal-tide pressure
% pB   [nx x nz]  Complex amplitude of left-going internal-tide pressure
% u0   [nx x ny]  Complex amplitude of surface internal-tide velocity
% u0A  [nx x ny]  Complex amplitude of right-going surface internal-tide velocity
% u0B  [nx x ny]  Complex amplitude of left-going surface internal-tide pressure
% p0   [nx x ny]  Complex amplitude of surface internal-tide pressure
% p0A  [nx x ny]  Complex amplitude of right-going surface internal-tide pressure
% p0B  [nx x ny]  Complex amplitude of left-going surface internal-tide pressure
%
% Sam Kelly, 16 MAY 2014 (smkelly@d.umn.edu)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nm=length(A0);
Nz=size(N2,2);
Nx=length(X);
H
X
nx=length(x);
ny=length(y);
nz=length(z);

ii=complex(0,1);
g=9.81;
Z=linspace(-max(H),0,Nz+1)';
Z=(Z(2:end)+Z(1:end-1))/2;
dz=mean(diff(Z));
B0=A0*0;

% Check if there's enough vertical resolution 
Nm0=Nz-dsearchn(Z,-min(H(H~=0)));
Nm0
Nz
%if Nm>Nm0-1
%   
%    disp(['ERROR: Not enough vertical resolution to match ',num2str(Nm),' modes'])
%    u=NaN; p=[]; h=[]; A=[]; B=[]; c=[]; k=[]; K=[]; uA=[]; uB=[]; pA=[]; pB=[];    
%    error('Die')
%    return
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate sructure functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing vertical modes');

N2=N2.';
N2=flipud(N2);
for i=1:Nx;
    if (H(i)~=0) & (H(i+1)~=0) 
        H(i)
        H(i+1)
        % Get bottom index for Right side
        ind.b(i)=dsearchn(Z,-H(i+1));
        indL=ind.b(i)
        % get the bottom index of Left side
        indR=dsearchn(Z,-H(i))
        
        % need the number to be OK for the shallowest...     
        Nnz=ceil((6*Nm+1)*(max([H(i),H(i+1)])/min([H(i),H(i+1)])))
        
        
        Zz=linspace(Z(min([indR, indL])),Z(end),Nnz)';        
        'Hey'
        size(N2)
        size(Zz)
        N2z = interp1(Z,N2(:,i),Zz(Zz>=Z(indR)));
        dzz{i}=median(diff(Zz));
        % Calculate modes: solve eigenvalue problem with depth-varying stratification
        [phiR{i} c{i}]=MODES(dzz{i},N2z,omega);
        % Keep necessary modes
        phiR{i}=phiR{i}(:,1:Nm)';
        c{i}=c{i}(1:Nm);
        phiL{1}=phiR{1};
        'Hi'
        % now get modes on the next step that match with this
        % note this can slightly change the bathymetry here.
        N2z = interp1(Z,N2(:,i),Zz(Zz>=Z(indL)));
        size(N2z)
        [phiL{i+1} c{i+1}]=MODES(dzz{i},N2z,omega);
        % Keep necessary modes
        phiL{i+1}=phiL{i+1}(:,1:Nm)';
        c{i+1}=c{i+1}(1:Nm);
        nzL = size(phiL{i+1},2)
        nzR = size(phiR{i},2)
        
        if nzL<Nnz
          phiL{i+1}=[zeros([Nm Nnz-nzL]) phiL{i+1}];
        elseif nzR<Nnz
          phiR{i}=[zeros([Nm Nnz-nzR]) phiR{i}];
        end
        % now swap! if getting deeper rather than shallower
        #if H(i)<H(i+1)
        #    phiL2{i+1}=phiR{i}
        #    phiR2{i}=phiL{i+1}
        #end
        %        c{i+1}=c{i+1}(1:Nm);
        zmode{i}=Zz;
        zmode{i+1}=Zz;
        
        %        zmode{1}=zmode{2};
        %aa=zmode{600}
        i
        % Convert from egienspeed to group speed and wavenumber
        c{i}=sqrt(1-f^2/omega^2)*c{i};
        c{i}(1)=sqrt(1-f^2/omega^2)*sqrt(g*H(i));
        K{i}=(1-f^2/omega^2)*omega./c{i};
        K{i}(1)=(1-f^2/omega^2)*omega/c{i}(1);
        c{i+1}=sqrt(1-f^2/omega^2)*c{i+1};
        c{i+1}(1)=sqrt(1-f^2/omega^2)*sqrt(g*H(i+1));
        K{i+1}=(1-f^2/omega^2)*omega./c{i+1};
        K{i+1}(1)=(1-f^2/omega^2)*omega/c{i+1}(1);
        i
    else
        phiR{i}=zeros([Nm Nz]);
        phiL{i+1}=zeros([Nm Nz]);
        c{i}=zeros([Nm 1]);
        K{i}=zeros([Nm 1]);
        K{i+1}=zeros([Nm 1]);
        c{i+1}=zeros([Nm 1]);
        
    end
    
   
    % Get across-step wavenumber
    if i==1
        [val ind.k]=max(A0);
        l=K{i}(ind.k)*sin(theta/180*pi);        
    end
    k{i}=(K{i}.^2-l^2).^(1/2);
    k{i+1}=(K{i+1}.^2-l^2).^(1/2);
    'Hi there'

    PROGRESS_BAR(i,1:Nx+1);
end
'Hi'
phiR{Nx+1}=phiL{Nx+1};

whos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up systems of equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating matrix coefficients');
G=zeros([Nm*Nx*2 Nm*(Nx+1)*2]);
F=zeros([Nm*Nx*2 1]);

% Cycle through steps 
for i=1:Nx  
    i
    Nx
    size(phiL{i+1})
    size(phiR{i})
    N20=mean(N2(:,i));

    % Viscous attenuation of incoming wave from left
    if i==1; rL=ones([Nm 1]); else
        dx=X(i)-X(i-1);
        rL=exp(-(dx./(c{i}.*k{i}./K{i})).*nu.*(1-f^2/omega^2).*(N20-omega^2)./c{i}.^2);
    end
    % Viscous attenuation of incoming wave from right
    if i==Nx; rR=ones([Nm 1]); else
        dx=X(i+1)-X(i);
        rR=exp(-(dx./(c{i+1}.*k{i+1}./K{i+1})).*nu.*(1-f^2/omega^2).*(N20-omega^2)./c{i+1}.^2);
    end    
    % Velocity and Pressure coefficients
    if H(i)>H(i+1) % PART I: Getting shallower        
        % left velocity constraints
        u1=diag(rL.*(k{i}+f*l/(ii*omega))./K{i}.*exp(ii*X(i)*k{i})); % left incoming wave (A)
        u2=diag((k{i}-f*l/(ii*omega))./K{i}.*exp(-ii*X(i)*k{i})); % left outgoing wave (B)      
        
        % right velocity constraints
        if H(i+1)==0 % If vertical wall                   
            u3=zeros(Nm); % right outgoing wave (A)            
            u4=zeros(Nm); % right incoming wave (B)
        else 
            for n=1:Nm
                for m=1:Nm
                    u3(n,m)=sum(phiL{i+1}(m,:).*phiR{i}(n,:)/H(i)*dzz{i},2).*(k{i+1}(m)+f*l/(ii*omega))./K{i+1}(m).*exp(ii*X(i)*k{i+1}(m));  % right outgoing wave (A)  
                    u4(n,m)=sum(phiL{i+1}(m,:).*phiR{i}(n,:)/H(i)*dzz{i},2).*rR(m).*(k{i+1}(m)-f*l/(ii*omega))./K{i+1}(m).*exp(-ii*X(i)*k{i+1}(m)); % right incoming wave (B)
                end
            end
        end
        '1'
    
        % left pressure constraints
        if H(i+1)==0 % If vertical wall (treat as if flat)
            p1=diag(rL.*c{i}.*exp(ii*X(i)*k{i})); % left incoming wave (A)
            p2=diag(c{i}.*exp(-ii*X(i)*k{i})); % left outgoing wave (B)
        else
            for n=1:Nm
                for m=1:Nm
                    p1(n,m)=sum(c{i}(m)*phiR{i}(m,:).*phiL{i+1}(n,:)/H(i+1)*dzz{i},2).*rL(m).*exp(ii*X(i)*k{i}(m)); % left incoming wave (A)
                    p2(n,m)=sum(c{i}(m)*phiR{i}(m,:).*phiL{i+1}(n,:)/H(i+1)*dzz{i},2).*exp(-ii*X(i)*k{i}(m)); % left outgoing wave (B)
                end
            end
        end
        '2'
        size(c{i})
        size(k{i})
        size(c{i+1})
        size(k{i+1})
        % right pressure constraints
        if H(i+1)==0 % If vertical wall (treat as if a mirror)
            p3=-diag(c{i}.*exp(ii*X(i)*k{i})); % right outgoing wave (A)
            p4=-diag(rL.*c{i}.*exp(-ii*X(i)*k{i})); % right incoming wave (B)
        else
            p3=diag(c{i+1}.*exp(ii*X(i)*k{i+1})); % right outgoing wave (A)
            p4=diag(rR.*c{i+1}.*exp(-ii*X(i)*k{i+1})); % right incoming wave (B)
        end       
     '3'   
    else % PART II: Getting deeper  
        'DEEPER'
        i
        size(phiL{i})
        size(phiR{i+1})
        % left velocity constraints
        if H(i)==0 % If vertical wall 
            u1=zeros(Nm); % left outgoing wave (A)
            u2=zeros(Nm); % left incoming wave (B)
        else
            for n=1:Nm
                for m=1:Nm
                    u1(n,m)=sum(phiL2{i}(m,:).*phiR2{i+1}(n,:)/H(i+1)*dzz{i},2).*rL(m).*(k{i}(m)+f*l/(ii*omega))./K{i}(m).*exp(ii*X(i)*k{i}(m)); % left incoming wave (A)
                    u2(n,m)=sum(phiL2{i}(m,:).*phiR2{i+1}(n,:)/H(i+1)*dzz{i},2).*(k{i}(m)-f*l/(ii*omega))./K{i}(m).*exp(-ii*X(i)*k{i}(m)); % left outgoing wave (B)
                end
            end
        end
        '4'
        % right velocity constraints
        u3=diag((k{i+1}+f*l/(ii*omega))./K{i+1}.*exp(ii*X(i)*k{i+1})); % right outgoing wave (A)
        u4=diag(rR.*(k{i+1}-f*l/(ii*omega))./K{i+1}.*exp(-ii*X(i)*k{i+1})); % right incoming wave (B)
        '5'
        % left pressure contraints
        if H(i)==0 % If vertical wall (treat as if a mirror)
            p1=-diag(rR.*c{i+1}.*exp(ii*X(i)*k{i+1})); % left incoming wave (A)
            p2=-diag(c{i+1}.*exp(-ii*X(i)*k{i+1})); % left outgoing wave (B)
        else            
            p1=diag(rL.*c{i}.*exp(ii*X(i)*k{i}));% left incoming wave (A)
            p2=diag(c{i}.*exp(-ii*X(i)*k{i}));% left outgoing wave (B)
        end
        '6'     
        % right pressure contraints
        if H(i)==0 % If vertical wall (treat as if flat)
            p3=diag(c{i+1}.*exp(ii*X(i)*k{i+1})); % right outgoing wave (A)
            p4=diag(rR.*c{i+1}.*exp(-ii*X(i)*k{i+1})); % right incoming wave (B)
        else
            for n=1:Nm
                for m=1:Nm
                    p3(n,m)=sum(c{i+1}(m)*phiR{i+1}(m,:).*phiL{i}(n,:)/H(i)*dzz{i},2).*exp(ii*X(i)*k{i+1}(m)); % right outgoing wave (B)
                    p4(n,m)=sum(c{i+1}(m)*phiR{i+1}(m,:).*phiL{i}(n,:)/H(i)*dzz{i},2).*rR(m).*exp(-ii*X(i)*k{i+1}(m)); % right incoming wave (B)
                end
            end
        end
    end
    'Hi'
    % Fill in the matrix of coefficient matrices
    ind.row=(i-1)*2*Nm+1:i*2*Nm;
    ind.col=(i-1)*2*Nm+[1:4*Nm];
    G(ind.row,ind.col)=[[u1  u2 -u3 -u4];...
                        [p1 -p2 -p3  p4]];
    
    % Fill in forcing vector
    if i==1
        F(1:2*Nm)=F(1:2*Nm)+[[-u1*A0];...
                   [-p1*A0]];
    end
    if i==Nx
        F(end+1-2*Nm:end)= F(end+1-2*Nm:end)+[[ u4*B0];...
                           [-p4*B0]];
    end
    
    % Track progress
    PROGRESS_BAR(i,1:Nx);
end

'Done'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=G(:,Nm+1:end-Nm); % Trim coefficient matrix

disp(['Inverting ',num2str(size(G,1)),' x ',num2str(size(G,2)),' matrix']);pause(.01);
AB=G\F; % Solve matrix problem
AB=[A0; AB; B0]; % Add forcing to output

% Parse solution
AB=reshape(AB.',[2*Nm Nx+1]);
A=AB(1:Nm,:);
B=AB(Nm+1:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Constructing the fields');pause(.01);
h=zeros([nx  1]);
u=zeros([nx nz]);
uA=zeros([nx nz]);
uB=zeros([nx nz]);
p=zeros([nx nz]);
pA=zeros([nx nz]);
pB=zeros([nx nz]);

u0=zeros([nx ny]);
u0A=zeros([nx ny]);
u0B=zeros([nx ny]);
p0=zeros([nx ny]);
p0A=zeros([nx ny]);
p0B=zeros([nx ny]);
'Hi'

% New structure functions for output
for i=1:Nx+1;
    i
    phi0{i}=zeros([Nm nz]);
    for n=1:Nm;
        n
        size(zmode{i})
        size(phiR{i}(n,:))
        pp =interp1(zmode{i},phiR{i}(n,:),-z);
        size(pp)
        phi0{i}(n,:)=pp;
    end
end
'Assemble'
% Assemble fields
for i=1:Nx+1
    N20=mean(N2(:,i));
        
    % Find locations between steps
    if i==1 
        ind.x=x<X(1);         
    elseif i==Nx+1
        ind.x=X(Nx)<=x; 
    else
        ind.x=X(i-1)<=x & x<X(i);
    end
    h(ind.x)=H(i);
   
    for n=2:Nm
        % Calculate viscous attenuation
        if i==1
            Aeps=1;  
            Beps=exp(-(X(i)-x(ind.x))./(c{i}(n).*k{i}(n)./K{i}(n))*nu.*(1-f^2/omega^2).*(N20-omega^2)./c{i}(n).^2);
        elseif i==Nx+1
            Aeps=exp(-(x(ind.x)-X(i-1))./(c{i}(n).*k{i}(n)./K{i}(n))*nu.*(1-f^2/omega^2).*(N20-omega^2)./c{i}(n).^2);
            Beps=1;
        else            
            Aeps=exp(-(x(ind.x)-X(i-1))./(c{i}(n).*k{i}(n)./K{i}(n))*nu.*(1-f^2/omega^2).*(N20-omega^2)./c{i}(n).^2);
            Beps=exp(-(X(i)-x(ind.x))./(c{i}(n).*k{i}(n)./K{i}(n))*nu.*(1-f^2/omega^2).*(N20-omega^2)./c{i}(n).^2);
        end
        % Assemble fields
        if 1
            % Calculate depth profiles
            uA(ind.x,:)=uA(ind.x,:)+(k{i}(n)+f*l/(ii*omega))./K{i}(n)*A(n,i)*Aeps.*exp( ii*x(ind.x)*k{i}(n))*phi0{i}(n,:);
            uB(ind.x,:)=uB(ind.x,:)+(k{i}(n)-f*l/(ii*omega))./K{i}(n)*B(n,i)*Beps.*exp(-ii*x(ind.x)*k{i}(n))*phi0{i}(n,:);
            u(ind.x,:)=uA(ind.x,:)+uB(ind.x,:);
            
            pA(ind.x,:)=pA(ind.x,:)+c{i}(n)*A(n,i)*Aeps.*exp( ii*x(ind.x)*k{i}(n))*phi0{i}(n,:);
            pB(ind.x,:)=pB(ind.x,:)-c{i}(n)*B(n,i)*Beps.*exp(-ii*x(ind.x)*k{i}(n))*phi0{i}(n,:);
            p(ind.x,:)=pA(ind.x,:)+pB(ind.x,:);  
                              
            % Calculate surface profiles
            u0A(ind.x,:)=u0A(ind.x,:)+(k{i}(n)+f*l/(ii*omega))./K{i}(n)*A(n,i)*Aeps.*exp( ii*x(ind.x)*k{i}(n))*exp(ii*y'*(K{i}(n).^2-k{i}(n)^2).^(1/2));
            u0B(ind.x,:)=u0B(ind.x,:)+(k{i}(n)-f*l/(ii*omega))./K{i}(n)*B(n,i)*Beps.*exp(-ii*x(ind.x)*k{i}(n))*exp(ii*y'*(K{i}(n).^2-k{i}(n)^2).^(1/2));
            u0(ind.x,:)=u0A(ind.x,:);+u0B(ind.x,:);            
            
            p0A(ind.x,:)=p0A(ind.x,:)+c{i}(n)*A(n,i)*Aeps.*exp( ii*x(ind.x)*k{i}(n))*exp(ii*y'*(K{i}(n).^2-k{i}(n)^2).^(1/2));
            p0B(ind.x,:)=p0B(ind.x,:)-c{i}(n)*B(n,i)*Beps.*exp(-ii*x(ind.x)*k{i}(n))*exp(ii*y'*(K{i}(n).^2-k{i}(n)^2).^(1/2));
            p0(ind.x,:)=p0A(ind.x,:)+p0B(ind.x,:); 
        else
            u(ind.x,:)=u(ind.x,:)+A(n,i)*Aeps.*exp( ii*x(ind.x)*k{i}(n))*phi0{i}(n,:);
            p(ind.x,:)=p(ind.x,:)+c{i}(n)*A(n,i)*Aeps.*exp( ii*x(ind.x)*k{i}(n))*phi0{i}(n,:);
        end
    end
end
k=cell2mat(k);
K=cell2mat(K);
c=cell2mat(c);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [phi C]=MODES(dz,N2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi C]=MODES(dz,N2,omega)
% USAGE: [phi C]=MODES(dz,N2,omega)
% Obtain vertical modes for arbitrary stratification
%
% INPUTS:
% dz  [1 x 1]   vertical spacing (positive)
% N2  [Nz x 1]  stratification 
%
% OUTPUTS:
% phi [Nz x Nz]  orthonormal pressure and velocity structure eigenfunctions 
% C   [Nz x  1]  eigenspeeds 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(N2);
H=dz*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second derivative matrix (2nd order accurate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D2=zeros(N-1,N-1);
for i=2:N-2
    D2(i,i-1)=1/dz^2;
    D2(i,i)=-2/dz^2;
    D2(i,i+1)=1/dz^2;
end

% Upper boundary condition: rigid lid
D2(1,1)=-2/dz^2; 
D2(1,2)=1/dz^2;

% Bottom boundary condition: flat bottom
D2(N-1,N-1)=-2/dz^2;
D2(N-1,N-2)=1/dz^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "A" Matrix
N2_tmp=(N2(1:end-1)+N2(2:end))/2;
A=diag(-(N2_tmp-omega^2)); % Non hydrostatic
%A=diag(-(N2_tmp));  % Hydrostatic

% Solve generalized eigenvalue problem
[phi k2]=eig(D2,A);

% Sort modes by eigenspeed
k2=diag(k2);
k2(k2<0)=Inf;
C=1./sqrt(k2);
C(C>1000)=0; % Remove modes with c faster 1000 m/s
[C,ind]=sort(C,1,'descend');
phi=phi(:,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive U and P structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add in boundary values
phi=[zeros([1 N-1]); phi; zeros([1 N-1])]; 

% Take derivative to get U and P structure
phi=-diff(phi,1)./dz;

% Add in surface mode
phi=[ones([N 1]) phi];
C=[1; C]; % Leave mode-0 eigenspeed for later

% Normalize 
A=repmat(sum(phi.^2.*dz,1)./H,[N 1]).^(1/2);
A(A==0)=Inf;
phi=phi./A;
phi(:,phi(N,:)<0)=-phi(:,phi(N,:)<0);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRESS_BAR(ind,all_inds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PROGRESS_BAR(ind,all_inds)
% vital function to show progress, needs update to hour glass or spinning wheel
global counters

% Initiate
if ind==all_inds(1)    
    disp(['0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100%'])
    counters.percent=0; 
    counters.N=1;
end

% Print progress
temp=round(counters.N/length(all_inds)*100);
while counters.percent<=temp
    fprintf('|');
    counters.percent=counters.percent+1;    
end
counters.N=counters.N+1;

% Closeout
if ind==all_inds(end)
    clear global counters
    fprintf('\n'); 
end

return

