%Direct PinT leapfrog scheme for 1D heat equation
%y_t=y_xx+f(x,t),  y(0,t)=y(1,t)=0; y(x,0)=y0(x);
clear
%%(B) Set the problem data: rhs, exact solutions, parameters
T=2; xa=0; xb=1;
y0=@(x) sin(pi*x);
y_sol=@(x,t) exp(-pi^2*t).*sin(pi*x); % require zero boundary condition
f=@(x,t) zeros(size(x,1),size(t,2));


maxL=9; % increase for larger mesh
Nx=2.^(5:maxL);Nt=T*2.^(3:maxL); levmax=length(Nx);
fprintf('(Nx,Nt) \t Error\t\t Order\t  CPU(eig,3-steps) \n');
y_err_0=0; iter=0; b0=[]; condV=[];
for s=1:levmax
    
    nx=Nx(s);m=nx-1; nt=Nt(s);
    dt=T/nt; h=(xb-xa)/nx; xx=xa+h:h:xb-h; %interior nodes in space
    [XX,TT] = meshgrid(xx,dt:dt:T);%not including initial time step
    
    Ix=speye(m,m); It=speye(nt,nt);
    Ax=(1/h^2)*gallery('tridiag',m,-1,2,-1); % central finite difference
    
    %     e_t=ones(nt,1);
    %     At= spdiags([-e_t/2 e_t/2],[-1 1],nt,nt)/dt; %time scheme
    %     At(end,end-1:end)=[-1 1]/dt; %fix last row for backward Euler
    %
    F=f(XX,TT); F(1,:)=F(1,:)+y0(xx)/(2*dt); %adjust rhs for first step
    %A=kron(Ix,At)+kron(Ax,It); %system matrix, no need to construct
    
    
    %y_h=A\F(:); %sparse direct solver
    %eigenvalue decomposition
    %[Vs,Ds] = eig(full(At),'vector'); %factorize At
    
    %use paper's formula to construct Vs
    %     Ds= eig(full(At)); %factorize At
    %     x=(-1i*Ds*dt);
    tic
    %Newton method for finding all eigenvalues
    [Vs,Ds,iVs,iter]=fasteigB(nt,dt,1e-8);
    teig=toc;
    %condV=[condV;norm(Vs)*norm(iVs)];
    %condV=[condV;svds(Vs,1)*svds(iVs,1)];
    tic
    R1=iVs*F;%O(nt^3)
    for j=1:nt %step (b)
        R1(j,:)=((Ds(j)*Ix+Ax)\R1(j,:).').'; %can use multigrid to speed up
    end
    y_h=Vs*R1;%step (c) :%O(nt^3)
    tdiag=toc;
    %measure error
    ysol=y_sol(XX,TT);%exact solution
    y_err=norm(y_h(:)-ysol(:),inf);%maxmum error
    fprintf('(%3d,%3d)&\t %1.2e& \t%1.2f\t&(%1.3f,%1.3f) &\t %d \n',...
        nx,nt,y_err,log2(y_err_0/y_err),teig,tdiag,iter)
    y_err_0=y_err;
end  
