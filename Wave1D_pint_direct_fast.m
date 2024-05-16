%new leapfrog scheme for 1D wave equation
%y_tt=y''+f(x,t),  y(0,t)=y(1,t)=0; y(x,0)=y0(x); y_t(x,0)=y1(x);
clear
global Ix Ax Vs Ds iVs
tic
maxL=10; % increase for larger mesh

%%(B) Set the problem data: rhs, exact solutions, parameters
T=2; xa=0; xb=1;
%Example 1
% y_sol=@(x,t) sin(pi*x).*cos(pi*t);
% y0fun=@(x) y_sol(x,0);
% y1fun=@(x) zeros(size(x));
% ffun=@(x,t) zeros(size(x));
%Example 2
%     y_sol=@(x,t) sin(pi*x).*sin(10*t);
%     y0fun=@(x) y_sol(x,0);
%     y1fun=@(x) 10*sin(pi*x);
%     ffun=@(x,t) (pi^2-100)*sin(pi*x).*sin(10*t);
%Example 3 (nonsmooth IC)
        y_sol=@(x,t) ysolfun(x,t);
        y0fun=@(x) ((x>=3/8)&(x<=5/8)).*(cos(4*pi*(x-0.5)).^2);
        y1fun=@(x) zeros(size(x));
        ffun=@(x,t) zeros(size(x));


Nx=2.^(3:maxL);%nx=ny
Nt=T*2.^(3:maxL);
levmax=length(Nx);
fprintf('(Nx,Nt) \t Error\t\t Order\t  CPU \n');
y_err_0=0;
for s=1:levmax
    nx=Nx(s);m=nx-1; nt=Nt(s);
    dt=T/nt; h=(xb-xa)/nx; xx=xa+h:h:xb-h; %interior nodes in space
    [XX,TT] = ndgrid(xx,dt:dt:T);%not including initial time step

    %one-shot, no time marching
    Ix=speye(m,m); It=speye(nt,nt);
    Ax=(1/h^2)*gallery('tridiag',m,-1,2,-1); % central finite difference
 
    F=ffun(XX,TT); F(:,1)=F(:,1)+y1fun(xx(:))/(2*dt);F(:,2)=F(:,2)-y0fun(xx(:))/(4*dt^2); %adjust rhs

    %A=kron(At,Ix)+kron(It,Ax); %system matrix, no need to construct
    tic 
    [Vs,Ds,iVs,iter]=fasteigB(nt,dt,1e-8);Ds=Ds.^2;
    y_h=PinTinvA(F(:));

    %measure error
    ysol=y_sol(XX,TT);%exact solution
    y_err=norm(y_h(:)-ysol(:),inf);
    fprintf('(%d,%d)&\t %1.2e& \t%1.2f\t&%1.3f \n',...
        nx,nt,y_err,log2(y_err_0/y_err),toc)
    %      fprintf('&%1.2e\n',...
    %         y_err)
    y_err_0=y_err;
end
 
function z=PinTinvA(r)
global Ix Ax Vs Ds iVs
nx=size(Ax,1); nt=length(Ds);
R=reshape(r,nx,nt);
R1=R*(iVs.'); %step (a)
for j=1:nt %parallel in time step (b)
    R1(:,j)=(Ds(j)*Ix+Ax)\R1(:,j); %can use multigrid to speed up
end
z=R1*(Vs.');%step (c)
z=real(z); %drop imaginary part, which should be close to zero
end

function z=ysolfun(x,t)
%series solution for %Example 3 (nonsmooth IC)
z=0;
for n=1:50
    if(n~=8)
        Rn=64*(cos(5*n*pi/8)-cos(3*n*pi/8))/(pi*(n^3-64*n));
        z=z+sin(n*pi*x).*cos(n*pi*t)*Rn;
    end
end
end