%Direct PinT leapfrog scheme for 1D nonlinear heat equation
%y_t-y_xx+nlfun(y)=f(x,t),  y(0,t)=y(1,t)=0; y(x,0)=y0(x);
clear
%%(B) Set the problem data: rhs, exact solutions, parameters
T=2; xa=0; xb=pi;
y0=@(x) x.*(x-pi);
y_sol=@(x,t) exp(-t).*x.*(x-pi);% require zero boundary condition
nlfun=@(z) -(z-z.^3)/3; nlfunjac=@(z) -(ones(size(z))-3*z.^2)/3; %nonlinear functions
%nlfun=@(z) sqrt(z); nlfunjac=@(z) 0.5./sqrt(z);
f=@(x,t) -exp(-t).*x.*(x-pi)-2*exp(-t)+nlfun(y_sol(x,t));
maxL=8; % increase for larger mesh
Nx=2.^(3:maxL);Nt=T*2.^(3:maxL); levmax=length(Nx);
fprintf('(Nx,Nt) \t\t Error\t\t Order\t  CPU\t Iter \n');
y_err_0=0;
for s=1:levmax
    tic
    nx=Nx(levmax);m=nx-1; nt=Nt(s);
    dt=T/nt; h=(xb-xa)/nx; xx=xa+h:h:xb-h; %interior nodes in space
    [XX,TT] = meshgrid(xx,dt:dt:T);%not including initial time step
    
    Ix=speye(m,m); It=speye(nt,nt);
    Ax=(1/h^2)*gallery('tridiag',m,-1,2,-1); % central finite difference
    
    e_t=ones(nt,1);
    At= spdiags([-e_t/2 e_t/2],[-1 1],nt,nt)/dt; %time scheme
    At(end,end-1:end)=[-1 1]/dt; %fix last row for backward Euler
    
    F=f(XX,TT); F(1,:)=F(1,:)+y0(xx)/(2*dt); %adjust rhs for first step
    %A=kron(Ix,At)+kron(Ax,It);%system matrix, no need to construct
    [Vs,Ds] = eig(full(At),'vector'); %factorize At
    
    %Newton iteration
    maxit=100;tol=1e-8;
    y_h=zeros(nt,m);%initial guess
    for iter=1:maxit
        fjac=spdiags(nlfunjac(mean(y_h)'),0,m,m); %mean of y_h
        %jac=spdiags(mean(nlfunjac(y_h))',0,m,m); %mean of nlfun'(y_h),bett
        Res=y_h*fjac-nlfun(y_h)+F;%residual   
        
        R1=Vs\Res; %step (a)
        for j=1:nt %step (b)
            R1(j,:)=((Ds(j)*Ix+Ax+fjac)\R1(j,:).').';  %parallel in time
        end
        y_h_new=Vs*R1;%step (c)
        if(norm(y_h_new(:)-y_h(:),inf)<tol)
            break;
        end
        y_h=y_h_new;
    end
    %measure error
    ysol=y_sol(XX,TT);%exact solution
    y_err=norm(y_h(:)-ysol(:),inf);%maxmum error
    fprintf('(%3d,%3d)&\t %1.2e& \t%1.2f\t&%1.3f&\t%d \n',...
        nx,nt,y_err,log2(y_err_0/y_err),toc,iter)
    y_err_0=y_err;
end
