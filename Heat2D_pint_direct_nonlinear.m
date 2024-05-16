%Direct PinT leapfrog scheme for 2D nonlinear heat equation
%y_t-Delta y+nlfun(y)=f(x,t),  y=0 on BC; y(x,0)=y0(x);
clear
%%(B) Set the problem data: rhs, exact solutions, parameters
T=2; xa=-1; xb=1;
y0=@(x,y) (x-xa).*(x-xb).*(y-xa).*(y-xb);
y_sol=@(x,y,t) exp(-t).*(x-xa).*(x-xb).*(y-xa).*(y-xb);% require zero boundary condition
nlfun=@(z) (z.^3-z); nlfunjac=@(z) (3*z.^2-ones(size(z))); %nonlinear functions
%nlfun=@(z) sqrt(z); nlfunjac=@(z) 0.5./sqrt(z);
f=@(x,y,t) -exp(-t).*(x-xa).*(x-xb).*(y-xa).*(y-xb)...
    -2*exp(-t).*((x-xa).*(x-xb)+(y-xa).*(y-xb))+nlfun(y_sol(x,y,t));
maxL=7; % increase for larger mesh
Nx=2.^(3:maxL);Nt=T*2.^(3:maxL); levmax=length(Nx);
fprintf('(Nx,Nx,Nt) \t Error\t\t Order\t  CPU\t Iter \n');
y_err_0=0;
for s=1:levmax
    tic
    nx=Nx(s);m=nx-1; nt=Nt(s);
    dt=T/nt; h=(xb-xa)/nx; xx=xa+h:h:xb-h; %interior nodes in space
    [X,Y] = meshgrid(xx,xx);
    [XX,YY,TT] = meshgrid(xx,xx,dt:dt:T);%not including initial time step
    
    Ix=speye(m^2,m^2); It=speye(nt,nt);
    Ax=(1/h^2)*gallery('poisson',m);% central finite difference
    
    e_t=ones(nt,1);
    At= spdiags([-e_t/2 e_t/2],[-1 1],nt,nt)/dt; %time scheme
    At(end,end-1:end)=[-1 1]/dt; %fix last row for backward Euler
    
    F=f(XX,YY,TT); 
    F(:,:,1)=F(:,:,1)+y0(X,Y)/(2*dt); %adjust rhs for first step
      
    F=reshape(F,m^2,nt)';  
    %A=kron(Ix,At)+kron(Ax,It);%system matrix, no need to construct
    [Vs,Ds] = eig(full(At),'vector'); %factorize At
    
    %Newton iteration
    maxit=50;tol=1e-8;
    y_h=zeros(nt,m^2);%initial guess
    for iter=1:maxit
        %fjac=spdiags(nlfunjac(mean(y_h)'),0,m^2,m^2); %mean of y_h in time
        fjac=spdiags(mean(nlfunjac(y_h))',0,m^2,m^2); %mean of nlfun'(y_h),bett
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
    ysol=y_sol(XX,YY,TT); ysol=reshape(ysol,m^2,nt)';  %exact solution
 
    y_err=norm(y_h(:)-ysol(:),inf);%maxmum error
    fprintf('(%3d,%3d,%3d)&\t %1.2e& \t%1.2f\t&%1.3f&\t%d \n',...
        nx,nx,nt,y_err,log2(y_err_0/y_err),toc,iter)
    y_err_0=y_err;
end
