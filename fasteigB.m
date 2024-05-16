function [V,D,iV,iter]=fasteigB(Nt,dt,tol)
% find the eigendecomposition of B
% V=eigen vector matrix; D=a column vector of eigenvalues; iV = inv(V)
if(nargin<3)
    tol=1e-10; %use smaller tol for better accuracy if needed
end
%Newton iterations vectorized to compute all Nt eigenvalues
[x,iter]=eigBNewton(Nt,tol); 

D=1i*x/dt;

V=(1i.^((0:Nt-1)')).*ChebyU((0:Nt-1)',x.'); %vectorized matrix
%normV=vecnorm(V,2); V=V./normV;%normalize V
iV=(ChebyVand2WangNew(x)).*(1i.^((0:Nt-1)'))'; %O(nt^2) fast inversion of U
%iV=normV'.*iV; 
end

function [x,iter,b0]=eigBNewton(nt,tol,b0)
%use initial guess b0 from coarse mesh can speed up convergence
%O(nt) algorithm for finding all eigenvalues, vectorized
p=@(n,theta) sin(n*theta)-1i*cos(n*theta).*sin(theta); % polynomial
pp=@(n,theta) n*cos(n*theta)-1i*cos(n*theta).*cos(theta)+1i*n*sin(n*theta).*sin(theta);
pg=@(n,theta) theta-p(n,theta)./pp(n,theta);

maxit=100; %tol=1e-8;

k=(1:(nt)/2)';
if(nargin<3||isempty(b0))
    theta20=0.5*k*pi*(1/nt+1/(nt+1))+1i/nt;
end
for iter=1:maxit
    theta2=pg(nt,theta20);
    if(norm(theta2-theta20)<tol*norm(theta2)) break; end
    theta20=theta2;
end

theta2(nt+1-k)=pi-conj(theta2(k)); 

x=cos(theta2);
%x2=[x2(1:end/2).'; x2(end:-1:end/2+1).'];
%x=x2(:);
end
 
function W=ChebyVand2WangNew(x)
    n=length(x); e=ones(n,1); %ys=cos((1:n)*pi/(n+1));    
    S=spdiags([-e 2*e -e],[-2 0 2],n,n); S(1,1)=3;S(end,end)=3;
    v=zeros(n,1); v(end-1:end)=[1i;2];    
    b=S\v; %Step 1: O(n) fast direct solver 
    %b=pentsolve(S,v); %fast O(n) Thomas algorithm
    %dS=2*e;dS(1)=3;dS(end)=3; b=Thomas5(dS,-e,v);%slightly slower
    p_n_x=(-n*ChebyT(n,x)+x.*ChebyU(n-1,x))./(1-x.^2)-1i*n*ChebyU(n-1,x);
    A=zeros(n,n);  %Step 2: O(n^2) vectorized recursion 
    %using linear system for solving A, slightly more accurate but slower
    S0=spdiags([e 0*e e],[-1 0 1],n,n); I=speye(n);
    for j=1:n 
        A(j,:)=((S0-2*x(j)*I)\(2*b/p_n_x(j))).'; 
        %A(j,:)=Thomas(-e,2*x(j)*e, -e, (-2/p_n_x(j))*b).';%slightly slower
    end  
    W=A*S/2;%Step 3: O(n^2) dense matrix*sparse matrix 
end

function y=ChebyT(n,x)
y= cos(n*acos(x));
end

function y=ChebyU(n,x)
y= sin((n+1).*acos(x))./sin(acos(x));
end



