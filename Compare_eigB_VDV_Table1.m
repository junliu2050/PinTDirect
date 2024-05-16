%% Compare eig and our fast Newton method for compute eigdecomposition 
% produce the results reported in Table 1
set(0, 'defaultaxesfontsize',20,'defaultaxeslinewidth',1,...
    'defaultlinelinewidth',2,'defaultpatchlinewidth',1.5,...
    'defaulttextfontsize',20,'defaulttextInterpreter','latex');

clear 
t1=[];t2=[]; t3=[]; Ntlist=2.^(5:10);
err1=[];err2=[]; err3=[];tol=1e-10;
fprintf('Nt\t CPU(eig)\t eig-res\t  eig-It  CPU(fast)\t fast-res\t eig-diff\n');
T=2;
for Nt=Ntlist  
    dt=T/Nt;
    tic
    e=ones(Nt,1); %E=eye(Nt);
    B=spdiags([-e/2 e/2],[-1 1],Nt,Nt)/dt;
    B(end,end-1)=-1/dt; B(end,end)=1/dt;%last row
    [V1,eigB1]=eig(full(B),'vector');  
    err2=[err2;norm(B-(V1.*(eigB1.'))/V1,'fro')/norm(B,'fro')];
    t1=[t1;toc];    
     
    tic 
    [V2,eigB2,iV2,iter]=fasteigB(Nt,dt,tol);  
    err3=[err3;norm(B-(V2.*(eigB2.'))*iV2,'fro')/norm(B,'fro')];
    t2=[t2;toc]; 
    
    
    eigB1=cplxpair(eigB1,tol);
    eigB2=cplxpair(eigB2,tol);
    
    err1=[err1;norm(eigB1-eigB2,'fro')./norm(eigB1,'fro')];
    fprintf(' %d&\t %1.3f &\t  %1.2e&\t %d &\t %1.3f &\t %1.2e&\t %1.2e \\\\\n',...
        Nt,t1(end),err2(end),iter,t2(end),err3(end),err1(end));
end

 
 