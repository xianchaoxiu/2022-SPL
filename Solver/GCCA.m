function [U,P1,P2] = GCCA(X1,X2,d)

    [n,d1]  = size(X1);
    [~,d2]  = size(X2);
    p       = 0.001;
    maxit   = 200;
    tol     = 1e-3;
    P1      = zeros(d1,d);
    P2      = zeros(d2,d);
    U       = zeros(n,d);
    X1T     = X1';
    X11     = inv(X1T * X1+p*eye(d1)) * X1T; 
    X2T     = X2';
    X22     = inv(X2T * X2+p*eye(d2)) * X2T; 
    
for iter = 1: maxit
    
        %  update U
        U_ori     = U;
        X         = X1 * P1 + X2 * P2; 
        [Q,~,V]   = svd(X,'econ');
        U         = Q*V';  
        
        %  update P1 and P2
        P1_ori   = P1;
        P2_ori   = P2;
        P1       = X11 * U; 
        P2       = X22 * U;   
       
        %  keep record 
       Obj(iter+1)   = norm(U - X1*P1,'fro')+norm(U - X2*P2,'fro');
       RelP1(iter+1) = norm(P1 - P1_ori,'fro')/(1+norm(P1,'fro'));
       RelP2(iter+1) = norm(P2 - P2_ori,'fro')/(1+norm(P2,'fro'));
       RelU(iter+1) = norm(U - U_ori,'fro')/(1+norm(U,'fro'));
%         fprintf('It:%2d, RelP1:%4.2e, RelP2:%4.2e, RelU:%4.2e, Obj:%4.2e,\n',...
%            iter,RelP1(iter+1),RelP2(iter+1),RelU(iter+1),Obj(iter+1));
    if max([RelP1(iter+1),RelP2(iter+1),RelU(iter+1)])<tol; break; end
    
end
 
end