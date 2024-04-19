function [U,P1,P2] = JSGCCA(X1,X2,r1,r2,d)

    [n,d1]  = size(X1);
    [~,d2]  = size(X2);
    maxit   = 200;
    tol     = 1e-3;
    beta    = 1;             %惩罚参数
    X1T     = X1';
    XTX1    = X1T * X1; 
    X2T     = X2';
    XTX2    = X2T * X2; 
    invX1   = inv(XTX1+beta*eye(d1)); 
    invX2   = inv(XTX2+beta*eye(d2));    
    P1      = zeros(d1,d);
    P2      = zeros(d2,d);
    Q1      = zeros(d1,d);
    Q2      = zeros(d2,d);
    Z1      = zeros(d1,d);
    Z2      = zeros(d2,d);
    U       = zeros(n,d);    
    
for iter = 1: maxit
    
        %  update U
        U_ori     = U;
        X         = X1 * P1 + X2 * P2; 
        [Q,~,V]   = svd(X,'econ');
        U         = Q*V';  
       
        %  update P1 and P2
        
        P1_ori   = P1;
        P2_ori   = P2;
        P1       = invX1*(X1T*U+Z1+beta*Q1);
        P2       = invX2*(X2T*U+Z2+beta*Q2);
         
        %  update Q1 and Q2
        Q1      = prox_l21(P1-Z1/beta,r1/beta);
        Q2      = prox_l21(P2-Z2/beta,r2/beta);        

        %  udate Z1 and Z2
        Z1     = Z1 - beta*(P1-Q1); 
        Z2     = Z2 - beta*(P2-Q2); 
        
        
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
  