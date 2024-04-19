function [U,P1,P2] = JSCGCCA(X1,X2,s1,s2,d)
 
    [n,d1]  = size(X1);
    [~,d2]  = size(X2);
    maxit   = 200;
    tol     = 1e-3;   
    P1      = zeros(d1,d);
    P2      = zeros(d2,d);
    U       = zeros(n,d);    
    
    Fnorm   = @(var)norm(var,'fro')^2;
    func1   = @(U0,P,T,J)CCA(X1,U0,P,T,J,d);
    func2   = @(U0,P,T,J)CCA(X2,U0,P,T,J,d);

for iter = 1: maxit
    
    %  update U
    U_old     = U;
    X         = X1 * P1 + X2 * P2;
    [Q,~,V]   = svd(X,'econ');
    U         = Q*V';  
    
    %  update P1 and P2
    P1_old    = P1;
    P2_old    = P2; 
    pars.P0   = P1;
    P1        = NHTP(@(P,T,J)func1(U,P,T,J),d1,d,s1,pars); 
    pars.P0   = P2;
    P2        = NHTP(@(P,T,J)func2(U,P,T,J),d2,d,s2,pars);
 
    % halting condition
    errU  = Fnorm(U_old - U)/(1+Fnorm(U_old));
    errP1 = Fnorm(P1 - P1_old)/(1+Fnorm(P1));
    errP2 = Fnorm(P2 - P2_old)/(1+Fnorm(P2));
    Obj   = norm(U - X1*P1,'fro') + norm(U - X2*P2,'fro');
%     fprintf('Iter = %3d   RelU = %4.2e   RelP1 = %4.2e   RelP2 = %4.2e  Obj = %4.2e  Time = %4.2f\n',...
%             iter,errU, errP1,errP2, Obj, toc(t0));
    if max([errU,errP1,errP2]) <tol; break; end
    
end

end
  