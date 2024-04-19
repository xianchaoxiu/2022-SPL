function [Pv,out] = NHTP(func,dv,d,sv,pars)
%This code aims at solving the sparsity constrained optimization
%     min_Pv \|XvPv-U\|_F^2), s.t. ||Pv||_0<=sv,
% where Pv\in\R^{dv\times d}, 
%       sv<<dv is an given integer
%       ||Pv||_2,0 is the number of nonzero rows of Pv
%     where  U\in R^{n x d}  
%       Xv \in R^{n x dv}
%       Pv \in R^{dv x d}
% Inputs: 
%        func  -- a function handle computing                  (required)
%                   (objective, gradient, sub-hessian)
%       

%        sv    -- the sparsity level, an integer in (0,dv)      (required)
%        pars  -- a structure
%                 pars.maxit = 1000        (default)
%                 pars.tol   = 1e-8        (default)
%                 pars.alpha = 2           (default)
%                 pars.Pv0    = zeros(dv,d)  (default)
% Outputs: 
%       Pv        -- the sparse solution
%       out.f     -- f(out.x)     
%       out.time  -- computational time  
%       out.iter  -- number of iterations


t0    = tic;
if nargin < 5; pars = []; end
if isfield(pars,'alpha'); alpha = pars.alpha; else; alpha = 1e-4;       end
if isfield(pars,'maxit'); maxit = pars.maxit; else; maxit = 1e3;        end
if isfield(pars,'P0');    Pv0   = pars.P0;    else; Pv0 = zeros(dv,d);  end

tol    = 1e-6;
sigma  = 1e-8;
gamma  = 1e-10;
beta   = 0.5; 
Pv     = Pv0;
 
[f0,gv] = func(Pv,[],[]);%compute gv use CCA(out2)
[~,Tv]  = maxk(sum((Pv-alpha* gv).^2,2),sv);%选取前sv行的指标集
Fnorm  = @(var)norm(var,'fro')^2;
Tv0     = Tv;
f       = f0;
% fprintf('\n Start to run the solver -- NHTP\n')
% fprintf(' ----------------------------------------------------\n')
% fprintf(' Iter         Error          Objective         Time\n')
% fprintf(' ----------------------------------------------------\n')

% main body
for iter  = 1 : maxit
   
    gvTv  = gv(Tv,:);
    PvTv  = Pv(Tv,:);
    temp1 = max(0,Fnorm(Pv)-Fnorm(PvTv)); 
    error = Fnorm(gvTv)+temp1;  
%     fprintf(' %3d        %1.2e         %1.4f         %.2fsec\n',...
%                 iter, error, f0, toc(t0)) 
    if error < tol || (iter>1 && abs(f-f0)<1e-8*f0); break; end
    
    % find the direction
    Jv      = setdiff(Tv0,Tv);     
    vecgvTv = reshape(gvTv,[],1);
    if isempty(Jv)
       vecDvTv = pinv(func(Pv,Tv,Jv))*(-vecgvTv );
     else
       [H1,H2] = func(Pv,Tv,Jv);
       vecDvTv = pinv(H1)*( H2(Pv(Jv,:))-vecgvTv );
     end
    DvTv    = reshape(vecDvTv,[sv,d]);
    marker  = (trace(DvTv'*gvTv) >= -gamma*(Fnorm(DvTv)+temp1) + temp1/4/alpha );   
    if norm(DvTv,'fro')^2 > 1e16 || marker
       DvTv = - gvTv; 
    end   
    Dv      = -Pv;
    Dv(Tv,:)= DvTv;
    
    % Armijio line search
    alphak  = 1;
    temp2   = sigma*trace(gv'*Dv); 
    Pv      = zeros(dv,d);
    f0      = f;
    for j   = 1 : 8
        Pv(Tv,:) = Pv0(Tv,:) + alphak*DvTv;  
        f        = func(Pv,[],[]);   
        if f    <= f0 + alphak*temp2; break; end
        alphak   = alphak*beta;  
    end
    
    Tv0      = Tv; 
    Pv0      = Pv;
    [f,gv]   = func(Pv,[],[]);
    [~,Tv]   = maxk(sum((Pv - alpha * gv).^2,2),sv);
    
    if mod(iter,5)==0; alpha = max(alpha/1.5,1e-5); end
  
end

% fprintf(' ----------------------------------------------------\n')
 
out.f    = f;
out.iter = iter;
out.time = toc(t0);
Pv       =  Pv;

end

