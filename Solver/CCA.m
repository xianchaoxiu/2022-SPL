function [out1,out2] = CCA(Xv,U,Pv,Tv,Jv,d)

if isempty(Tv) && isempty(Jv)
    Tx    = find( sum(abs(Pv),2)>0 );
    XvPvU = Xv(:,Tx)*Pv(Tx,:)-U; 
    out1  = norm(XvPvU,'fro')^2/2; % objective
    if nargout == 2    
      out2 = Xv'*XvPvU;
    end
else
    XvTv  = Xv(:,Tv); 
    XvTvt = XvTv';
    out1  = kron(eye(d),XvTvt* XvTv); %Hessian TvTv
    if nargout == 2
       out2 = @(D) (reshape(XvTvt*(Xv(:,Jv)*D),[],1) ); 
    end    
end 

end

