function [f]=merger_focs(p_to_solve)
global cdid ns mc deltajt0 market_no...
    x2 own_dummy theta1_repl theta2_est price theti thetj alphai_merger

% new price
price_new                  = price;
price_new(cdid==market_no) = p_to_solve;

% new delta
x22         = x2;
x22(:,2)    = price_new;
deltajt_new = deltajt0-price*theta1_repl(1)+price_new*theta1_repl(1);

theta2w     = full(sparse(theti,thetj,theta2_est));

% new exp(mu)
expmu_new   = exp(mufunc(x22,theta2w));

% new market shares
sjt   = mktsh(exp(deltajt_new) ,expmu_new);
sijt  = ind_sh(exp(deltajt_new),expmu_new);

% keep market shares relevant for the market under consideration
sjt   =   sjt(cdid==market_no,:);
sijt  =  sijt(cdid==market_no,:);

% keep only relevnat alphais
alphai_m = alphai_merger(cdid==market_no,:);
pjt      = p_to_solve;
deriv    = zeros(size(pjt,1),size(pjt,1));

% derivative matrix
for j=1:size(pjt,1)    
    deriv(j,j)=(1/ns)*sum(alphai_m(j,:).*sijt(j,:).*(ones(1,ns)-sijt(j,:)));        
    for k=1:size(pjt,1)
        if k~=j
            deriv(j,k)=-(1/ns)*sum(alphai_m(j,:).*sijt(j,:).*(sijt(k,:)));            
        end
    end
end

% FOCs
omega =(deriv).*(own_dummy*own_dummy');
f     = (p_to_solve-mc(cdid==market_no))+inv(omega')*sjt;






