function [focs]=bert_focs(price_to_solve_for)

global price cdid mkt_no x2 vfull ns mc deltajt0 theta1 theta2 own_dummy_post

mvalold_opt             = deltajt;
price_opt               = price;

price_opt(cdid==mkt_no) = price_to_solve_for;
price_coeff             = theta1(1);

% Update mean utility for price
mvalold_opt             = mvalold_opt - price*price_coeff + price_opt*price_coeff;

theta2w=zeros(size(theta2,1),size(theta2,1));
theta2w(:,1)=theta2;

x2_temp=x2;
x2_temp(:,1)=price_opt;

s_ijt_opt               = ind_sh(exp(mvalold_opt),exp(mufunc(x2_temp,theta2w)));
s_jt_opt                = mktsh(exp(mvalold_opt),exp(mufunc(x2_temp,theta2w)));

% Limit to only relevant market
s_jt_opt                =   s_jt_opt(cdid==mkt_no,:);
s_ijt_opt               =  s_ijt_opt(cdid==mkt_no,:);
vfull_opt               =      vfull(cdid==mkt_no,:);
price_opt               =  price_opt(cdid==mkt_no,:);
mc_bert_opt             =         mc(cdid==mkt_no,:);

%2217 x ns matrix of price random coefficients
alphai = vfull(:,1:ns).*repmat(theta2(1),size(vfull,1),ns)+...
    repmat(theta1(1),size(vfull,1),ns);

%Derivative matrix and the implied margins for the relevant market
alpha_i = alphai(cdid==mkt_no,:);
s_ijt   = s_ijt_opt;

deriv=[];
for j=1:size(s_ijt,1)
    deriv(j,j)= (1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
    for k=1:size(s_ijt,1)
        if k~=j
            deriv(j,k)= -(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
        end
    end
end

% Calculate focs
om   = deriv.*(own_dummy_post(cdid==mkt_no,:)*own_dummy_post(cdid==mkt_no,:)');
focs = (price_opt-mc_bert_opt)+ inv(om')*s_jt_opt;








