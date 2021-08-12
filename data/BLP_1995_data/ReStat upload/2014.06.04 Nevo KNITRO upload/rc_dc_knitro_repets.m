clear all
close all
format long
clc

optrout      =4;
code_path    =pwd;
results_path =[code_path,'/optimization results/'];
logs_path    =[code_path,'/optimization logs/'];
add_path     =[code_path,'/optimization routines/'];
addpath(add_path);

global invA ns x1 x2 s_jt IV vfull dfull theta1 theti thetj cdid cdindex

% load data. see description in readme.txt
load ps2
load iv
IV = [iv(:,2:21) x1(:,2:25)];
clear iv

ns = 20;       % number of simulated "indviduals" per market %
nmkt = 94;     % number of markets = (# of cities)*(# of quarters)  %
nbrn = 24;     % number of brands per market. if the numebr differs by market this requires some "accounting" vector %
% this vector relates each observation to the market it is in %
cdid = kron([1:nmkt]',ones(nbrn,1));
% this vector provides for each index the of the last observation %
% in the data used here all brands appear in all markets. if this %
% is not the case the two vectors, cdid and cdindex, have to be   %
% created in a different fashion but the rest of the program works fine.%
cdindex = [nbrn:nbrn:nbrn*nmkt]';


% starting values. zero elements in the following matrix correspond to %
% coeff that will not be max over,i.e are fixed at zero. %
theta2w=    [0.3302   5.4819         0    0.2037         0;
    2.4526  15.8935    -1.2000        0    2.6342;
    0.0163  -0.2506         0    0.0511         0;
    0.2441   1.2650         0   -0.8091         0];

% create a vector of the non-zero elements in the above matrix, and the %
% corresponding row and column indices. this facilitates passing values %
% to the functions below. %
[theti, thetj, theta2]=find(theta2w);

horz=['    mean         sigma        income      income^2      age       child'];
vert=['constant  ';
    'price     ';
    'sugar     ';
    'mushy     '];

% create weight matrix
invA = inv([IV'*IV]);

% Logit results and save the mean utility as initial values for the search below

% compute the outside good market share by market
temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
outshr = 1.0 - sum1(cdid,:);

y = log(s_jt) - log(outshr);
mid = x1'*IV*invA*IV';
t = inv(mid*x1)*mid*y;
mvalold = x1*t;
oldt2 = zeros(size(theta2));
ESS = y'*y-2*t'*x1'*y+t'*x1'*x1*t;
s2  = ESS/(size(x1,1)-size(x1,2));

mvalold = exp(mvalold);
mvalold0 = mvalold;

save mvalold mvalold oldt2
clear mid y outshr t oldt2 mvalold temp sum1

vfull = v(cdid,:);
dfull = demogr(cdid,:);

counts2    = [];
perturbs2  = [];
gmmresids  = [];
deltas     = [];
fvals      = [];
exit_infos = [];
theta2s    = [];
theta1s    = [];
std_errors = [];
coeffs2    = [];
tocs       = [];

if optrout==1
    ktropts = optimset(...
        'Algorithm','interior-point',...
        'Hessian','bfgs',...
        'Display','iter',...
        'GradObj','on',...
        'TolCon',1e-3,...
        'TolFun',1e-3,...
        'TolX'  ,1e-3);
end
if optrout==2
    ktropts = optimset(...
        'Algorithm','interior-point',...
        'Hessian','bfgs',...
        'Display','iter',...
        'GradObj','on',...
        'TolCon',1e-3,...
        'TolFun',1e-3,...
        'TolX'  ,1e-3);
end
if optrout==3
    ktropts = optimset(...
        'Algorithm','interior-point',...
        'Hessian','bfgs',...
        'Display','iter',...
        'GradObj','on',...
        'TolCon',1e-6,...
        'TolFun',1e-6,...
        'TolX'  ,1e-6);
end
if optrout==4
    ktropts = optimset(...
        'Algorithm','interior-point',...
        'Hessian','bfgs',...
        'Display','iter',...
        'GradObj','on',...
        'TolCon',1e-6,...
        'TolFun',1e-6,...
        'TolX'  ,1e-6);
end
if optrout==5
    ktropts = optimset(...
        'Algorithm','interior-point',...
        'Hessian','bfgs',...
        'Display','iter',...
        'GradObj','on',...
        'TolCon',1e-6,...
        'TolFun',1e-6,...
        'TolX'  ,1e-6);
end

matfile=[results_path,'knitro_',num2str(optrout),'.mat'];
logfile=[logs_path   ,'knitro_',num2str(optrout),'.log'];

diary(logfile);

for perturb=1:50
    tic    
    mylb = 0;
    myub = [];
    
    
    % starting values. zero elements in the following matrix correspond to %
    % coeff that will not be max over,i.e are fixed at zero. %
    theta2w=    [0.3302   5.4819         0    0.2037         0;
        2.4526  15.8935    -1.2000        0    2.6342;
        0.0163  -0.2506         0    0.0511         0;
        0.2441   1.2650         0   -0.8091         0];
    
    % create a vector of the non-zero elements in the above matrix, and the %
    % corresponding row and column indices. this facilitates passing values %
    % to the functions below. %
    [theti, thetj, theta2]=find(theta2w);
    
    % initialize theta2s
    randn('state',1000*perturb);
    theta2   = normrnd(0,1,size(theta2));
    theta2w0 = full(sparse(theti,thetj,theta2));
    
    % initialize deltas
    randn('state',1000*perturb);
    mvalold = exp(mvalold0+normrnd(0,sqrt(s2),size(x1,1),1));
    oldt2   = zeros(size(theta2));
    
    fprintf('============================================\n');
    fprintf(' KNITRO iteration: %5i\n',perturb);
    fprintf('============================================\n');
    
    if optrout==1
        [theta2,fval,exit_info,tmp,lambda] = ktrlink(@gmmobj,theta2,...
            [],[],[],[],mylb,myub,[],ktropts);
    end
    if optrout==2
        [theta2,fval,exit_info,tmp,lambda] = ktrlink(@gmmobj_dfs,theta2,...
            [],[],[],[],mylb,myub,[],ktropts);
    end
    if optrout==3
        [theta2,fval,exit_info,tmp,lambda] = ktrlink(@gmmobj_dfs,theta2,...
            [],[],[],[],mylb,myub,[],ktropts);
    end
    if optrout==4
        [theta2,fval,exit_info,tmp,lambda] = ktrlink(@gmmobj,theta2,...
            [],[],[],[],mylb,myub,[],ktropts);
    end
    if optrout==5
        [theta2,fval,exit_info,tmp,lambda] = ktrlink(@gmmobj_dfs,theta2,...
            [],[],[],[],[],[],[],ktropts);
    end
    
    
    counts2    = [counts2;tmp.funcCount];
    perturbs2  = [perturbs2;perturb];
    load gmmresid
    gmmresids  = [gmmresids,gmmresid];
    deltas     = [deltas,log(mvalold)];
    fvals      = [fvals;fval];
    exit_infos = [exit_infos;exit_info];
    
    theta2s    = [theta2s;theta2'];
    theta1s    = [theta1s;theta1'];
    
    fprintf('\nObj. function : \t');
    printm(fval);
    
    fprintf('\ntheta1        : \t');
    printm(theta1');
    fprintf('\ntheta2        : \t');
    printm(theta2');
    
    % standard errors
    vcov       = var_cov(theta2);
    se         = full(sqrt(diag(vcov)));
    std_errors = [std_errors;se'];
    
    toc_tmp=toc;
    comp_t = toc_tmp/60;
    
    tocs =[tocs;toc_tmp];
    
    % computing the s.e.
    vcov = var_cov(theta2);
    se = sqrt(diag(vcov));
    
    theta2w = full(sparse(theti,thetj,theta2));
    t = size(se,1) - size(theta2,1);
    se2w = full(sparse(theti,thetj,se(t+1:size(se,1))));
    
    % the MD estimates
    omega = inv(vcov(2:25,2:25));
    xmd = [x2(1:24,1) x2(1:24,3:4)];
    ymd = theta1(2:25);
    
    beta = inv(xmd'*omega*xmd)*xmd'*omega*ymd;
    resmd = ymd - xmd*beta;
    semd = sqrt(diag(inv(xmd'*omega*xmd)));
    mcoef = [beta(1); theta1(1); beta(2:3)];
    semcoef = [semd(1); se(1); semd];
    
    coeffs   = [mcoef theta2w];
    coeffs2  = [coeffs2;coeffs(:)'];
    
    Rsq = 1-((resmd-mean(resmd))'*(resmd-mean(resmd)))/((ymd-mean(ymd))'*(ymd-mean(ymd)));
    Rsq_G = 1-(resmd'*omega*resmd)/((ymd-mean(ymd))'*omega*(ymd-mean(ymd)));
    Chisq = size(id,1)*resmd'*omega*resmd;
    
    fprintf('\n\n');
    disp(horz)
    disp('  ')
    
    semcoef=full(semcoef);
    
    
    for i=1:size(theta2w,1)
        disp(vert(i,:))
        printm([mcoef(i) theta2w(i,:)])
        printm([semcoef(i) se2w(i,:)])
    end
    
    fprintf('\n\n');
    
    disp(['GMM objective                       :  ' num2str(fval)])
    disp(['MD R-squared                        :  ' num2str(Rsq)])
    disp(['MD weighted R-squared               :  ' num2str(Rsq_G)])
    disp(['# of objective function evaluations :  ' num2str(tmp.funcCount)])
    disp(['run time (minutes)                  :  ' num2str(comp_t)])
end

diary off

save (matfile, 'perturbs2', 'fvals', 'theta1s', 'theta2s','exit_infos',...
    'deltas' ,'gmmresids' ,'std_errors','counts2','tocs','coeffs2');