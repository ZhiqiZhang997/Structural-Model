%**************************************************************************
% To accompany Knittel and Metaxoglou
%**************************************************************************

clear all
close all
close hidden
warning off all
format long
clc

km       = maxNumCompThreads(1);
mydims   = 5;

%**************************************************************************
% Define globals
%**************************************************************************
global ns x1 x2 s_jt vfull dfull theta1 theti thetj...
    cdid cdindex IV1 invA1 nbrn...
    mvalold_logit v mvalold oldt2 gmmresid...
    mymaxfunevals mvalolds mvalold0 mvalold00 mvalolds2...
    fvals_track fcnevals ppp mval_track deltajt0 price_pre mkt_no theta2_est...
    theta1_repl mc own_dummy_post

load BLP_data cdid cdindex share outshr price firmid id const hpwt air mpd space mpg trend...
    product model_id own_dummies

load BLP_data_str model_name

% *************************************************************************
% Define paths for input and output
% Controls for full simulation and calculation checks
% *************************************************************************
flag_check      = 0; % if equal to 1, post and 1pct calculations are identical
flag_best       = 1; % if equal to 1, limit analysis to params for min fval
flag_full_simul = 0; % if equal to 1, full merger simulation analysis
p_inc           = 1.01;
code_path       = pwd;
results_path    = [code_path,'\Optimization results\' ];
merger_path     = [code_path,'\Merger results\'   ];
logs_path       = merger_path;
add_path        = [code_path,'\Optimization routines\'];
outfile         = [logs_path,'elast_boot.txt'    ];
fid             = fopen(outfile,'w'); fclose(fid);
tic
diary(outfile)
addpath(add_path);

% *************************************************************************
% Demand instruments
% *************************************************************************
sum_other=[];
sum_rival=[];

X = [const,(hpwt),air,(mpg),space];
for i=1:size(id,1)
    other_ind      = (firmid==firmid(i)  & cdid==cdid(i) & id~=id(i));
    rival_ind      = (firmid~=firmid(i)  & cdid==cdid(i));
    total_ind      = (cdid==cdid(i));
    
    sum_other(i,:) = sum(X(other_ind==1,:));
    sum_rival(i,:) = sum(X(rival_ind==1,:));
    sum_total(i,:) = sum(X(total_ind==1,:));
    
    if size(X(other_ind==1,:),1)==0
        sum_other(i,:)=zeros(1,size(X,2));
    end
    
    if size(X(other_ind==1,:),1)==1
        sum_other(i,:)=X(other_ind==1,:);
    end
end

IV1=[X,sum_other,sum_rival];

%Load N(0,I) drwas
load v

% %ownership structure matrix
own=own_dummies; clear owndummies

%variables in demand without random coefficients
x1=[price,X]; clear X

%variables in demand with random coefficients
x2=x1(:,1:size(x1,2)-1);

%#of indivduals and markets
ns =  size(v,2)/5;
nmkt = 20;

%#of brands per market
nbrn=zeros(nmkt,1);
nbrn(1)=cdindex(1);
for i=2:max(nmkt)
    nbrn(i)=sum(cdid==i);
end

%demographics and N(0,I) unobservables
demogr=zeros(size(v));
vfull = v(cdid,:);
dfull = demogr(cdid,:);

% *************************************************************************
% Logit regressions
% *************************************************************************
theta2=ones(5,1);
theta2w=zeros(5,5);
theta2w(:,1)=theta2;
[theti, thetj, theta2]=find(theta2w);

invA1 = inv(IV1'*IV1);

s_jt=share;
y = log(s_jt) - log(outshr);

mid = x1'*IV1*invA1*IV1';
theta1 = inv(mid*x1)*mid*y;
mvalold_logit = x1*theta1;

n=size(x1,1);
k=size(x1,2);

ESS=y'*y-2*theta1'*x1'*y+theta1'*x1'*x1*theta1;
s2=ESS/(n-k);
A=(x1'*(IV1*inv(IV1'*IV1)*IV1')*x1);
se=sqrt(diag(s2*inv(A)));

optrout_no_all  = [];
perturbs_all    = [];
counts_all      = [];
exit_infos_all  = [];
tocs_all        = [];

%**************************************************************************
% Loop over optimization routines
%**************************************************************************
optrouts = [1,2,3,4,5,7,8,9,10,16]';
optrouts = [2]';

for km_optrout=1:1:size(optrouts,1)
    
    optrout = optrouts(km_optrout,1);
    
    if optrout~=6
        
        cd(results_path)
        
        if optrout<=9
            matfile           = ['blp_0',num2str(optrout),'_data_optim'];
            elast_boot_file   = [merger_path,'blp_elast_boot_0',num2str(optrout),'.txt'];
        else
            matfile           = ['blp_',num2str(optrout),'_data_optim'];
            elast_boot_file   = [merger_path,'blp_elast_boot_',num2str(optrout),'.txt'];
        end
        
        load (matfile, 'theta1s', 'theta2s','deltas','fvals','perturbs2',...
            'gmmresids','exit_infos','counts2','tocs');
        
        cd(code_path)
        
        optrout_no      = repmat(optrout,size(fvals,1),1);
        optrout_no_all  = [optrout_no_all;optrout_no];
        perturbs_all    = [perturbs_all;perturbs2];
        counts_all      = [counts_all;counts2];
        exit_infos_all  = [exit_infos_all;exit_infos];
        tocs_all        = [tocs_all;tocs];
        
        %******************************************************************
        % Limit analysis to starting values produced min fval
        %******************************************************************
        if flag_best==1
            [min_fval,min_fval_ind]=min(fvals);
            theta1s   = theta1s(min_fval_ind,:);
            theta2s   = theta2s(min_fval_ind,:);
            deltas    = deltas(:,min_fval_ind);
            gmmresids = gmmresids(:,min_fval_ind);
            fvals     = fvals(min_fval_ind,:);
            perturbs2 = perturbs2(min_fval_ind,:);
        end
        
        merger_results = [];
        
        % *****************************************************************
        % Loop over starting values
        % *****************************************************************
        for jj_optrout=1:1:size(fvals,1)
            
            % Track information collected in main.m
            perturb      = perturbs2(jj_optrout  ,:) ;
            theta1_est   = theta1s(jj_optrout    ,:)';
            theta2_est   = theta2s(jj_optrout    ,:)';
            delta_est    = deltas(:,jj_optrout   ,:) ;
            fval_est     = fvals(jj_optrout      ,:) ;
            min_index    = fval_est==min(min(fvals)) ;
            gmmresid_est = gmmresids(:,jj_optrout,:) ;
            
            temp1        = gmmresid_est'*IV1;
            fval_est2    = temp1*invA1*temp1';
            
            % Initialize NFP with exp(y);
            mvalold      = exp(y);
            oldt2        = zeros(size(theta2_est));
            
            % Account for DFS v. Nevo NFP
            % Both gmmobj2 and gmmobj2_dfs use 1e-16
            if (optrout>=13 && optrout<=15)
                [fval_repl,grad_repl] = gmmobj2_dfs(theta2_est);
            else
                [fval_repl,grad_repl] = gmmobj2(theta2_est);
            end
            
            % Gradient norm inf
            grad_norm_inf = max(abs(grad_repl));
            
            % Make sure that theta2,theta1,
            % delta, gmmresid and fval are internally consistent
            vcov           = var_cov(theta2_est);
            se             = full(sqrt(diag(vcov)));
            std_errors     = se';
            
            std_errors_theta1 = se(1:6);
            std_errors_theta2 = se(7:11);
            
            theta1_repl    = theta1;
            gmmresid_repl  = gmmresid;
            
            temp1          = gmmresid_repl'*IV1;
            fval_repl2     = temp1*invA1*temp1';
            delta_repl     = log(mvalolds(:,2));
            delta_repl2    = gmmresid_repl+x1*theta1;
            
            temp1          = x1'*IV1;
            temp2          = delta_repl'*IV1;
            theta1_repl2   = inv(temp1*invA1*temp1')*temp1*invA1*temp2';
            gmmresid_repl2 = delta_repl - x1*theta1;
            
            %**************************************************************
            % Use theta2_est and theta1_repl
            %**************************************************************
            theta2 = theta2_est;
            theta1 = theta1_repl;
            
            theta2w     = zeros(mydims,mydims);
            theta2w(:,1)= theta2;
            expmu       = exp(mufunc(x2,theta2w));
            
            % Replicate market shares
            s_jt_repl  =  mktsh(exp(delta_repl),expmu);
            s_ijt_repl = ind_sh(exp(delta_repl),expmu);
            
            fprintf('===============================\n')
            fprintf('        Merger Diagnostics     \n');
            fprintf('===============================\n')
            fprintf('optrout          : %12i\n'   ,optrout                              );
            fprintf('perturb          : %12i\n'   ,perturb                              );
            fprintf('fval est         : %12.4f\n',fval_est                              );
            fprintf('fval est check   : %12.4f\n',fval_est2                             );
            fprintf('fval repl        : %12.4f\n',fval_repl                             );
            fprintf('fval repl check  : %12.4f\n',fval_repl2                            );
            fprintf('theta1 check     : %12.4f\n',max(abs(theta1_repl-theta1_est      )));
            fprintf('delta check      : %12.4f\n',max(abs(delta_repl-delta_repl2      )));
            fprintf('gmmresid check   : %12.4f\n',max(abs(gmmresid_repl-gmmresid_repl2)));
            fprintf('s_jt check       : %12.4f\n',max(abs(s_jt_repl-s_jt              )));
            fprintf('grad_norm_inf    : %12.4f\n',grad_norm_inf                         );
            fprintf('===============================\n')
            fprintf('fval check     2 : %12.4f\n',max(abs(fval_repl-fval_est          )));
            fprintf('theta1 check   2 : %12.4f\n',max(abs(theta1_repl-theta1_est      )));
            fprintf('delta check    2 : %12.4f\n',max(abs(delta_repl-delta_est        )));
            fprintf('gmmresid check 2 : %12.4f\n',max(abs(gmmresid_repl-gmmresid_est  )));
            fprintf('===============================\n')
            
            % Market shares and mean utility levels pre-merger
            sijt_pre       = s_ijt_repl;
            sjt_pre        = s_jt_repl;
            deltajt_pre    = delta_repl;
            
            % Price random coefficient
            vfull1 = vfull(:,1:ns);
            alpha_i=[];
            for i=1:size(vfull1,1)
                alpha_i(i,:) = vfull1(i,:).*(kron(theta2(1),ones(1,ns)))+...
                    (kron(theta1(1),ones(1,ns)));
            end
            
            alphai    = alpha_i;
            
            % Derive matrices of price derivatives and elasticities
            deriv_all = zeros(max(nbrn),max(nbrn),nmkt);
            elast_all = zeros(max(nbrn),max(nbrn),nmkt);
            
            % Loop over markets
            for i=1:nmkt
                
                % Use replicated market shares
                ind     = cdid==i;
                pjt     = price(ind==1,1);
                sjt     = sjt_pre(ind==1,1);
                alpha_i = alphai(ind==1,:);
                s_ijt   = sijt_pre(ind==1,:);
                
                elast = zeros(size(pjt,1),size(pjt,1));
                deriv = zeros(size(pjt,1),size(pjt,1));
                
                % Loop over products in the same market
                for j=1:size(pjt,1)
                    deriv(j,j) =                  (1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
                    elast(j,j) = (pjt(j)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
                    
                    for k=1:size(pjt,1)
                        
                        if k~=j
                            elast(j,k) = -(pjt(k)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                            deriv(j,k) =                  -(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                        end
                        
                    end
                    
                end
                
                elast_all(1:size(elast,1),1:size(elast,2),i) = elast;
                deriv_all(1:size(deriv,1),1:size(deriv,2),i) = deriv;
                
            end
            
            % Track elasticities
            elast_own   = [];
            elast_cross = [];
            for i=1:nmkt
                temp        = diag(elast_all(1:nbrn(i),1:nbrn(i),i));
                elast_own   = [elast_own;temp];
                elast_cross = [elast_cross;elast_all(1:nbrn(i),:,i)];
            end
            
            fprintf('mean elast_own   : %12.4f\n',mean(elast_own)  );
            fprintf('std elast_own    : %12.4f\n',std(elast_own)   );
            fprintf('min elast_own    : %12.4f\n',min(elast_own)   );
            fprintf('max elast_own    : %12.4f\n',max(elast_own)   );
            fprintf('median elast_own : %12.4f\n',median(elast_own));
            
        end
        
        cd(merger_path);
        
        cd(code_path);
        
    end
end

%**************************************************************************
% bootstrap the own elasticity
%**************************************************************************
randn('seed',19751208);
draws           = 10^4;
vcov_draws      = vcov(7:11,7:11);
theta2_draws_01 = randn(size(theta2_est,1),draws);
theta2_draws    = repmat(theta2_est',draws,1)+(chol(vcov_draws)'*theta2_draws_01)';

s_jt_draws      = NaN*ones(size(price,1),draws);
delta_draws     = NaN*ones(size(price,1),draws);
elast_own_draws = NaN*ones(size(price,1),draws);
elast_boot_results =[];

for d=1:draws
    fprintf('iteration:%5i\n',d);
    [fval,grad] = gmmobj2(theta2_draws(d,:)');
    delta_draws(:,d) = log(mvalold);
    
    theta2      = theta2_draws(d,:);
    theta2w     = zeros(mydims,mydims);
    theta2w(:,1)= theta2;
    expmu       = exp(mufunc(x2,theta2w));
    
    % Replicate market shares
    s_jt_draws(:,d)  =  mktsh(exp(delta_draws(:,d)),expmu);
    
    % Price random coefficient
    vfull1 = vfull(:,1:ns);
    alpha_i=[];
    for i=1:size(vfull1,1)
        alpha_i(i,:) = vfull1(i,:).*(kron(theta2(1),ones(1,ns)))+...
            (kron(theta1(1),ones(1,ns)));
    end
    
    alphai    = alpha_i;
    
    % Derive matrices of price derivatives and elasticities
    deriv_all = zeros(max(nbrn),max(nbrn),nmkt);
    elast_all = zeros(max(nbrn),max(nbrn),nmkt);
    
    % Loop over markets
    for i=1:nmkt
        
        % Use replicated market shares
        ind     = cdid==i;
        pjt     = price(ind==1,1);
        sjt     = sjt_pre(ind==1,1);
        alpha_i = alphai(ind==1,:);
        s_ijt   = sijt_pre(ind==1,:);
        
        elast = zeros(size(pjt,1),size(pjt,1));
        deriv = zeros(size(pjt,1),size(pjt,1));
        
        % Loop over products in the same market
        for j=1:size(pjt,1)
            deriv(j,j) =                  (1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
            elast(j,j) = (pjt(j)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
            
            for k=1:size(pjt,1)
                
                if k~=j
                    elast(j,k) = -(pjt(k)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                    deriv(j,k) =                  -(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                end
                
            end
            
        end
        
        elast_all(1:size(elast,1),1:size(elast,2),i) = elast;
        deriv_all(1:size(deriv,1),1:size(deriv,2),i) = deriv;
        
    end
    
    % Track elasticities
    elast_own   = [];
    elast_cross = [];
    for i=1:nmkt
        temp        = diag(elast_all(1:nbrn(i),1:nbrn(i),i));
        elast_own   = [elast_own;temp];
        elast_cross = [elast_cross;elast_all(1:nbrn(i),:,i)];
    end
    
    elast_own_draws(:,d)=elast_own;
    
    
    market=cdid;
    elast_boot_results_temp = [repmat(d,size(market,1),1),id,elast_own_draws(:,d)];
    ind1                    = id==268;
    ind2                    = id==3772;
    ind                     = ind1+ind2;
    elast_boot_results      = [elast_boot_results;elast_boot_results_temp(ind==1,:)];
    
end
diary off

head_xls = cellstr({'draw','id','elast_own'});
data_xls = num2cell(elast_boot_results);
data_xls = [head_xls;data_xls];
xlswrite([merger_path,'elast_boot.xlsx'],data_xls);
