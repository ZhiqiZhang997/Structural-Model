%**************************************************************************
% To accompany Knittel and Metaxoglou
%**************************************************************************
clear all
close all
close hidden
warning off all
format long
clc
tic
km       = maxNumCompThreads(1);
mydims   = 5;

%**************************************************************************
% Define globals
%**************************************************************************
global ns x1 x2 s_jt vfull dfull theta1 theti thetj...
    cdid cdindex IV1 invA1 nbrn...
    mvalold_logit v mvalold oldt2 gmmresid...
    mymaxfunevals mvalolds mvalold0 mvalold00 mvalolds2...
    fvals_track fcnevals ppp mval_track

load BLP_data cdid cdindex share outshr price firmid id const hpwt air mpd space mpg trend...
    product model_id own_dummies

load BLP_data_str model_name

% *************************************************************************
% Define paths for codes, optimization results and logs
% *************************************************************************
code_path    = pwd;
results_path = [code_path,   '\optimization results\'];
logs_path    = [code_path,   '\optimization logs\'];
add_path     = [code_path,   '\optimization routines\'];
excel_file   = [results_path,'Optimization results diags.xlsx'];
outfile      = [logs_path,   'Optimization results diags.txt'];
fid          = fopen(outfile,'w'); fclose(fid);

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
% own=own_dummies; clear owndummies

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

% *************************************************************************
% Tables of diagnostics
% *************************************************************************
optrout_no_all = [];
perturbs_all   = [];
counts_all     = [];
exit_infos_all = [];
tocs_all       = [];

fvals_all      = [];
thetas_all     = [];
grad_all       = [];
hess1_eigs     = [];
hess2_eigs     = [];
std_errors_all = [];
conv_crit_all  = [];

% *************************************************************************
% Loop over optimization routines
% *************************************************************************
optrouts =[1;3;4;5;14;15];

for kmkm = 1:1:size(optrouts,1)
    
    optrout=optrouts(kmkm,1);
    
    if optrout<=9
        matfile = ['blp_0',num2str(optrout),'_data_optim'];
    else
        matfile = ['blp_',num2str(optrout),'_data_optim'];
    end
    
    cd(results_path)
    load (matfile, 'theta1s', 'theta2s','deltas','fvals','perturbs2',...
        'gmmresids','exit_infos','counts2','tocs');
    cd(code_path)
    
    % Track information collected in main.m
    optrout_no      = repmat(optrout,size(fvals,1),1);
    optrout_no_all  = [optrout_no_all;   optrout_no];
    perturbs_all    = [perturbs_all;      perturbs2];
    counts_all      = [counts_all;          counts2];
    exit_infos_all  = [exit_infos_all;   exit_infos];
    tocs_all        = [tocs_all;               tocs];
    
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
        
        theta2w        = zeros(mydims,mydims);
        theta2w(:,1)   = theta2_est;
        expmu          = exp(mufunc(x2,theta2w));
        s_jt_repl      = mktsh(exp(delta_repl),expmu);
        
        fprintf('===============================\n')
        fprintf('           Diagnostics         \n');
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
        
        % Hessian eigenvalues
        hess1_eigs_tmp = NaN*ones(1,mydims)';
        hess2_eigs_tmp = NaN*ones(1,mydims)';
        
        % Convergence criterion of the form g'inv(H)g
        conv_crit1 = NaN;
        conv_crit2 = NaN;
        
        % GMM objective function value when calculating hessian
        fval_hess  = NaN;
        
        % *************************************************************
        % Numerical Hessian calculations start here
        % *************************************************************
        flag_hess = 0;
        if (grad_norm_inf<=0.1 || min_index==1)
            flag_hess = 1;
        end
        
        if flag_hess == 1
            
            fprintf('Hessian calculations\n');
            
            mvalold0   = mvalolds(:,1);
            mvalold00  = mvalolds(:,2);
            
            hess1 = eye(mydims);
            
            mytolfun   = 1e-3;
            mytolx     = 1e-3;
            
            if ((optrout==12) || (optrout==14) || (optrout==15))
                disp('Tight tolerance');
                mytolfun = 1e-6;
                mytolx   = 1e-6;
            end
            
            % Use a single iteration of Quasi-Newton to calculate
            % the numerical Hessian
            options = optimset(...
                'MaxIter'    , 0,...
                'GradObj'    , 'on',...
                'HessUpdate' , 'bfgs',...
                'LargeScale' , 'off',...
                'TolFun'     ,  mytolfun,...
                'TolX'       ,  mytolx,...
                'Display'    ,'iter');
            
            mvalold = mvalold00;
            if ((optrout==13) || (optrout==14) || (optrout==15))
                [theta2_hess,fval_hess,exit_info_hess,tmp_hess,g_hess,hess2]=fminunc('gmmobj2_dfs',theta2_est,options);
            else
                [theta2_hess,fval_hess,exit_info_hess,tmp_hess,g_hess,hess2]=fminunc('gmmobj2',theta2_est,options);
            end
            
            % Hessian eigenvalues
            hess1_eigs_tmp = eig(reshape(hess1,mydims,mydims));
            hess2_eigs_tmp = eig(reshape(hess2,mydims,mydims));
            
            % Convergence criterion of the form g'inv(H)g
            conv_crit1 = grad_repl'*inv(reshape(hess1,mydims,mydims))*grad_repl;
            conv_crit2 = grad_repl'*inv(reshape(hess2,mydims,mydims))*grad_repl;
            
            fprintf('Hessian eigenvalues 1\n');
            printm(hess1_eigs_tmp);
            fprintf('Hessian eigenvalues 2\n');
            printm(hess2_eigs_tmp);
            
        end
        
        fvals_all      = [fvals_all      ; fval_repl,fval_est,fval_hess];
        thetas_all     = [thetas_all     ; [theta1_repl',theta2_est']  ];
        grad_all       = [grad_all       ; grad_repl'                  ];
        std_errors_all = [std_errors_all ; std_errors                  ];
        hess1_eigs     = [hess1_eigs     ; hess1_eigs_tmp'             ];
        hess2_eigs     = [hess2_eigs     ; hess2_eigs_tmp'             ];
        conv_crit_all  = [conv_crit_all  ; [grad_norm_inf,conv_crit1,conv_crit2,max(abs(s_jt_repl-s_jt))]];
        
    end
end


toc
diary off

%**************************************************************************
% Create cells for the various matrices to be saved in an Excel file
%**************************************************************************
optrout_no_all_cell     = num2cell(optrout_no_all);
perturbs_all_cell       = num2cell(perturbs_all);
counts_all_cell         = num2cell(counts_all);
exit_infos_all_cell     = num2cell(exit_infos_all);
tocs_all_cell           = num2cell(tocs_all);

fvals_all_cell          = num2cell(fvals_all);
thetas_all_cell         = num2cell(thetas_all);
grad_all_cell           = num2cell(grad_all);
std_errors_all_cell     = num2cell(std_errors_all);

hess1_eigs_cell         = num2cell(hess1_eigs);
hess2_eigs_cell         = num2cell(hess2_eigs);
conv_crit_all_cell      = num2cell(conv_crit_all);

fvals_header           =  {'optmethod','stvalue','fcn_evals','exit_info','toc','fval','fval_est','fval_hess','grad_norm_inf','convcrit1','convcrit2','share_norm_inf'};
thetas_header           = {'optmethod','stvalue','price_mean','const_mean','hpwt_mean','air_mean','mpg_mean','space_mean','price_sigma','const_sigma','hpwt_sigma','air_sigma','mpg_sigma'};
stderrs_header          = {'optmethod','stvalue','price_mean_se','const_mean_se','hpwt_mean_se','air_mean_se','mpg_mean_se','space_mean_se','price_sigma_se','const_sigma_se','hpwt_sigma_se','air_sigma_se','mpg_sigma_se'};
gradients_header        = {'optmethod','stvalue','price_sigma','const_sigma','hpwt_sigma','air_sigma','mpg_sigma'};
hessians_header         = {'optmethod','stvalue','eig1','eig2','eig3','eig4','eig5'};

fvals_data              = [optrout_no_all_cell,perturbs_all_cell,counts_all_cell,exit_infos_all_cell,tocs_all_cell,fvals_all_cell,conv_crit_all_cell];
fvals_excel             = [fvals_header;fvals_data];

thetas_data             = [optrout_no_all_cell,perturbs_all_cell,thetas_all_cell];
thetas_excel            = [thetas_header;thetas_data];

stderrs_data            = [optrout_no_all_cell,perturbs_all_cell,std_errors_all_cell];
stderrs_excel           = [stderrs_header;stderrs_data];

gradients_data          = [optrout_no_all_cell,perturbs_all_cell,grad_all_cell];
gradients_excel         = [gradients_header;gradients_data];

hessians_data           = [optrout_no_all_cell,perturbs_all_cell,hess1_eigs_cell];
hessians_excel          = [hessians_header;hessians_data];

hessians2_data          = [optrout_no_all_cell,perturbs_all_cell,hess2_eigs_cell];
hessians2_excel         = [hessians_header;hessians2_data];

xlswrite(excel_file,fvals_excel,          'fvals');
xlswrite(excel_file,thetas_excel,        'thetas');
xlswrite(excel_file,stderrs_excel,      'stderrs');
xlswrite(excel_file,gradients_excel,  'gradients');
xlswrite(excel_file,hessians_excel,    'hessians');
xlswrite(excel_file,hessians2_excel,  'hessians2');