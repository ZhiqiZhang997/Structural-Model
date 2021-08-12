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
mydims   = 13;

%**************************************************************************
% Define globals
%**************************************************************************
global invA ns x1 x2 v demogr s_jt IV vfull dfull...
    theta1 theti thetj cdid cdindex optrout...
    mymaxfunevals mvalold...
    oldt2 gmmresid mvalolds mvalold0 mvalold00 mvalolds2...
    fcnevals grad_est delta_est delta_est0 jacob_est mvalold_logit...
    mval_track y

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
% Load data
% *************************************************************************
load ps2
load iv
IV = [iv(:,2:21) x1(:,2:25)];
clear iv

ns   = 20;
nmkt = 94;
nbrn = 24;

cdid    = kron([1:nmkt]',ones(nbrn,1));
cdindex = [nbrn:nbrn:nbrn*nmkt]';

% *********************************************************************
% Logit IV regression
% *********************************************************************
dfull = demogr(cdid,:);
vfull = v(cdid,:);

invA = inv([IV'*IV]);

temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
outshr = 1.0 - sum1(cdid,:);

y   = log(s_jt) - log(outshr);
mid = x1'*IV*invA*IV';
t   = inv(mid*x1)*mid*y;
mvalold_logit = x1*t;

n   = size(x1,1);
k   = size(x1,2);
ESS = y'*y-2*t'*x1'*y+t'*x1'*x1*t;
s2  = ESS/(n-k);
A   = (x1'*(IV*inv(IV'*IV)*IV')*x1);
se  = sqrt(diag(s2*inv(A)));

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
for optrout=1:11
    
    if (optrout~=6)
        
        if optrout<=9
            matfile=['nevo_0',num2str(optrout),'_data_optim'];
        else
            matfile=['nevo_',num2str(optrout),'_data_optim'];
        end
        
        cd(results_path);
        load (matfile, 'theta1s', 'theta2s','deltas','fvals','perturbs2',...
            'gmmresids','exit_infos','counts2','tocs');
        cd(code_path);
                
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
            
            % Nevo (JEMS,2000) estimates
            theta2w = [0.377,3.089, 0,1.186,0;
                1.848, 16.598, -0.659,0,11.625;
                0.004,-0.193,0,0.029,0;
                0.081, 1.468,0,-1.514,0];
            
            % Initialize theta2s
            [theti, thetj, theta2] = find(theta2w);
            
            % Track information collected in main.m
            perturb      = perturbs2(jj_optrout  ,:) ;
            theta1_est   = theta1s(jj_optrout    ,:)';
            theta2_est   = theta2s(jj_optrout    ,:)';
            delta_est    = deltas(:,jj_optrout   ,:) ;
            fval_est     = fvals(jj_optrout      ,:) ;
            min_index    = fval_est==min(min(fvals)) ;            
            gmmresid_est = gmmresids(:,jj_optrout,:) ;
            
            temp1        = gmmresid_est'*IV;
            fval_est2    = temp1*invA*temp1';
            
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
            
            temp1          = gmmresid_repl'*IV;
            fval_repl2     = temp1*invA*temp1';
            delta_repl     = log(mvalolds(:,2));
            delta_repl2    = gmmresid_repl+x1*theta1;
            
            temp1          = x1'*IV;
            temp2          = delta_repl'*IV;
            theta1_repl2   = inv(temp1*invA*temp1')*temp1*invA*temp2';
            gmmresid_repl2 = delta_repl - x1*theta1;
            
            theta2w        = full(sparse(theti,thetj,theta2_est));
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
%                 mvalold = mvalold0;
%                 if (optrout>=13 && optrout<=15)
%                     hess1 = nevo_hessian_mat_new('gmmobj2_dfs',theta2_est);
%                 else
%                     hess1 = nevo_hessian_mat_new('gmmobj2',theta2_est);
%                 end
                
                mytolfun   = 1e-3;
                mytolx     = 1e-3;                
                
                if ((optrout==14) || (optrout==15) || (optrout==16))
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
                if ((optrout==13) || (optrout==14) || (optrout==16))
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

fvals_header           = {'optmethod','stvalue','fcn_evals','exit_info','toc','fval','fval_est','fval_hess','grad_norm_inf','convcrit1','convcrit2','share_norm_inf'};
thetas_header          = {'optmethod','stvalue','price','brand1','brand2','brand3','brand4','brand5','brand6',...
    'brand7','brand8','brand9','brand10','brand11','brand12',...
    'brand13','brand14','brand15','brand16','brand17','brand18',...
    'brand19','brand20','brand21','brand22','brand23','brand24',...
    'const_sigma','price_sigma','sugar_sigma','mushy_sigma',...
    'const_inc','price_inc','sugar_inc','mushy_inc','price_inc2',...
    'const_age','sugar_age','mushy_age','price_child'};

stderrs_header          = {'optmethod','stvalue','price_se','brand1_se','brand2_se','brand3_se','brand4_se','brand5_se','brand6_se',...
    'brand7_se','brand8_se','brand9_se','brand10_se','brand11_se','brand12_se',...
    'brand13_se','brand14_se','brand15_se','brand16_se','brand17_se','brand18_se',...
    'brand19_se','brand20_se','brand21_se','brand22_se','brand23_se','brand24_se',...
    'const_sigma_se','price_sigma_se','sugar_sigma_se','mushy_sigma_se',...
    'const_inc_se','price_inc_se','sugar_inc_se','mushy_inc_se','price_inc2_se',...
    'const_age_se','sugar_age_se','mushy_age_se','price_child_se'};

gradients_header        = {'optmethod','stvalue','const_sigma','price_sigma','sugar_sigma','mushy_sigma',...
    'const_inc','price_inc','sugar_inc','mushy_inc','price_inc2',...
    'const_age','sugar_age','mushy_age','price_child'};

hessians_header         = {'optmethod','stvalue','eig1','eig2','eig3','eig4','eig5','eig6','eig7','eig8','eig9','eig10','eig11','eig12','eig13'};

fvals_data              = [optrout_no_all_cell,perturbs_all_cell,counts_all_cell,exit_infos_all_cell,tocs_all_cell,fvals_all_cell,conv_crit_all_cell];
fvals_excel             = [fvals_header;fvals_data];

thetas_data             = [optrout_no_all_cell,perturbs_all_cell,thetas_all_cell];
thetas_excel            = [thetas_header;thetas_data];

std_errors_data         = [optrout_no_all_cell,perturbs_all_cell,std_errors_all_cell];
std_errors_excel        = [stderrs_header;std_errors_data];

gradients_data          = [optrout_no_all_cell,perturbs_all_cell,grad_all_cell];
gradients_excel         = [gradients_header;gradients_data];

hessians_data           = [optrout_no_all_cell,perturbs_all_cell,hess1_eigs_cell];
hessians_excel          = [hessians_header;hessians_data];

hessians2_data          = [optrout_no_all_cell,perturbs_all_cell,hess2_eigs_cell];
hessians2_excel         = [hessians_header;hessians2_data];

xlswrite(excel_file,fvals_excel,          'fvals');
xlswrite(excel_file,thetas_excel,        'thetas');
xlswrite(excel_file,std_errors_excel, 'stderrors');
xlswrite(excel_file,gradients_excel,  'gradients');
xlswrite(excel_file,hessians_excel,    'hessians');
xlswrite(excel_file,hessians2_excel,  'hessians2');