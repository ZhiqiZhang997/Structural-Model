%**************************************************************************
% To accompany Knittel and Metaxoglou (2008)
% Estimation of Random Coefficient Demand Models:
% Challenges, Difficulties and Warnings
% Knittel      : crknittel@ucdavis.edu
% Metaxoglou   : konstantinos.metaxoglou@bateswhite.com
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
for optrout=5:5
    
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
        for jj_optrout=42
            
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
                % mvalold = mvalold0;
                % if (optrout>=13 && optrout<=15)
                % hess1 = nevo_hessian_mat_new('gmmobj2_dfs',theta2_est);
                % else
                % hess1 = nevo_hessian_mat_new('gmmobj2',theta2_est);
                % end
                
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
                
                
                [V, D] = eig(hess2);
                D_orig = D;
                V_orig = V;
                % find the maximum eigenvalue
                maxEV = max(abs(D(:)));
                % normalize all eigenvalues
                D = D/maxEV;
                % get the column1 by taking the diagonal entries and putting them in a vector
                column1 = diag(D);
                column2 = NaN * ones(size(V));
                % iterate through each eigenvector
                for hh = 1:size(V, 2)
                    % find the index of the largest abs. value in eigenvector i, corresponding to eigenvalue i
                    [maxx, indx] = max(abs(V(:, hh)));
                    % Norm-2 of eigenvectors is equal to 1
                    printm(norm(V(:,hh),2));
                    % Length of hess2 x eigenvector i is equal to eigenvalue i
                    printm(norm(hess2*V_orig(:,hh))-D_orig(hh,hh));
                    % It should also be the case that the inner product of
                    % the eigenvectors is equal to zero b/c they are
                    % perpendicular
                    % get the corresponding value of that index in column 2
                    column2(hh, 1) = V(indx, hh);
                    column2(hh, 2) = indx;
                end
            end
        end
    end
end
toc

theta2s_header={'const_sigma','price_sigma','sugar_sigma','mushy_sigma',...
    'const_inc','price_inc','sugar_inc','mushy_inc','price_inc2',...
    'const_age','sugar_age','mushy_age','price_child'};

eigvalue_normal = column1;
eigvec_normal   = column2(:,1);
eigvec_elem     = column2(:,2);
temp            = theta2s_header';
temp            = temp(eigvec_elem,:);

head_xls = {'var','eigvalue_norm','extrem_elem','extrem_val',};
data_xls = [cellstr(temp),num2cell([eigvalue_normal,eigvec_elem,eigvec_normal])];
data_xls = [head_xls;data_xls];

cd(results_path)

xlswrite('hessian_diags.xlsx',data_xls,'hessian_eigsystem');

head_xls = cellstr(theta2s_header);
data_xls = num2cell(hess2);
data_xls = [head_xls;data_xls];
xlswrite('hessian_diags.xlsx',data_xls,'hessian');

head_xls = {'var','hess_diagonal'};
data_xls = [cellstr(theta2s_header'),num2cell(diag(hess2))];
data_xls = [head_xls;data_xls];
xlswrite('hessian_diags.xlsx',data_xls,'hessian_diagonal');

head_xls = {'hessian_cond_num'};
data_xls = num2cell(cond(hess2));
data_xls = [head_xls;data_xls];
xlswrite('hessian_diags.xlsx',data_xls,'hessian_cond_num');

data_xls = num2cell(V);
xlswrite('hessian_diags.xlsx',data_xls,'hessian_eigvectors');





