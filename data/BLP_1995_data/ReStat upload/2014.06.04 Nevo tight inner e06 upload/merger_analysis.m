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
mydims   = 13;

%**************************************************************************
% Define globals
%**************************************************************************
global invA ns x1 x2 s_jt IV vfull dfull...
    theta1 theti thetj cdid cdindex optrout...
    omega mvalold oldt2 gmmresid v demogr...
    beta mc deltajt0 market_no nmkt...
    nbrn own_dummy alphaib rcoeffs2 mvalolds price theta1_repl theta2_est...
    alphai_merger

%**************************************************************************
% Define paths for input and output
% Controls for full simulation and calculation checks
%**************************************************************************
flag_check      = 0; % if equal to 1, post and 1pct calculations are identical
flag_best       = 0; % if equal to 1, limit analysis to params for min fval
flag_full_simul = 0; % if equal to 1, full merger simulation analysis
p_inc           = 1.01;
code_path       = pwd;
results_path    = [code_path,'\Optimization results\' ];
merger_path     = [code_path,'\Merger results\'];
logs_path       = merger_path;
add_path        = [code_path,'\Optimization routines\'];
outfile         = [logs_path,'Merger analysis.txt'    ];
fid             = fopen(outfile,'w'); fclose(fid);
tic
diary(outfile)
addpath(add_path);

%**************************************************************************
% Load data
%**************************************************************************
load ps2
load iv
IV = [iv(:,2:21) x1(:,2:25)];
clear iv

ns   = 20;
nmkt = 94;
nbrn = 24;

cdid    = kron([1:nmkt]',ones(nbrn,1));
cdindex = [nbrn:nbrn:nbrn*nmkt]';

%**************************************************************************
% Logit IV regression
%**************************************************************************
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

merger = importdata('own_matrix.xls');
pre_merger = merger.data.pre_merger;
merger1 = merger.data.merger1;
merger2 = merger.data.merger2;
merger3 = merger.data.merger3;

price      = full(x1(:,1));

optrout_no_all  = [];
perturbs_all    = [];
counts_all      = [];
exit_infos_all  = [];
tocs_all        = [];

%**************************************************************************
% Loop over optimization routines
%**************************************************************************
for optrout=16:16
    
    if (optrout~=6)
        
        cd(results_path);
        
        if optrout<=9
            matfile         = ['nevo_0',num2str(optrout),'_data_optim.mat'];
            merger_file     = [merger_path,'nevo_merger_results_0',num2str(optrout),'.txt'];
            simul_file_xls = [merger_path,'nevo_merger_results_0',num2str(optrout),'.xlsx'];
        else
            matfile         = ['nevo_',num2str(optrout),'_data_optim.mat'];
            merger_file     = [merger_path,'nevo_merger_results_',num2str(optrout),'.txt'];
            simul_file_xls = [merger_path,'nevo_merger_results_',num2str(optrout),'.xlsx'];
        end
        
        load (matfile, 'theta1s', 'theta2s','deltas','fvals','perturbs2',...
            'gmmresids','exit_infos','counts2','tocs','coeffs2');
        
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
            coeffs2   = coeffs2(min_fval_ind,:);
            deltas    = deltas(:,min_fval_ind);
            gmmresids = gmmresids(:,min_fval_ind);
            fvals     = fvals(min_fval_ind,:);
            perturbs2 = perturbs2(min_fval_ind,:);
        end
        
        merger_results = [];
        coeffs20       = coeffs2;
        
        %******************************************************************
        % Loop over starting values
        %******************************************************************
        for jj_optrout=1:size(fvals,1);
            
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
            % gmmobj2 and gmmobj2_dfs use 1e-14!
            if ((optrout==13) || (optrout==14) || (optrout==16))
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
            fprintf('     Merger Diagnostics        \n');
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
            
            %**************************************************************
            % Use theta2_est and theta1_repl
            %**************************************************************
            theta2 = theta2_est;
            theta1 = theta1_repl;
            
            % MD estimation
            theta2w = full(sparse(theti,thetj,theta2));
            t       = size(se,1) - size(theta2,1);
            se2w    = full(sparse(theti,thetj,se(t+1:size(se,1))));
            
            omega   = inv(vcov(2:25,2:25));
            xmd     = [x2(1:24,1) x2(1:24,3:4)];
            ymd     = theta1(2:25);
            
            beta    = inv(xmd'*omega*xmd)*xmd'*omega*ymd;
            resmd   = ymd - xmd*beta;
            semd    = sqrt(diag(inv(xmd'*omega*xmd)));
            mcoef   = [beta(1); theta1(1); beta(2:3)];
            semcoef = [semd(1); se(1); semd];
            
            coeffs2_repl = [mcoef theta2w];
            coeffs2_est  = reshape(coeffs20(jj_optrout,:),4,6);
            
            fprintf('coeffs2 check    : %12.4f\n',max(abs(coeffs2_repl(:)-coeffs2_est(:))));
            
            % The coefficients are as follows: constant price sugar mushy
            beta  = coeffs2_repl(:,1);
            sigma = coeffs2_repl(:,2);
            pai   = coeffs2_repl(:,3:6);
            
            % Construct the various random coefficients
            rcoeffs  = NaN*zeros(size(cdid,1),ns);
            rcoeffs2 = NaN*zeros(size(cdid,1),ns);
            
            for i=1:size(rcoeffs,1)
                k=1;
                for j=1:ns
                    tmpd = dfull(i,j:ns:80);
                    tmpv = vfull(i,j:ns:80);
                    rcoeff_tmp       = pai*tmpd'+diag(sigma)*tmpv';
                    rcoeffs(i,k:k+3) = (beta+rcoeff_tmp)';
                    rcoeffs2(i,k:k+3)= (rcoeff_tmp)';
                    k=k+4;
                end
            end
            
            % Various random coefficients
            rcoeffs_const = rcoeffs(:,1:4:80);
            rcoeffs_price = rcoeffs(:,2:4:80);
            rcoeffs_sugar = rcoeffs(:,3:4:80);
            rcoeffs_mushy = rcoeffs(:,4:4:80);
            
            alphai        = rcoeffs_price;
            alphai_merger = alphai;
            
            % Replicate market shares
            theta2w    = full(sparse(theti,thetj,theta2));
            expmu      = exp(mufunc(x2,theta2w));
            s_jt_repl  = mktsh(exp(delta_repl),expmu);
            s_ijt_repl = ind_sh(exp(delta_repl),expmu);
            
            % Market shares and mean utility levels pre-merger
            s_ijt_pre   = s_ijt_repl;
            s_jt_pre    =  s_jt_repl;
            deltajt_pre = delta_repl;
            
            % Derive matrices of price derivatives and elasticities
            deriv_all = zeros(max(nbrn),max(nbrn),nmkt);
            elast_all = zeros(max(nbrn),max(nbrn),nmkt);
            
            % Loop over markets
            for i=1:max(cdid)
                l=1;
                ind     = cdid==i;
                
                % Use replicated market shares
                pjt     =        x1(ind==1,1);
                sjt     =  s_jt_pre(ind==1,1);
                alpha_i =    alphai(ind==1,:);
                s_ijt   = s_ijt_pre(ind==1,:);
                
                elast   = zeros(size(pjt,1),size(pjt,1));
                deriv   = zeros(size(pjt,1),size(pjt,1));
                
                % Loop over products in the same market
                for j=1:size(pjt,1)
                    deriv(j,j)     =                  (1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
                    elast(j,j)     = (pjt(j)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(ones(1,ns)-s_ijt(j,:)));
                    deriv_own(i,j) = deriv(j,j);
                    elast_own(i,j) = elast(j,j);
                    
                    for k=1:size(pjt,1)
                        if k~=j
                            deriv(j,k)       =                  -(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                            elast(j,k)       = -(pjt(k)./sjt(j))*(1/ns)*sum(alpha_i(j,:).*s_ijt(j,:).*(s_ijt(k,:)));
                            elast_cross(i,l) = elast(j,k);
                            deriv_cross(i,l) = deriv(j,k);
                            l=l+1;
                        end
                    end
                    
                end
                elast_all(1:size(elast,1),1:size(elast,2),i) = elast;
                deriv_all(1:size(deriv,1),1:size(deriv,2),i) = deriv;
            end
            
            % Track elasticities
            temp = [];
            temp2 = [];
            for j=1:nmkt
                temp  = [temp; (elast_all(:,:,j))];
                temp2 = [temp2; diag(elast_all(:,:,j))];
            end
            
            elast_all = temp;
            elast_own = temp2;
            
            % Consumer surplus pre-merger
            [n k]   = size(x2);
            j       = size(theta2w,2)-1;
            mu      = zeros(n,ns);
            theta2w = full(sparse(theti,thetj,theta2));
            
            for i = 1:ns
                v_i     = vfull(:,i:ns:k*ns);
                d_i     = dfull(:,i:ns:j*ns);
                mu(:,i) = (x2.*v_i*theta2w(:,1))+x2.*(d_i*  theta2w(:,2:j+1)')*ones(k,1);
            end
            
            mu    = mu+repmat(deltajt_pre,1,ns);
            V_exp = exp(mu);
            for i=1:nmkt
                alphai_tmp  = -alphai(cdid==i,:);
                alphai_tmp  = alphai_tmp(1,:);
                tmp(i,:)    = log(sum(V_exp(cdid==i,:))+1)./alphai_tmp;
                CV_pre(i,:) = tmp(i,:);
            end
            
            % Market power pre-merger
            own_dummy_pre = pre_merger;
            price_pre     = price;
            
            mm = [];
            for i=1:max(cdid)
                p = price_pre(cdid==i,:);
                s = s_jt_pre(cdid==i,:);
                om = deriv_all(:,:,i).*(own_dummy_pre*own_dummy_pre');
                m  = -inv(om')*s;
                mm = [mm;m];
            end
            
            margin_pre = mm;
            mc         = price_pre-margin_pre;
            profit_pre = margin_pre.*s_jt_pre;
            
            % Approximate post-merger price
            own_dummy_post = merger2;
            mm = [];
            
            for i=1:max(cdid)
                p  = price_pre(cdid==i,:);
                s  = s_jt_pre(cdid==i,:);
                om = deriv_all(:,:,i).*(own_dummy_post*own_dummy_post');
                m  = -inv(om')*s;
                mm = [mm;m];
            end
            
            price_approx = mc+mm;            
            price_post   = price_approx;
            
            % Price used to simulate 1% price increase
            price_agg    = price_pre*p_inc;
            
            % All the economic variables of interest should be identical
            % for pre, approx and post if flag_check==1
            if flag_check == 1
                price_approx   = price_pre;
                price_post     = price_pre;
                price_agg      = price_pre;
                own_dummy_post = own_dummy_pre;
            end
            
            %**************************************************************
            % Full merger simulation
            %**************************************************************
            if flag_full_simul==1
                diags = [];
                for market_no=1:max(cdid)
                    tic
                    fprintf('=================================\n');
                    fprintf('Merger simulation market: %5i\n',market_no);
                    fprintf('=================================\n');
                    
                    options = optimset('Display','iter',...
                        'TolFun',1e-16,...
                        'TolX',  1e-16);
                                        
                    deltajt0       = deltajt_pre;
                    p_to_solve     = price_pre(cdid==market_no,:);
                    own_dummy      = own_dummy_post;
                    xxx            = merger_focs(p_to_solve);                   
                    [price_post_tmp,fval,exitflag,output,jacobians] = fsolve('merger_focs',p_to_solve,options);                    
                    
                    fprintf('\n norm_inf_focs: %12.4f\t  norm_inf_price: %12.6f\n',max(abs(xxx)),max(abs(price_post_tmp-price_pre(cdid==market_no))));
                    
                    toc_tmp    = toc;
                    price_post(cdid==market_no)=price_post_tmp;
                    diags_tmp  = [max(abs(fval)),toc_tmp,rank(jacobians),exitflag,output.iterations,output.funcCount,output.firstorderopt];
                    diags      = [diags;diags_tmp];
                end
                head_diags      = cellstr({'max_abs_fval','toc','jacob_rank','exit_flag','iters','count','conv_crit'});
                cd(merger_path);
                xlswrite(simul_file_xls,[head_diags;num2cell(diags)]);
                cd(code_path);
                fprintf('\n');
            end
            
            % Mean utility levels post-merger
            % Mean utility levels 1% price increase
            deltajt_post = deltajt_pre-price_pre*theta1(1)+price_post*theta1(1);
            deltajt_agg  = deltajt_pre-price_pre*theta1(1)+price_agg*theta1(1);
            
            x2_post = x2;
            x2_post(:,2) = price_post;
            
            x2_agg = x2;
            x2_agg(:,2)  = price_agg;
            
            theta2w    = full(sparse(theti,thetj,theta2));
            expmu      = exp(mufunc(x2_post,theta2w));
            s_jt_post  = mktsh(exp(deltajt_post),expmu);
            s_ijt_post = ind_sh(exp(deltajt_post),expmu);
            
            theta2w    = full(sparse(theti,thetj,theta2));
            expmu      = exp(mufunc(x2_agg,theta2w));
            s_jt_agg   = mktsh(exp(deltajt_agg),expmu);
            s_ijt_agg  = ind_sh(exp(deltajt_agg),expmu);
            
            % Deviations from mean utility post-merger
            % Deviations from mean utility 1% price increase
            [n k]   = size(x2);
            j       = size(theta2w,2)-1;
            mu      = zeros(n,ns);
            mu_agg  = zeros(n,ns);
            theta2w = full(sparse(theti,thetj,theta2));
            
            for i = 1:ns
                v_i         = vfull(:,i:ns:k*ns);
                d_i         = dfull(:,i:ns:j*ns);
                mu(:,i)     = (x2_post.*v_i*theta2w(:,1))+x2_post.*(d_i*theta2w(:,2:j+1)')*ones(k,1);
                mu_agg(:,i) = (x2_agg.*v_i*theta2w(:,1))+x2_agg.*(d_i*theta2w(:,2:j+1)')*ones(k,1);
            end
            
            mu        = mu+repmat(deltajt_post,1,ns);
            mu_agg    = mu_agg+repmat(deltajt_agg,1,ns);
            
            V_exp     = exp(mu);
            V_exp_agg = exp(mu_agg);
            
            % Consumer surplus post-merger
            % Consumer surplus 1% price increase
            for i=1:nmkt
                alphai_tmp   = -alphai(cdid==i,:);
                alphai_tmp   =  alphai_tmp(1,:);
                
                tmp(i,:)     = log(sum(V_exp(cdid==i,:))+1)./alphai_tmp;
                CV_post(i,:) = tmp(i,:);
                
                tmp(i,:)     = log(sum(V_exp_agg(cdid==i,:))+1)./alphai_tmp;
                CV_agg(i,:)  = tmp(i,:);
            end
            
            mean_CV      = mean((CV_post-CV_pre)')';
            mean_CV_agg  = mean((CV_agg-CV_pre)')' ;
            
            % Comply with dimensions of remaining vectors
            mean_CV_aux     = kron(mean_CV    ,ones(24,1));
            mean_CV_agg_aux = kron(mean_CV_agg,ones(24,1));
            
            % Market power post-merger
            % Market power 1% price increase
            margin_post = price_post-mc;
            margin_agg  = price_agg-mc;
            profit_post = margin_post.*s_jt_post;
            profit_agg  = margin_agg.*s_jt_agg;
            
            optrout_aux = repmat(optrout,2256,1);
            market      = kron((1:1:94)',ones(nbrn,1));
            brand       = repmat((1:1:24)',94,1);
            
            % Aggregate elasticity calculation
            agg_el = [];
            for my_mkt=1:1:max(cdid)
                tmp_s      = [sum(s_jt_agg(cdid==my_mkt)),sum(s_jt_pre(cdid==my_mkt))];
                numer      = 100*(tmp_s(:,1)-tmp_s(:,2))/(0.5*(tmp_s(:,1)+tmp_s(:,2)));
                denom      = 100*(p_inc-1)/(0.5*(p_inc+1));
                agg_el_tmp = numer/denom;
                agg_el     = [agg_el; repmat(agg_el_tmp,size(s_jt_agg(cdid==my_mkt),1),1)];
            end
            
            jj_merger_aux  = repmat(perturb,2256,1);
            merger_results = [merger_results;...
                [optrout_aux,jj_merger_aux,market,brand,s_jt,...
                price_pre, price_post,price_agg,...
                s_jt_pre,  s_jt_post, s_jt_agg,...
                mc,elast_own...
                profit_pre,profit_post,profit_agg,...
                mean_CV_aux,mean_CV_agg_aux,...
                agg_el]];
            
            fprintf('min elast_own    : %12.4f\n',min(elast_own)   );
            fprintf('min mc           : %12.4f\n',min(mc)          );
            fprintf('mean elast_own   : %12.4f\n',mean(elast_own)  );
            fprintf('mean mc          : %12.4f\n',mean(mc)         );
            fprintf('max elast_own    : %12.4f\n',max(elast_own)   );
            fprintf('max mc           : %12.4f\n',max(mc)          );
            fprintf('std elast_own    : %12.4f\n',std(elast_own)   );
            fprintf('std mc           : %12.4f\n',std(mc)          );
            fprintf('median elast_own : %12.4f\n',median(elast_own));
            fprintf('median mc        : %12.4f\n',median(mc)       );
            
            if flag_check ==1
                fprintf('   \nChecks on shares prices and deltas\n');
                fprintf('     Pre          Post      1-percent\n');
                printm(std([  s_jt_pre,     s_jt_post,    s_jt_agg]));
                printm(std([  price_pre,   price_post,   price_agg]));
                printm(std([deltajt_pre, deltajt_post, deltajt_agg]));
            end
        end
        
        cd(merger_path)
        save(merger_file,'merger_results','-ASCII');
        cd(code_path)
    end
end
fprintf('\n\n');
toc
diary off