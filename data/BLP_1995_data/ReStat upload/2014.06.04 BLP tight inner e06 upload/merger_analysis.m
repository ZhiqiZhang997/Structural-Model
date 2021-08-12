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
% Set theta2_est=0 in line 184 for simple logit results
% *************************************************************************
flag_check      = 0; % if equal to 1, post and 1pct calculations are identical
flag_best       = 0; % if equal to 1, limit analysis to params for min fval
flag_full_simul = 0; % if equal to 1, full merger simulation analysis
p_inc           = 1.01;
code_path       = pwd;
results_path    = [code_path,'\Optimization results\' ];
merger_path     = [code_path,'\Merger results\'   ];
logs_path       = merger_path;
add_path        = [code_path,'\Optimization routines\'];
outfile         = [logs_path,'Merger analysis.txt'    ];
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
optrouts = [1,3,4,5,14,15]';

for km_optrout=1:1:size(optrouts,1)
    
    optrout = optrouts(km_optrout,1);
    
    if optrout~=6
        
        cd(results_path)
        
        if optrout<=9
            matfile           = ['blp_0',num2str(optrout),'_data_optim'];
            merger_file       = [merger_path,'blp_merger_results_0',num2str(optrout),'.txt'];
            simul_file_xls    = [merger_path,'blp_merger_simul_xls_0',num2str(optrout),'.xlsx'];
        else
            matfile           = ['blp_',num2str(optrout),'_data_optim'];
            merger_file       = [merger_path,'blp_merger_results_',num2str(optrout),'.txt'];
            simul_file_xls    = [merger_path,'blp_merger_simul_xls',num2str(optrout),'.xlsx'];
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
            
            %Consumer surplus pre-merger
            exp_V  = ind_eg(exp(deltajt_pre),exp(mufunc(x2,theta2w)));
            tmp    = [];
            CV_pre = [];
            
            for i=1:nmkt
                alphai_tmp = -alphai(cdid==i,:);
                alphai_tmp = alphai_tmp(1,:);
                tmp(i,:)   = log(sum(exp_V(cdid==i,:))+1)./alphai_tmp;
                CV_pre(i,:)= tmp(i,:);
            end
            
            %Market power pre_merger
            own_dummy_pre = own;
            price_pre     = price;
            
            mm=[];
            for i=1:max(cdid)
                p  = price_pre(cdid==i,:);
                s  = sjt_pre(cdid==i,:);
                nn = nbrn(i);
                om = deriv_all(1:nn,1:nn,i).*(own_dummy_pre(cdid==i,:)*own_dummy_pre(cdid==i,:)');
                m  = -inv(om')*s;
                mm = [mm;m];
            end
            
            margin_pre    = mm;
            mc            = price_pre-margin_pre;
            profit_pre    = margin_pre.*sjt_pre;
            
            optrout_aux   = repmat(optrout,size(price,1),1);
            startval_aux  = repmat(perturb,size(price,1),1);
            fval_aux      = repmat(fval_repl,size(price,1),1);
            market        = cdid;
            
            % Approximate post-merger prices
            tmp            = own(:,16)+own(:,19);
            own_dummy_post = [own(:,1:15),tmp,own(:,17:18),own(:,20:26)];
            
            mm=[];
            for i=1:max(cdid)
                p  = price_pre(cdid==i,:);
                s  = sjt_pre(cdid==i,:);
                nn = nbrn(i);
                om = deriv_all(1:nn,1:nn,i).*(own_dummy_post(cdid==i,:)*own_dummy_post(cdid==i,:)');
                m  = -inv(om')*s;
                mm = [mm;m];
            end
            
            price_approx = mc+mm;
            price_post   = price_approx;
            margin_post  = mm;
            
            % Price used to simulate 1% price increase
            price_agg  = price_pre*p_inc;
            margin_agg = price_agg-mc;
            
            % All the economic variables of interest should be identical
            % for pre, approx and post if flag_check==1
            if flag_check == 1
                price_approx   = price_pre;
                price_post     = price_pre;
                price_agg      = price_pre;
                own_dummy_post = own_dummy_pre;
            end
            
            % *************************************************************
            % Full merger simulation
            % *************************************************************
            if flag_full_simul==1
                diags = [];
                for mkt_no=1:max(cdid)
                    tic
                    fprintf('=================================\n');
                    fprintf('Merger simulation market: %5i\n',mkt_no);
                    fprintf('=================================\n');
                    
                    options = optimset('Display','iter',...
                        'MaxFunEvals',200000,...
                        'MaxIter'    ,1000,...
                        'TolFun'     ,1e-16,...
                        'TolX'       ,1e-16);
                    
                    deltajt0       = deltajt_pre;
                    p_to_solve     = price_approx(cdid==mkt_no,:);
                    own_dummy      = own_dummy_post;
                    xxx            = merger_focs(p_to_solve);
                    
                    [price_post_tmp,fval,exitflag,output,jacobians] = fsolve('merger_focs_new',p_to_solve,options);
                    
                    fprintf('\n norm_inf_focs: %12.4f\t  norm_inf_price: %12.6f\n',max(abs(xxx)),max(abs(price_post_tmp-price_pre(cdid==mkt_no))));
                    
                    toc_tmp    = toc;
                    price_post(cdid==mkt_no)=price_post_tmp;
                    diags_tmp  = [max(abs(fval)),toc_tmp,rank(jacobians),exitflag,output.iterations,output.funcCount,output.firstorderopt];
                    diags      = [diags;diags_tmp];
                    
                end
                margin_post     = price_post-mc;
                head_diags      = cellstr({'max_abs_fval','toc','jacob_rank','exit_flag','iters','count','conv_crit'});
                cd(merger_path);
                xlswrite(simul_file_xls,[head_diags;num2cell(diags)]);
                cd(code_path);
            end
            
            % Mean utility levels post-merger
            % Mean utility levels 1% price increase
            deltajt_post = deltajt_pre-price_pre*theta1(1)+price_post*theta1(1);
            deltajt_agg  = deltajt_pre-price_pre*theta1(1)+price_agg*theta1(1);
            
            x2_post = x2;
            x2_agg  = x2;
            
            x2_post(:,1) = price_post;
            x2_agg(:,1)  = price_agg;
            
            %calculate implied market shares
            theta2w      = zeros(5,5);
            theta2w(:,1) = theta2;
            
            % Deviations from mean utility post-merger
            % Deviations from mean utility 1% price increase
            [n k]  = size(x2_post);
            j      = size(theta2w,2)-1;
            mu     = zeros(n,ns);
            mu_agg = zeros(n,ns);
            
            for i = 1:ns
                v_i        = vfull(:,i:ns:k*ns);
                d_i        = dfull(:,i:ns:j*ns);
                mu(:,i)    = (x2_post.*v_i*theta2w(:,1));
                mu_agg(:,i)= ( x2_agg.*v_i*theta2w(:,1));
            end
            
            expmu       = exp(mu);
            expmu_agg   = exp(mu_agg);
            
            expmval     = exp(deltajt_post);
            expmval_agg = exp(deltajt_agg);
            
            sijt_post   = ind_sh(expmval,expmu);
            sjt_post    = (1/ns)*sum(sijt_post')';
            
            sijt_agg    = ind_sh(expmval_agg,expmu_agg);
            sjt_agg     = (1/ns)*sum(sijt_agg')';
            
            exp_V       = ind_eg(expmval,expmu);
            exp_V_agg   = ind_eg(expmval_agg,expmu_agg);
            
            tmp     = [];
            CV_post = [];
            CV_agg  = [];
            
            % Consumer surplus post-merger
            % Consumer surplus 1% price increase
            for i=1:nmkt
                alphai_tmp   = -alphai(cdid==i,:);
                alphai_tmp   = alphai_tmp(1,:);
                tmp(i,:)     = log(sum(exp_V(cdid==i,:))+1)./alphai_tmp;
                CV_post(i,:) = tmp(i,:);
                tmp(i,:)     = log(sum(exp_V_agg(cdid==i,:))+1)./alphai_tmp;
                CV_agg(i,:)  = tmp(i,:);
            end
            
            mean_CV     = mean((CV_post-CV_pre)')';
            mean_CV_agg = mean((CV_agg-CV_pre)')';
            
            mean_CV_aux     = [];
            mean_CV_agg_aux = [];
            
            % Comply with dimensions of remaining vectors
            for i=1:size(nbrn,1)
                tmp             = nbrn(i);
                mean_CV_aux     = [mean_CV_aux     ; repmat(mean_CV(i,1)    ,tmp,1)];
                mean_CV_agg_aux = [mean_CV_agg_aux ; repmat(mean_CV_agg(i,1),tmp,1)];
            end
            
            % Profit post-merger
            % Profit 1% price increase
            profit_post = margin_post.*sjt_post;
            profit_agg  = margin_agg.*sjt_agg;
            
            % Aggregate elasticity calculation
            agg_el = [];
            for my_mkt=1:1:max(cdid)
                tmp_s      = [sum(sjt_agg(cdid==my_mkt)),sum(sjt_pre(cdid==my_mkt))];
                numer      = 100*(tmp_s(:,1)-tmp_s(:,2))/(0.5*(tmp_s(:,1)+tmp_s(:,2)));
                denom      = 100*(p_inc-1)/(0.5*(p_inc+1));
                agg_el_tmp = numer/denom;
                agg_el     = [agg_el; repmat(agg_el_tmp,size(sjt_agg(cdid==my_mkt),1),1)];
                
            end
            
            merger_results = [merger_results;...
                [optrout_aux,startval_aux,...
                market,model_id,firmid,id,product,share...
                price_pre,price_post,price_agg,...
                sjt_pre,  sjt_post, sjt_agg,...
                mc, elast_own,...
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
                printm(std([    sjt_pre,     sjt_post,     sjt_agg]));
                printm(std([  price_pre,   price_post,   price_agg]));
                printm(std([deltajt_pre, deltajt_post, deltajt_agg]));
            end
            
        end
        
        cd(merger_path);
        save(merger_file,'merger_results','-ASCII');
        cd(code_path);
        
    end
end

toc

