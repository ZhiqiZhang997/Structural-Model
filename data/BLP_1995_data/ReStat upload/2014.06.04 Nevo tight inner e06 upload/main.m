%**************************************************************************
% To accompany Knittel and Metaxoglou
%**************************************************************************

clear all
close all
close hidden
format long
warning off all
clc

km = maxNumCompThreads(1);

%**************************************************************************
%Define globals
%**************************************************************************
global invA ns x1 x2 v demogr s_jt IV vfull dfull...
    theta1 theti thetj cdid cdindex optrout...
    mymaxfunevals mvalold...
    oldt2 gmmresid mvalolds mvalold0 mvalold00 mvalolds2...
    fcnevals grad_est delta_est delta_est0 jacob_est mvalold_logit

% *************************************************************************
% Define paths for codes, optimization results and logs
% *************************************************************************
code_path    =pwd;
results_path =[code_path,'\optimization results\'];
logs_path    =[code_path,'\optimization logs\'];
add_path     =[code_path,'\optimization routines\'];
addpath(add_path);

% *************************************************************************
% Loop over optimization routines
% *************************************************************************
for optrout=3:5
    
    mytolx        = 1e-6;
    mytolfun      = 1e-6;
    mymaxiters    = 5*10^5;
    mymaxfunevals = 4000;
    perturbs      = (1:1:50)';
    
    if optrout==1
        perturbs=perturbs(perturbs~=6,:);
        perturbs=perturbs(perturbs~=22,:);
        perturbs=perturbs(perturbs~=39,:);
    end
    if optrout==2
        perturbs=perturbs(perturbs~=6,:);
        perturbs=perturbs(perturbs~=28,:);
        perturbs=perturbs(perturbs~=36,:);
    end
    
    if optrout==12 || optrout==14
        perturbs(perturbs==20)=1;
        perturbs(perturbs==19)=1;
    end
    
    if optrout==15
        perturbs(perturbs==19)=1;
    end
    
    if optrout<=9
        outfile=[logs_path,['nevo_0',num2str(optrout),'_optim_log.txt']];
        matfile=['nevo_0',num2str(optrout),'_data_optim.mat'];
    else
        outfile=[logs_path,['nevo_',num2str(optrout),'_optim_log.txt']];
        matfile=['nevo_',num2str(optrout),'_data_optim.mat'];
    end
    
    fid = fopen(outfile,'w'); fclose(fid);
    
    % *********************************************************************
    % Initialize log files and matrices containing various results
    % *********************************************************************
    perturbs2  =[];                  %store perturbation set number
    fvals      =[];                  %store GMM values
    exit_infos =[];                  %store exit info
    counts2    =[];                  %store function evaluations
    tocs       =[];                  %store time of completion
    theta1s    =[];                  %store theta1s
    theta2s    =[];                  %store theta2s
    coeffs2    =[];                  %store MD estimates
    deltas     =[];                  %store deltas
    gmmresids  =[];                  %store gmm residuals
    std_errors =[];                  %store std.errors
    gradients  =[];                  %store analytical gradients
    
    % *********************************************************************
    % Load data
    % *********************************************************************
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
    t    = inv(mid*x1)*mid*y;
    mvalold_logit = x1*t;
    
    n   = size(x1,1);
    k   = size(x1,2);
    ESS = y'*y-2*t'*x1'*y+t'*x1'*x1*t;
    s2  = ESS/(n-k);
    A   = (x1'*(IV*inv(IV'*IV)*IV')*x1);
    se  = sqrt(diag(s2*inv(A)));
    
    % *********************************************************************
    % Start optimization routine with perturb_no different starting values
    % for theta2: normrnd(0,1,size(theta2));
    % for delta:  delta_logit+normrnd(0,stddev(delta_logit),2256,1)
    % *********************************************************************
    diary(outfile)
    
    for ppp=1:size(perturbs,1)
        tic
        perturb=perturbs(ppp,1);
        fprintf('\n');
        fprintf('========================================================================\n')
        fprintf('                       optimization routine :%2i\n',optrout);
        fprintf('                     set of starting values :%2i\n',perturb);
        fprintf('========================================================================\n')
        
        
        % Nevo (JEMS,2000) estimates
        theta2w = [0.377,3.089, 0,1.186,0;
            1.848, 16.598, -0.659,0,11.625;
            0.004,-0.193,0,0.029,0;
            0.081, 1.468,0,-1.514,0];
        
        % initialize theta2s
        [theti, thetj, theta2]=find(theta2w);
        randn('state',1000*perturb);
        theta2   = normrnd(0,1,size(theta2));
        theta2w0 = full(sparse(theti,thetj,theta2));
        
        % initialize deltas
        randn('state',1000*perturb);
        mvalold = exp(mvalold_logit+normrnd(0,sqrt(s2),size(x1,1),1));
        oldt2   = zeros(size(theta2));
                
        fprintf('\ntheta2 starting values:\n')
        printm(theta2');
        fprintf('\n\n')
        
        % initialize counter of function evaluations
        fcnevals=0;
        
        % perform_optimization
        perform_optimization
        
        % collect optimization results
        delta_est0 = delta_est;
        counts2    = [counts2;counts];
        perturbs2  = [perturbs2;perturb];
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
        
        coeffs   = [mcoef theta2w];
        coeffs2  = [coeffs2;coeffs(:)'];
                
        % gradients
        g           = grad_est;
        gradients   = [gradients;    g'];
        
        fprintf('\ncoeffs\n');
        printm(coeffs);
        
        fprintf('\ngradient-analytical\n');
        fprintf('%18.4f\n',g);
        
        
        toc_tmp     = toc;
        tocs        = [tocs;toc];
        
    end %perturbations loop
    
    %**********************************************************************
    % Save results
    %**********************************************************************
    cd(results_path)
    fprintf('\n');
    fprintf('Saving optimization results...\n');
    save (matfile, 'perturbs2', 'fvals', 'theta1s', 'theta2s','exit_infos',...
        'gradients','deltas' ,'gmmresids' ,'std_errors','counts2','tocs','coeffs2');
    cd(code_path)
    diary off
    
end %optimization routines loop