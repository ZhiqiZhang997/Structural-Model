% *****************************************************************
% Quasi-Newton 1
% *****************************************************************
if optrout==1
    options = optimset(...
        'GradObj','on',...
        'HessUpdate','bfgs',...
        'LargeScale','off',...
        'MaxFunEvals',mymaxfunevals,...
        'TolFun',mytolfun,...
        'TolX',mytolx,...
        'MaxIter',mymaxiters,...
        'Display','iter');
    [theta2,fval,exit_info,tmp]=fminunc('gmmobj2',theta2,options);
    counts=tmp.funcCount;
end

% *****************************************************************
% Nelder-Mead Simplex
% *****************************************************************
if optrout==2
    options=optimset('TolFun',mytolfun,...
        'TolX',mytolx,...
        'Display','iter',....
        'MaxIter',mymaxiters,...
        'MaxFunEvals',mymaxfunevals);
    [theta2,fval,exit_info,tmp] = fminsearch('gmmobj2',theta2,options);
    counts=tmp.funcCount;
end

% *****************************************************************
% Solvopt
% *****************************************************************
if optrout==3
    warning('off','MATLAB:dispatcher:InexactMatch');
    options=Soptions;
    options(2)=mytolx;
    options(3)=mytolfun;
    options(4)=mymaxiters;
    options(5)=1;
    [theta2,fval,tmp] = SOLVOPT(theta2,'gmmobj2','gradobj',options);
    exit_info=tmp(9);
    counts=tmp(10);
end

% *****************************************************************
% Conjugate gradient
% *****************************************************************
if optrout==4
    opts=[];
    opts(1)=1;
    opts(2)=1;
    opts(3)=1;
    opts(4)=mytolfun;
    opts(5)=mytolx;
    opts(6)=mymaxfunevals;
    opts(7:9)=[1e-4 1e-6 10];
    [theta2,tmp,perf,neval] = conj_grad('gmmobj3',0,theta2,opts);
    theta2=theta2(:,end);
    fval=tmp(1);
    exit_info=tmp(6);
    counts=(neval);
end

% *****************************************************************
% Quasi-Newton 2
% *****************************************************************
if optrout==5
    opts=[];
    opts(1)=1;
    opts(2)=mytolfun;
    opts(3)=mytolx;
    opts(4)=mymaxfunevals;
    [theta2, tmp, perf, D,neval] = ucminf('gmmobj3',0,theta2,opts);
    theta2=theta2(:,end);
    fval=tmp(1);
    exit_info=tmp(6);
    counts=(neval);
end

% *****************************************************************
% GA-JBES
% *****************************************************************
if optrout==6;
    uu=-10^4*ones(size(theta2));
    vv=10^4*ones(size(theta2));
    parspace=[uu,vv];
    [ga_out,ga_beta,ga_funeval]=gamin('gmmobj2',parspace,theta2);
    theta2=ga_beta;
    fval=ga_out(3);
    exit_info=1;
    counts=ga_funeval;
end;

% *****************************************************************
% Simulated annealing
% *****************************************************************
if optrout==7
    neps   = 4;
    eps    = 1.0E-4;
    rt     = .9;
    seed   = perturb*1000;
    ns_sa = 20;
    nt     = 5;
    maxevl = mymaxfunevals;
    iprint = 2;
    npar   = size(theta2,1);
    lb     = theta2-10*theta2;
    ub     = theta2+10*theta2;
    c      = 2.0*ones(npar,1);
    t      = 5.0;
    vm     = 1.0*ones(npar,1);
    randn('state',seed);
    rand('state',seed);
    fprintf('-------------------------------------------------\n');
    fprintf('              Simulated Annealing parameters     \n');
    fprintf('-------------------------------------------------\n');
    fprintf('number of parameters                   : %8.i\n',npar);
    fprintf('Temperature                            : %8.4f\n',t);
    fprintf('Temperature reduction factror          : %8.4f\n',rt);
    fprintf('epsilon                                : %8.4f\n',eps);
    fprintf('Iterations for the step to be adjusted : %8i\n',ns_sa);
    fprintf('Iterations before reducing temperature : %8i\n',nt);
    fprintf('Number of last iterations              : %8i\n',neps);
    fprintf('Max number of function evaluations     : %8i\n',maxevl);
    fprintf('Options for output printing            : %8i\n',iprint);
    fprintf('Random generator Seed                  : %8i\n',seed);
    disp('-------------------------------------------------');
    [est, fval, nacc, nfcnev, nobds, ier, t, vm] =sim_anneal('gmmobj2',...
        theta2,0,rt,eps,ns_sa,nt,neps,maxevl,lb,ub,c,iprint,t,vm);
    theta2=est;
    exit_info=1;
    fval=abs(fval);
    counts=nfcnev;
end

% *****************************************************************
% MADS
% *****************************************************************
if optrout==8
    options=psoptimset('Cache','on',...
        'CompletePoll','off',...
        'CompleteSearch','off',...
        'Display','iter',...
        'MaxFunevals',mymaxfunevals,...
        'MaxIter',mymaxiters,...
        'SearchMethod','MADSPositiveBasisNp1',...
        'PollMethod','MADSPositiveBasis2N',...
        'TolFun',mytolfun,...
        'TolX',mytolx,...
        'PlotFcn',[]);
    [theta2, fval, exit_info, tmp] =...
        patternsearch(@gmmobj2,theta2,[],[],[],[],[],[],options);
    counts=tmp.funccount;
end

% *****************************************************************
% GPS
% *****************************************************************
if optrout==9
    options=psoptimset('Cache','on',...
        'CompletePoll','off',...
        'CompleteSearch','off',...
        'Display','iter',...
        'MaxFunevals',mymaxfunevals,...
        'MaxIter',mymaxiters,...
        'SearchMethod','GPSPositiveBasisNp1',...
        'PollMethod','GPSPositiveBasis2N',...
        'TolFun',mytolfun,...
        'TolX',mytolx,...
        'PlotFcn',[]);
    [theta2, fval, exit_info, tmp] =...
        patternsearch(@gmmobj2,theta2,[],[],[],[],[],[],options);
    counts=tmp.funccount;
end

% *****************************************************************
% GA-GADS
% *****************************************************************
if optrout==10
    options=gaoptimset(...
        'TolFun'        ,mytolfun,...
        'Generations'   ,mymaxfunevals,...
        'StallGenLimit' ,50,...
        'Display'       ,'iter');
    rand('state',1000*perturb);
    randn('state',1000*perturb);
    [theta2,fval,reason,tmp,population,scores] =...
        ga(@gmmobj2,size(theta2,1),[],[],[],[],[],[],[],options);
    counts=tmp.funccount;
    theta2=theta2';
    exit_info=1;
end

% *****************************************************************
% Simulated Annealing - GADS
% *****************************************************************
if optrout==11
    options=saoptimset(...
        'AnnealingFcn',@annealingboltz,...
        'TolFun'      ,mytolfun,...
        'MaxIter'     ,mymaxiters,...
        'MaxFunEvals' ,mymaxfunevals,...
        'InitialTemperature',5,...
        'ReannealInterval'  ,100,...
        'TemperatureFcn'    ,@temperatureboltz,...
        'Display'           ,'iter');
    
    rand('state',1000*perturb);
    randn('state',1000*perturb);
    
    [theta2,fval,exit_flag,tmp] =...
        simulannealbnd(@gmmobj2,theta2,[],[],options);
    exit_info=1;
    counts=tmp.funccount;
end

