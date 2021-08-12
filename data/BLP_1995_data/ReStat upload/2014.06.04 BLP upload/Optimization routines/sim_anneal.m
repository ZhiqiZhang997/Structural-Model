function [xopt, fopt, nacc, nfcnev, nobds, ier, t, vm] = ...
    sim_anneal(func,x,max_sa,rt,eps,ns,nt,neps,maxevl,lb,ub,c,iprint,t,vm)

% *************************************************************************
% Set initial values
% *************************************************************************
n      = size(x,1);
xp     = zeros(n,1);
nacc   = 0;
nobds  = 0;
nfcnev = 0;
ier    = 99;
xopt   = x;
nacp   = zeros(n,1);
fstar  = 1e20*ones(neps,1);
fopts  = NaN*ones(maxevl,1);
xopts  = NaN*ones(maxevl,size(x,1));

% *************************************************************************
% Evaluate the function with input x and return value as f
% *************************************************************************
f = feval(func,x);

% *************************************************************************
%if the function is to be minimized, switch the sign of the function.
% Note that all intermediate and final output switches the sign back
%to eliminate any possible confusion for the user.
% *************************************************************************
if (max_sa == 0)
    f = -f;
end
nfcnev   = nfcnev + 1;
fopt     = f;
fstar(1) = f;
fopts(nfcnev,1) = fopt;
xopts(nfcnev,:) = x';

% *************************************************************************
%   Start the main loop. Note that it terminates if (i) the algorithm
%   succesfully optimizes the function or (ii) there are too many
%   function evaluations (more than MAXEVL).
% *************************************************************************
quit=0;
while quit==0
    nup    = 0;
    nrej   = 0;
    nnew   = 0;
    ndown  = 0;
    lnobds = 0;
    
    m=1;
    while m <= nt;
        j=1;
        while j <= ns;
            h=1;
            while h <= n;
                
                % *********************************************************
                % Generate xp, the trial value of x
                % *********************************************************
                i=1;
                while i <= n;
                    if (i == h);
                        xp(i) = x(i) + (rand(1,1)*2.- 1.) * vm(i);
                    end
                    
                    if (i ~= h);
                        xp(i) = x(i);
                    end
                    
                    % *****************************************************
                    % If xp is out of bounds,
                    % select a point in bounds for the trial
                    % *****************************************************
                    if((xp(i) < lb(i)) || (xp(i) > ub(i)));
                        xp(i) = lb(i) + (ub(i) - lb(i))*rand(1,1);
                        lnobds = lnobds + 1;
                        nobds  = nobds + 1;
                        %print some output
                        sim_anneal_prt1;
                    end
                    i = i+1;
                end
                
                % *********************************************************
                %  Evaluate the function with the trial point xp
                % *********************************************************
                fp = feval(func,xp);
                if(max_sa == 0)
                    fp = -fp;
                end
                nfcnev = nfcnev + 1;
                fopts(nfcnev,1)=fp;
                xopts(nfcnev,:)=xp';
                % *********************************************************
                % If too many function evaluations occur,
                % terminate the algorithm
                % *********************************************************
                if(nfcnev >= maxevl);
                    disp('Too many function evaluations');
                    disp('consider increasing maxevl or eps');
                    disp('or decreasing nt or rt');
                    disp('These results are likely to be poor');
                    return
                    if (max_sa == 0);
                        fopt = -fopt;
                    end
                    ier = 1;
                end
                
                % *********************************************************
                % Accept the new point if the function value increases
                % *********************************************************
                if (fp >= f);
                    x       = xp;
                    f       = fp;
                    nacc    = nacc + 1;
                    nacp(h) = nacp(h) + 1;
                    nup     = nup + 1;
                end
                
                % *********************************************************
                % If greater than any other point, record as new optimum
                % *********************************************************
                if (fp > fopt)
                    xopt = xp;
                    fopt = fp;
                    nnew = nnew + 1;
                end
                
                % *********************************************************
                %  if not moving to the appropriate direction,
                %  use the Metropolis criterion for acceptance or rejection
                % *********************************************************
                if (fp<=fopt)
                    p  = exp((fp - f)/t);
                    pp = rand(1,1);
                    
                    % Acceptance
                    if (pp < p)
                        x       = xp;
                        f       = fp;
                        nacc    = nacc + 1;
                        nacp(h) = nacp(h) + 1;
                        ndown   = ndown + 1;
                    end
                    
                    % Rejection
                    if (pp>=p)
                        nrej = nrej + 1;
                    end
                end
                h = h+1;
            end
            j = j+1;
        end
        
        % *******************************************************************
        % Adjust vm so that approximately half of all evaluations are accepted.
        % *******************************************************************
        i=1;
        while i <= n;
            ratio = nacp(i) /ns;
            if (ratio > .6);
                vm(i) = vm(i)*(1. + c(i)*(ratio - .6)/.4);
            elseif (ratio < .4);
                vm(i) = vm(i)/(1. + c(i)*((.4 - ratio)/.4));
            end
            if (vm(i) > (ub(i)-lb(i)));
                vm(i) = ub(i) - lb(i);
            end
            i = i+1;
        end
        
        %print some output
        sim_anneal_prt2;
        
        nacp = zeros(n,1);
        
        m = m+1;
    end;
    % *******************************************************************
    % print some intermediate output
    % *******************************************************************
    sim_anneal_prt3;
    
    % *******************************************************************
    % Check termination criteria
    % *******************************************************************
    quit = 0;
    fstar(1) = f;
    if ((fopt - fstar(1)) <= eps)
        quit = 1;
    end
    if ( sum(abs(f-fstar)> eps) > 0 )
        quit = 0;
    end;
    
    % *******************************************************************
    % Terminate SA if appropriate
    % *******************************************************************
    if quit==1
        x = xopt;
        ier = 0;
        if (max_sa == 0)
            fopt = -fopt;
        end
        if(iprint >= 1)
            disp('SA achieved termination criteria. ier=0');
        end
    end
    
    % *******************************************************************
    % If termination criteria are not met, prepare for another loop.
    % *******************************************************************
    %temperature update
    t = rt*t;
    i = neps;
    while i >= 2;
        fstar(i) = fstar(i-1);
        i = i-1;
    end;
    f = fopt;
    x = xopt;
    
    fopts = -fopts(1:nfcnev,:);
    xopts = xopts(1:nfcnev,:);
end
