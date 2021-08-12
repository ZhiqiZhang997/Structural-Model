function [f,g] = gmmobj3(theta2,par)

global invA theta1 x1 IV fcnevals delta gmmresid grad_est delta_est

if size(theta2,1)==1
    theta2=theta2';
end

par       = 0;
delta     = meanval(theta2);
delta_est = delta;

if max(isnan(delta)) == 1
    f=1e+10;
    g=-999999*ones(size(theta2));
    gmmresid=1e+10*ones(size(delta));
    fcnevals=fcnevals+1;
else
	temp1    = x1'*IV;
	temp2    = delta'*IV;
    theta1   = inv(temp1*invA*temp1')*temp1*invA*temp2';
	gmmresid = delta - x1*theta1;
	temp1    = gmmresid'*IV;
	f        = temp1*invA*temp1';
    fcnevals = fcnevals+1;
    g        = gradobj(theta2);
    grad_est = g;
end

