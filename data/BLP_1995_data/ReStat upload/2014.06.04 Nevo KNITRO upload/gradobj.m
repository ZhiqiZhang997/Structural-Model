function df = gradobj(theta2)
% This function computes the gradient of the objective function

global invA IV 

load gmmresid 
load mvalold
temp = jacob(mvalold,theta2)';
df = 2*temp*IV*invA*IV'*gmmresid;

