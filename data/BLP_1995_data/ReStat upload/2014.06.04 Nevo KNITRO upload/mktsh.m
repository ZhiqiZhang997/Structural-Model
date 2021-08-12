function f = mktsh(mval, expmu)
% This function computes the market share for each product

global ns 
f = sum((ind_sh(mval,expmu))')/ns;
f = f';