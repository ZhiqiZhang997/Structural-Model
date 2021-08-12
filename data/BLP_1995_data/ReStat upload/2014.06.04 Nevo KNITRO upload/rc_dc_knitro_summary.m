clear all
close all
format long
warning off all
clc

optrout      =3;
code_path    =pwd;
results_path =[code_path,'\optimization results\'];
logs_path    =[code_path,'\optimization logs\'];
add_path     =[code_path,'\optimization routines\'];
addpath(add_path);

info1_all = [];
info2_all = [];
info3_all = [];
info4_all = [];

xlsfile = [results_path,'knitro_','all','.xlsx'];

h_cell1 = cellstr({'scenario','stvalue','fval'});

h_cell2 = {'scenario','stvalue','price','brand1','brand2','brand3','brand4','brand5','brand6',...
    'brand7','brand8','brand9','brand10','brand11','brand12',...
    'brand13','brand14','brand15','brand16','brand17','brand18',...
    'brand19','brand20','brand21','brand22','brand23','brand24',...
    'const_sigma','price_sigma','sugar_sigma','mushy_sigma',...
    'const_inc','price_inc','sugar_inc','mushy_inc','price_inc2',...
    'const_age','sugar_age','mushy_age','price_child'};

for optrout=1:3
    matfile = [results_path,'knitro_',num2str(optrout),'.mat'];

    load (matfile, 'perturbs2', 'fvals', 'theta1s', 'theta2s','exit_infos',...
        'deltas' ,'gmmresids' ,'std_errors','counts2','tocs','coeffs2');

    info1_all = [info1_all;[repmat(optrout,50,1),perturbs2,fvals]];
    info2_all = [info2_all;[repmat(optrout,50,1),perturbs2,theta1s,theta2s]];
    info3_all = [info3_all;[repmat(optrout,50,1),perturbs2,std_errors]];
end

d_cell1 = num2cell(info1_all);
d_cell2 = num2cell(info2_all);
d_cell3 = num2cell(info3_all);

xlswrite(xlsfile,[h_cell1;d_cell1],'fval');
xlswrite(xlsfile,[h_cell2;d_cell2],'thetas');
xlswrite(xlsfile,[h_cell2;d_cell3],'std_errors');



