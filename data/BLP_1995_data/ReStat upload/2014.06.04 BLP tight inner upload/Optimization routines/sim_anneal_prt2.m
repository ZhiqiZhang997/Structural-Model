if(iprint >= 2);
    disp('-------------------------------------------------');
    disp('intermediate results after step length adjustment');
    disp('-------------------------------------------------');
    fprintf('new step vm            : \t');
    printm(vm');
    fprintf('current optimum theta2 : \t');
    printm(xopt');
    fprintf('current theta2         : \t');
    printm(x');
end
