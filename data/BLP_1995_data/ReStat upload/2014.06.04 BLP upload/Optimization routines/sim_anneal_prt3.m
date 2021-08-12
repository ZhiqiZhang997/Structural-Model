if(iprint >= 1)
    totmov = nup+ndown+nrej;
    fprintf('------------------------------------------------------\n');
    fprintf('intermediate results before next temperature reduction\n');
    fprintf('------------------------------------------------------\n');    
    
    if max_sa==1
        fprintf('temperature                :%18.3f\n',t);
        fprintf('function value             :%18.4f\n',abs(fopt));
        fprintf('total moves                :%18.0f\n',totmov);
        fprintf('downhill                   :%18.0f\n',nup);
        fprintf('accepted uphill            :%18.0f\n',ndown);
        fprintf('rejected uphill            :%18.0f\n',nrej);
        fprintf('trials out of bounds       :%18.0f\n',lnobds);
        fprintf('new minima at this temp    :%18.0f\n',nnew);
    end
    if max_sa==0
        fprintf('temperature                :%18.3f\n',t);
        fprintf('function value             :%18.4f\n',abs(fopt));
        fprintf('total moves                :%18.0f\n',totmov);
        fprintf('downhill                   :%18.0f\n',nup);
        fprintf('accepted uphill            :%18.0f\n',ndown);
        fprintf('rejected uphill            :%18.0f\n',nrej);
        fprintf('trials out of bounds       :%18.0f\n',lnobds);
        fprintf('new minima at this temp    :%18.0f\n',nnew);
    end

end
