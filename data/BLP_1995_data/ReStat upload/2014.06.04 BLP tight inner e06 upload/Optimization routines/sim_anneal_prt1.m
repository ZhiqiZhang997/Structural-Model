if(iprint >= 3);
    disp('');
    disp('current x');
    disp(x)
    if (max_sa==1)
        disp('current f');
        disp(f);
    end
    if (max_sa==0)
        disp('current f');
        disp(-f)
    end
    disp('trial x');
    disp(xp)
    disp('point rejected since out of bounds');
end