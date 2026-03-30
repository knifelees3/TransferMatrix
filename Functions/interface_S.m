function S = interface_S(n1, n2, wl, kphi, TEM)
    k0=2*pi/wl;
    kz1 = sqrt((k0*n1)^2 - kphi^2);
    kz2 = sqrt((k0*n2)^2 - kphi^2);
    
    if TEM == "TE"
        % TE模式
        r = (kz1 - kz2) / (kz1 + kz2);
        t = 2*kz1 / (kz1 + kz2);
        tp = 2*kz2 / (kz1 + kz2);
    elseif TEM == "TM"
        % TM模式
        r = (n2^2*kz1 - n1^2*kz2) / (n2^2*kz1 + n1^2*kz2);
        t = (2*n1*n2*kz1) / (n2^2*kz1 + n1^2*kz2);
        tp = (2*n1*n2*kz2) / (n2^2*kz1 + n1^2*kz2);
    else
        error('TEM must be "TE" or "TM"');
    end
    
    rp = -r;  % 反向传播的反射系数
    
    S = [r, tp; t, rp];
end