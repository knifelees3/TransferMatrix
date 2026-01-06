function S = propagation_S(d, n, wl, kphi)
    k0=2*pi/wl;
    kz = sqrt((k0*n)^2 - kphi^2);
    P = exp(1i * kz * d);
    % 传播矩阵：s11=s22=0 (无反射)，s12=s21=P (双向相位延迟)
    S = [0, P; P, 0];
end