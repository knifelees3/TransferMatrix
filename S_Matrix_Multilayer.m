function S_total = S_Matrix_Multilayer(d_list, n_layer, wl, kphi, TEM)
    S_total = interface_S(n_layer(1), n_layer(2), wl, kphi, TEM);
    for i = 2:length(n_layer)-1
        S_prop = propagation_S(d_list(i-1), n_layer(i), wl, kphi);
        S_total = cascade_S(S_total, S_prop);
        S_int = interface_S(n_layer(i), n_layer(i+1), wl, kphi, TEM);
        S_total = cascade_S(S_total, S_int);
    end
end