function S_ab = cascade_S(Sa, Sb)
    denom = 1 - Sa(2,2)*Sb(1,1);
    s11 = Sa(1,1) + Sa(1,2)*Sb(1,1)*Sa(2,1) / denom;
    s12 = Sa(1,2)*Sb(1,2) / denom;
    s21 = Sb(2,1)*Sa(2,1) / denom;
    s22 = Sb(2,2) + Sb(2,1)*Sa(2,2)*Sb(1,2) / denom;
    S_ab = [s11, s12; s21, s22];
end