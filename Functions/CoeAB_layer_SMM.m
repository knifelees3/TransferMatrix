function ABCoe = CoeAB_layer_SMM(d_layer, n_layer, wavelength, kphi, TEM, Ain, Bin)
%% The parameters descriptions:
% d_layer: thickness of each layer, excluding the upper and lower inifinite
% space, this is very important, here d_layer has only num_layer-2
% dimensions
% n_layer: refractive index of each layer
% wavelength: wavelength
% kphi: the transverse vector
% TEM: the polarizations: "TE" or "TM"
% Ain: The input amplitudes from Left
% Bin: The input amplitudes from Right
% ============================================================
% Scattering-matrix method
% Field definition:
%   ABCoe{i} = [Ai; Bi]
%   Field is evaluated at the RIGHT interface of layer i (left side)
%   Fully consistent with CoeAB_layer_TMM
% ============================================================

    num_layer = length(n_layer);

    %% ---------- Interface & propagation matrices ----------
    S_int  = cell(num_layer-1,1);   % interface i -> i+1
    S_prop = cell(num_layer-2,1);   % propagation in layer 2...N-1

    for i = 1:num_layer-1
        S_int{i} = interface_S(n_layer(i), n_layer(i+1), wavelength, kphi, TEM);
    end
    for i = 1:num_layer-2
        S_prop{i} = propagation_S(d_layer(i), n_layer(i+1), wavelength, kphi);
    end

    %% ---------- S_left_R: left boundary -> layer i RIGHT interface ----------
    S_left = cell(num_layer,1);

    curr_SL = [0,1;1,0];   % left reference plane

    % layer 1 right interface = layer 1 left interface (semi-infinite)
    S_left{1} = curr_SL;

    for i = 1:num_layer-1
        % interface i -> i+1
        curr_SL = cascade_S(curr_SL, S_int{i});

        % propagation inside layer i+1 (if finite)
        if i <= num_layer-2
            curr_SL = cascade_S(curr_SL, S_prop{i});
        end

        % now at RIGHT interface of layer i+1
        S_left{i+1} = curr_SL;
    end

    %% ---------- S_right: layer i RIGHT interface -> right boundary ----------
    S_right = cell(num_layer,1);

    curr_SR = [0,1;1,0];   % right reference plane
    S_right{num_layer} = curr_SR;

    for i = num_layer-1:-1:1
        % from right boundary back to right interface of layer i
        if i <= num_layer-2
            curr_SR = cascade_S(S_prop{i}, curr_SR);
        end
        curr_SR = cascade_S(S_int{i}, curr_SR);

        S_right{i} = curr_SR;
    end

    %% ---------- Solve fields ----------
    ABCoe = cell(num_layer,1);

    for i = 1:num_layer
        SL = S_left{i};
        SR = S_right{i};

        % Ai = SL21*Ain + SL22*Bi
        % Bi = SR11*Ai + SR12*Bin
        denom = 1 - SL(2,2)*SR(1,1);

        Ai = ( SL(2,1)*Ain + SL(2,2)*SR(1,2)*Bin ) / denom;
        Bi = ( SR(1,1)*SL(2,1)*Ain + SR(1,2)*Bin ) / denom;

        ABCoe{i} = [Ai; Bi];
    end
end
