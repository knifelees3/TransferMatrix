% To compute the field in multi layer structures
% The material in non-magnetic
function ABCoe=CoeAB_layer_TMM(d_layer,n_layer,wavelength,kphi,TEM)

%% The parameters descriptions:
% d_layer: thickness of each layer, including the upper and lower inifinite
% space
% n_layer: refractive index of each layer
% wavelength: wavelength
% kphi: the transverse vector
% TEM: the polarizations: "TE" or "TM"

num_layer=size(n_layer(:),1);

k0=2*pi/wavelength;

k_layer=k0*n_layer;
costheta_mat=zeros(num_layer,1);
kz_layer=zeros(num_layer,1);
% Prepare the D matrix and 
D=cell(num_layer,1);
P=cell(num_layer-1,1);
M=cell(num_layer,1);
ABCoe=cell(num_layer,1);

% Initial k vector
costheta_mat(1)=sqrt(k_layer(1)^2-kphi^2)/k_layer(1);
kz_layer(1)=sqrt(k_layer(1)^2-kphi^2);
M{1}=[[1,0];[0,1]];

if TEM=="TE" % 电场方向在面外
    D{1}=[[1,1];[n_layer(1)*costheta_mat(1),-n_layer(1)*costheta_mat(1)]];
elseif TEM=="TM"% 电场方向在面内
    D{1}=[[costheta_mat(1),costheta_mat(1)];[n_layer(1),-n_layer(1)]]; 
else
    printf('Please assign the correct polarizations!')
end

%% Calculate the M matrix
for l=2:num_layer

kz_layer(l)=sqrt(k_layer(l)^2-kphi^2);

costheta_mat(l)=kz_layer(l)/k_layer(l);

if TEM=="TE"
    D{l}=[[1,1];[n_layer(l)*costheta_mat(l),-n_layer(l)*costheta_mat(l)]];
elseif TEM=="TM"
    D{l}=[[costheta_mat(l),costheta_mat(l)];[n_layer(l),-n_layer(l)]]; 
else
    printf('Please assign the correct polarizations!')
end

if l<num_layer
    phase=kz_layer(l)*(d_layer(l-1));
    if abs(imag(phase))>50
    phase=real(phase)+1i*50;
    end
    P{l-1}=[[exp(-1i*phase),0];[0,exp(1i*phase)]];
    M{l}=M{l-1}*(diag([1,1])/D{l-1})*D{l}*P{l-1};
else
M{l}=M{l-1}*(diag([1,1])/D{l-1})*D{l};
end
end

% Calculate the solution according to the boundary conditions
% A1=1,B1=0;
AN=1/(M{num_layer}(1,1));
B1=M{num_layer}(2,1)*AN;
% Obtain the coefficients
AB1=transpose([1,B1]);
for l=1:num_layer
    a=M{l}(1,1);b=M{l}(1,2);
    c=M{l}(2,1);d=M{l}(2,2);
    detM=a*d-b*c;
    if detM ~= 0
    a_inv=d/detM;
    b_inv=-b/detM;
    c_inv=-c/detM;
    d_inv=a/detM;
    invM=[[a_inv,b_inv];[c_inv,d_inv]];
    ABCoe{l}=invM*AB1;
    else
        l
     ABCoe{l}=[0,0]';% 这样写是不是合理需要去思考
%     detM=1e-13;
%     a_inv=d/detM;
%     b_inv=-b/detM;
%     c_inv=-c/detM;
%     d_inv=a/detM;
%     invM=[[a_inv,b_inv];[c_inv,d_inv]];
%     ABCoe{l}=invM*AB1;
%     invM*AB1
    end
end

end




