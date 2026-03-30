
% To compute the field in multi layer structures
% The material in non-magnetic
function [ABCoe,M,Mcas]=CoeAB_layer_TMM(d_layer,n_layer,wavelength,kphi,TEM, Ain, Bin)

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

num_layer=size(n_layer(:),1);

k0=2*pi/wavelength;

k_layer=k0*n_layer;
costheta_mat=zeros(num_layer,1);
kz_layer=zeros(num_layer,1);

% Prepare the D matrix and 
D=cell(num_layer,1);
P=cell(num_layer-2,1);
M=cell(num_layer,1);
ABCoe=cell(num_layer,1);

% k vector in each layers
for l=1:num_layer
costheta_mat(l)=sqrt(k_layer(l)^2-kphi^2)/k_layer(l);
kz_layer(l)=sqrt(k_layer(l)^2-kphi^2);
end

% initial the first D matrix
if TEM=="TE" % 电场方向在面外
    D{1}=[[1,1];[n_layer(1)*costheta_mat(1),-n_layer(1)*costheta_mat(1)]];
elseif TEM=="TM"% 电场方向在面内
    D{1}=[[costheta_mat(1),costheta_mat(1)];[n_layer(1),-n_layer(1)]]; 
else
    printf('Please assign the correct polarizations!')
end
M{1}=[[1,0];[0,1]];% The first Matrix

%% Calculate the M matrix
% we have num_layer-1 Mayrix
for l=1:num_layer-1

    ld=l+1;
if TEM=="TE"
    D{ld}=[[1,1];[n_layer(ld)*costheta_mat(ld),-n_layer(ld)*costheta_mat(ld)]];
elseif TEM=="TM"
    D{ld}=[[costheta_mat(ld),costheta_mat(ld)];[n_layer(ld),-n_layer(ld)]]; 
else
    printf('Please assign the correct polarizations!')
end

if l<num_layer-1
    phase=kz_layer(ld)*(d_layer(l)); % 注意这里，之前写的是l-1,是有问题的。
    if abs(imag(phase))>50
    phase=real(phase)+1i*50;
    end
    P{l}=[[exp(-1i*phase),0];[0,exp(1i*phase)]];
    M{ld}=M{l}*(diag([1,1])/D{l})*D{ld}*P{l};
    Mcas{ld}=(diag([1,1])/D{l})*D{ld}*P{l};
else
M{ld}=M{l}*(diag([1,1])/D{l})*D{ld};
Mcas{ld}=(diag([1,1])/D{l})*D{ld};
end
end

% Calculate the solution according to the boundary conditions
Mtot = M{num_layer};
A1 = Ain;
BNp1 = Bin; %left most layer left propogating
B1 = -((Mtot(1,2)*Mtot(2,1)-Mtot(1,1)*Mtot(2,2))*BNp1 - Mtot(2,1)*A1) / Mtot(1,1);
AB1 = [A1; B1];
ANp1 = (A1 - Mtot(1,2)*BNp1) / Mtot(1,1);%right most layer right propogating
ABCoe{end}=[ANp1;BNp1];

for l=1:num_layer-1
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

    %ABCoe{l} = M{l} \ AB1;
    else
    ABCoe{l}=[0,0]';% 这样写是不是合理需要去思考
    % detM=1e-10;
    % a_inv=d/detM;
    % b_inv=-b/detM;
    % c_inv=-c/detM;
    % d_inv=a/detM;
    % invM=[[a_inv,b_inv];[c_inv,d_inv]];
    % ABCoe{l}=invM*AB1;
    % %invM*AB1
    end
end

end






