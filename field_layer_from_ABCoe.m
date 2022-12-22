% 对于给定的坐标范围以及系数，计算E(0,0,z)
function [z_layer_list,field]=field_layer_from_ABCoe(AB_Coe,d_layer,dz,zmin,zmax,kz_layer)
num_layer=size(d_layer(:),1)+2;
z_layer=cell(num_layer-1,1);
%%

field={};
z_layer_list={};

z_layer{1}=0;
for l=2:num_layer-1
z_layer{l}=z_layer{l-1}+d_layer(l-1);
end

%% Obtain the field coefficients
count=0;


l=1;
kz=kz_layer(l);
% (1) (zmin,l,l+1,zmax)
if zmin<=z_layer{l} && zmax>=z_layer{l}
count=count+1;
z_layer_list{count}=zmin:dz:z_layer{l};
field{count}=AB_Coe{l}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l}))+...
    AB_Coe{l}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l}));
elseif zmin<=z_layer{l} && zmax<z_layer{l}
count=count+1;
z_layer_list{count}=zmin:dz:zmax;
field{count}=AB_Coe{l}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l}))+...
    AB_Coe{l}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l}));   
end

%% l in gap

for l=1:2
    kz=kz_layer(l+1);
    % (1) (l,zmin,zmax,l+1)
    if zmin>z_layer{l} && zmin<=z_layer{l+1} && zmax>=z_layer{l} && zmax<z_layer{l+1}
    count=count+1;
    z_layer_list{count}=zmin:dz:zmax;
    field{count}=AB_Coe{l+1}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l}))+...
        AB_Coe{l+1}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l}));
    end
    % (2) (zmin,l,l+1,zmax)
    if zmin<=z_layer{l} && zmin<=z_layer{l+1} && zmax>=z_layer{l} && zmax>=z_layer{l+1}
    count=count+1;
    z_layer_list{count}=z_layer{l}:dz:z_layer{l+1};
    field{count}=AB_Coe{l+1}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l+1}))+...
        AB_Coe{l+1}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l+1}));   

    end
    % (3) (zmin,l,zmax,l+1)
    if zmin<=z_layer{l} && zmin<=z_layer{l+1} && zmax>=z_layer{l} && zmax<z_layer{l+1}
    count=count+1;
    z_layer_list{count}=z_layer{l}:dz:zmax;
    field{count}=AB_Coe{l+1}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l+1}))+...
        AB_Coe{l+1}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l+1}));
    end


    % (4) (l,zmin,l+1,zmax)
    if zmin>z_layer{l} && zmin<z_layer{l+1} && zmax>=z_layer{l} && zmax>=z_layer{l+1}
    count=count+1;
    z_layer_list{count}=zmin:dz:z_layer{l+1};
    field{count}=AB_Coe{l+1}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l+1}))+...
        AB_Coe{l+1}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l+1}));
    end
end
%%
% l in final
l=num_layer-1;
kz=kz_layer(l+1);
if zmin>z_layer{l}
    count=count+1;
z_layer_list{count}=zmin:dz:zmax;
field{count}=AB_Coe{l+1}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l}))+...
    AB_Coe{l+1}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l}));
end

if zmin<=z_layer{l} && zmax>=z_layer{l}
    count=count+1;
z_layer_list{count}=z_layer{l}:dz:zmax;
field{count}=AB_Coe{l+1}(1)*exp(1i*kz*(z_layer_list{count}-z_layer{l}))+...
    AB_Coe{l+1}(2)*exp(-1i*kz*(z_layer_list{count}-z_layer{l}));
end
end