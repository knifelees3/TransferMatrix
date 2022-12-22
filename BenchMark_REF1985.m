% The input wavelength
lambda0=431e-9;  
% The geometry structures
layer=[
  500e-9 1.0
  850*1e-9 1.65+0.02i
  100*1e-9 1.47
  500e-9 3.5+0.08i
];
% The generated boundary conditions
nL=1;
nR=layer(4,2);
l1=layer(1,1);
l2=layer(2,1);
l3=layer(3,1);
% the dz
dz=1e-9;
%%
TEM="TM";
n_layer=layer(:,2);
num_layer=4;
z_layer=zeros(num_layer-1,1);
for l=2:num_layer-1
z_layer(l)=z_layer(l-1)+layer(l,1);
end
wavelength=lambda0;
theta=0;
d_layer=layer(2:end-1,1);
kphi=2*pi/wavelength*sin(theta);

dz=1e-9;
zmin=-500e-9;
zmax=2200e-9;
k0=2*pi/wavelength;

k_layer=k0*n_layer;
kz_layer=k_layer;

AB_Coe=CoeAB_layer_TMM(d_layer,n_layer,wavelength,kphi,TEM);
[z_layer_list,field]=field_layer_from_ABCoe(AB_Coe,d_layer,dz,zmin,zmax,kz_layer);
%%
figure()
for l=1:4
plot(z_layer_list{l}*1e6,(abs(field{l})).^2,'b-','linewidth',2);
hold on;
end
xlabel('z (nm)');
ylabel('Intensity');
set(gca,'FontName','times new roman','Fontsize',15,'XColor','k','YColor','k','LineWidth',1.3);
