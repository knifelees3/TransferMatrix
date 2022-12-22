n1=3.5;
n2=1.45;

nm=1e-9;
WL0=1300*nm;num_swep=601;
WLMat=linspace(600,2200,num_swep)*nm;
d1=WL0/4/n1;
d2=WL0/4/n2;

% The geometry structures
layer=[
  500e-9 1.0
  d1 n1
  d2 n2
  d1 n1
  d2 n2
  d1 n1
  d2 n2
  d1 n1
  d2 n2
  d1 n1
  d2 n2
  500e-9 1
];

TEM="TM";
n_layer=layer(:,2);
num_layer=12;
z_layer=zeros(num_layer-1,1);
for l=2:num_layer-1
z_layer(l)=z_layer(l-1)+layer(l,1);
end
theta=0;
d_layer=layer(2:end-1,1);
kphi=2*pi./WLMat*sin(theta);

RT=zeros(num_swep,1);

for l=1:num_swep
    AB_Coe=CoeAB_layer_TMM(d_layer,n_layer,WLMat(l),kphi(l),"TE");
    RT(l,1)=AB_Coe{1}(2);
    RT(l,2)=AB_Coe{end}(1);
end

%%
figure()
plot(WLMat,abs(RT(:,1)).^2,'r');
hold on
plot(WLMat,abs(RT(:,2)).^2,'b');
hold on
hold on
plot(WLMat,abs(RT(:,1)).^2+abs(RT(:,2)).^2,'k');
legend('Reflections','Transmission');
xlabel('Wavelength (nm)');