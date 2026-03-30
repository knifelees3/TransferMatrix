clear; clc; close all;
addpath('../Functions');
%% ===== 基本参数 =====
lambda0 = 431e-9;
TEM = "TM";
theta = 0;

layer = [
  500e-9     1.0
  850e-9     1.65 + 0.02i
  100e-9     1.47
  500e-9     3.5  + 0.08i
];

n_layer = layer(:,2);
d_layer = layer(2:end-1,1);
num_layer = length(n_layer);

Ain = 1;
Bin = 0;

%% ===== 波矢 =====
k0 = 2*pi/lambda0;
kphi = k0*sin(theta);

k_layer = k0 * n_layer;
kz_layer = sqrt(k_layer.^2 - kphi^2);

%% ===== TMM：求 ABCoe =====
ABCoe = CoeAB_layer_TMM( ...
    d_layer, n_layer, lambda0, kphi, TEM, Ain, Bin);

%% ===== 空间网格 =====
zmin = -500e-9;
zmax = 2200e-9;
dz   = 1e-9;
z_vec = zmin:dz:zmax;

%% ===== 用你已有的函数算场 =====
Ez = Get_Field_From_ABCoe(z_vec, ABCoe, d_layer, kz_layer);

%% ===== 作图 =====
figure;
plot(z_vec*1e9, abs(Ez).^2, 'b-', 'LineWidth',2);
xlabel('z (nm)');
ylabel('|E|^2');
set(gca,'FontName','Times New Roman', ...
        'FontSize',15, ...
        'LineWidth',1.3, ...
        'TickDir','out');
grid on;

%% ============================================================
%  保存图片
% ============================================================
print(gcf,'BenchMark_Mack1985.png','-dpng','-r300');
