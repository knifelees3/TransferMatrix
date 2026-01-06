clear; clc; close all;

%% ===== 1. 结构参数 =====
c = 1;
n0 = 1.0;
ns = 2.0;
ds = 2.0;
dg = 1.0;
N_pairs = 5;

TEM  = "TE";
kphi = 0;

% 层结构：air | (ns, n0)* | air
d_layer = [repmat([ds, dg], 1, N_pairs-1), ds];
n_layer = [n0, repmat([ns, n0], 1, N_pairs-1), ns, n0];

num_layer = length(n_layer);

%% ===== 2. 频率点（选一个共振附近）=====
omega = 3.0;
wavelength = 2*pi / omega;

k0 = omega / c;
kz_layer = zeros(num_layer,1);
for l = 1:num_layer
    kz_layer(l) = sqrt((k0*n_layer(l))^2 - kphi^2);
end

%% ===== 3. 计算 ABCoe =====
Ain = 1; Bin = 0;

[ABCoe_TMM,~,~] = CoeAB_layer_TMM( ...
    d_layer, n_layer, wavelength, kphi, TEM, Ain, Bin);

ABCoe_SMM = CoeAB_layer_SMM( ...
    d_layer, n_layer, wavelength, kphi, TEM, Ain, Bin);

%% ===== 4. 构造空间坐标 =====
Ltot = sum(d_layer);
z_vec = linspace(-1, Ltot+1, 4000);

%% ===== 5. 重构空间场 =====
Ez_TMM = Get_Field_From_ABCoe(z_vec, ABCoe_TMM, d_layer, kz_layer);
Ez_SMM = Get_Field_From_ABCoe(z_vec, ABCoe_SMM, d_layer, kz_layer);

%% ===== 6. 误差评估 =====
err = max(abs(Ez_TMM - Ez_SMM));
fprintf('Max |Ez_TMM - Ez_SMM| = %.3e\n', err);

%% ===== 7. 可视化 =====
figure('Color','w','Position',[100,100,1000,400]);

subplot(1,2,1);
plot(z_vec, real(Ez_TMM), 'b-', 'LineWidth',1.5); hold on;
plot(z_vec, real(Ez_SMM), 'r--', 'LineWidth',1.2);
xlabel('z'); ylabel('Re(E)');
title('Real(E(z))');
legend('TMM','SMM'); grid on;

subplot(1,2,2);
plot(z_vec, abs(Ez_TMM), 'b-', 'LineWidth',1.5); hold on;
plot(z_vec, abs(Ez_SMM), 'r--', 'LineWidth',1.2);
xlabel('z'); ylabel('|E|');
title('|E(z)|');
legend('TMM','SMM'); grid on;
%% ============================================================
% 6. 保存图片
% ============================================================
print(gcf,'BraggMirror_TMM_vs_SMM_FieldDistribution.png','-dpng','-r300');