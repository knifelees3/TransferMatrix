clear; clc; close all;

%% ============================================================
% 1. 物理结构参数
% ============================================================
c  = 1;
ns = 3; 
n0 = 1;

ds = 2;
dg = 1;

N_pairs = 20;              % Bragg 周期数
TEM  = "TE";
kphi = 0;

d_layer = [repmat([ds, dg], 1, N_pairs-1), ds];
n_layer = [1, repmat([ns, n0], 1, N_pairs-1), ns, 1];

num_layer = length(n_layer);

%% ============================================================
% 2. 频率扫描
% ============================================================
w_scan = linspace(39, 40, 1000);

T_TMM = zeros(size(w_scan));
R_TMM = zeros(size(w_scan));
T_SMM = zeros(size(w_scan));
R_SMM = zeros(size(w_scan));

fprintf('Computing Bragg mirror spectrum (TMM vs SMM)...\n');

%% ============================================================
% 3. 主循环：TMM / SMM 并行计算
% ============================================================
for iw = 1:length(w_scan)

    omega = w_scan(iw);
    wavelength = 2*pi / omega;

    Ain = 1;
    Bin = 0;

    % ------------------ TMM ------------------
    ABCoe_TMM = CoeAB_layer_TMM( ...
        d_layer, n_layer, wavelength, kphi, TEM, Ain, Bin);

    r_TMM = ABCoe_TMM{1}(2);        % B1
    t_TMM = ABCoe_TMM{end}(1);      % A_{N+1}

    T_TMM(iw) = abs(t_TMM)^2;
    R_TMM(iw) = abs(r_TMM)^2;

    % ------------------ SMM ------------------
    ABCoe_SMM = CoeAB_layer_SMM( ...
        d_layer, n_layer, wavelength, kphi, TEM, Ain, Bin);

    r_SMM = ABCoe_SMM{1}(2);
    t_SMM = ABCoe_SMM{end}(1);

    T_SMM(iw) = abs(t_SMM)^2;
    R_SMM(iw) = abs(r_SMM)^2;
end

%% ============================================================
% 4. 数值一致性与能量守恒检查
% ============================================================
fprintf('Max |T_TMM - T_SMM| = %.3e\n', max(abs(T_TMM - T_SMM)));
fprintf('Max |R_TMM - R_SMM| = %.3e\n', max(abs(R_TMM - R_SMM)));
fprintf('Max |T+R-1| (TMM)  = %.3e\n', max(abs(T_TMM + R_TMM - 1)));

%% ============================================================
% 5. 可视化
% ============================================================
col_T = [0.10, 0.45, 0.85];
col_R = [0.85, 0.33, 0.10];

figure('Color','w','Position',[100,100,900,550]); hold on;

% Band gap shading
T_threshold = 0.05;
is_gap = T_TMM < T_threshold;
gap_idx = find(is_gap);

if ~isempty(gap_idx)
    d_idx = diff(gap_idx);
    breaks = [0, find(d_idx > 1), length(gap_idx)];
    for k = 1:length(breaks)-1
        idx = gap_idx(breaks(k)+1 : breaks(k+1));
        patch([w_scan(idx(1)) w_scan(idx(end)) ...
               w_scan(idx(end)) w_scan(idx(1))], ...
              [0 0 1.05 1.05], ...
              [0.85 0.85 0.85], ...
              'FaceAlpha',0.4,'EdgeColor','none');
    end
end

% Spectra
plot(w_scan, T_TMM, 'Color',col_T,'LineWidth',4);
plot(w_scan, T_SMM, '--','Color',col_T,'LineWidth',2.0,'HandleVisibility','off');

plot(w_scan, R_TMM, 'Color',col_R,'LineWidth',4);
plot(w_scan, R_SMM, '--','Color',col_R,'LineWidth',2.0,'HandleVisibility','off');

xlabel('\omega','FontSize',22);
ylabel('Intensity','FontSize',22);
title('Bragg Mirror Spectrum: TMM vs SMM','FontSize',22);

legend({'Band gap','Transmission','Reflection'}, ...
    'Location','best','Box','off','FontSize',16);

set(gca,'FontSize',20,'LineWidth',1.4,'TickDir','out');
ylim([0 1.05]);
grid on;

%% ============================================================
% 6. 保存图片
% ============================================================
print(gcf,'BraggMirror_TMM_vs_SMM.png','-dpng','-r300');
