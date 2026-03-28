%% ======================== 全局参数 =======================================
ds         = 0.5;          % 采样间隔 (mm)
D_full     = 100;          % 母口径直径 (mm)
R_full     = D_full / 2;
R_sub      = 30;           % 子孔径半径 (mm)

overlap_target = 0.50;
rot_deg        = 2.0;
theta_true     = deg2rad(rot_deg);
phi_layout_deg = 25;
phi_layout     = deg2rad(phi_layout_deg);

target_pv_f1   = 2.5;     

% 绘图参数
viewAz        = 36;
viewEl        = 24;
wallAlpha     = 0.10;
projLineWidth = 1.8;
surfLineWidth = 2.0;
gridStepLocal = 0.35;      % mm
nContour      = 8;
useNormalizedHeatmap = true;

%% ======================== 构造母面型 =====================================
xv = -R_full : ds : R_full;
yv = -R_full : ds : R_full;
[Xm, Ym] = meshgrid(xv, yv);
Rm = hypot(Xm, Ym);
mask_full = (Rm <= R_full);

Xbar = 6 * Xm / D_full;
Ybar = 6 * Ym / D_full;
Z_base = peaks(Xbar, Ybar);
Z_base(~mask_full) = NaN;

% 去平面项
coef0 = fit_plane_ls(Xm(mask_full), Ym(mask_full), Z_base(mask_full));
Z_base = Z_base - (coef0(1)*Xm + coef0(2)*Ym + coef0(3));
Z_base(~mask_full) = NaN;

% 归一化
vmax = max(abs(Z_base(mask_full)));
if vmax <= 0, vmax = 1; end
Z_base(mask_full) = Z_base(mask_full) / vmax;

% F1: 不平滑，直接缩放到目标 PV（mm 量级）
Z_f1 = Z_base;
coef1 = fit_plane_ls(Xm(mask_full), Ym(mask_full), Z_f1(mask_full));
Z_f1(mask_full) = Z_f1(mask_full) - (coef1(1)*Xm(mask_full) + coef1(2)*Ym(mask_full) + coef1(3));
Z_f1(~mask_full) = NaN;

pv0 = max(Z_f1(mask_full)) - min(Z_f1(mask_full));
if pv0 > 0
    Z_f1(mask_full) = Z_f1(mask_full) * (target_pv_f1 / pv0);
end

truth_X    = Xm;
truth_Y    = Ym;
truth_Z    = Z_f1;         % 单位保持为 mm
truth_mask = mask_full & isfinite(Z_f1);

fprintf('当前 Peaks 母面型 PV = %.4f mm\n', ...
    max(truth_Z(truth_mask)) - min(truth_Z(truth_mask)));

%% ======================== 子孔径布局 =====================================
sep = solve_circle_sep_from_overlap(R_sub, overlap_target);
u_layout = [cos(phi_layout); sin(phi_layout)];
cA = -0.5 * sep * u_layout;
cB =  0.5 * sep * u_layout;

%% ======================== 子孔径采样 =====================================
F1_interp = scatteredInterpolant( ...
    Xm(truth_mask), Ym(truth_mask), Z_f1(truth_mask), 'natural', 'none');

subA = sample_circular_subap_local(F1_interp, cA, 0,          R_sub, ds);
subB = sample_circular_subap_local(F1_interp, cB, theta_true, R_sub, ds);

%% ======================== 预计算绘图数据 =================================
zTruth_mm = truth_Z;
maskT     = truth_mask;

zMin_mm  = min(zTruth_mm(maskT));
zMax_mm  = max(zTruth_mm(maskT));
zSpan_mm = zMax_mm - zMin_mm;
if zSpan_mm <= 0, zSpan_mm = 1; end
zRef_mm  = zMin_mm - 0.55 * zSpan_mm;

% 子孔径边界（世界坐标）
[bxA_w, byA_w] = aperture_boundary_world(subA);
[bxB_w, byB_w] = aperture_boundary_world(subB);
bzA_w = F1_interp(bxA_w, byA_w);
bzB_w = F1_interp(bxB_w, byB_w);

% 子孔径局部网格
[XA, YA, ZA] = subap_to_local_grid(subA, gridStepLocal);
[XB, YB, ZB] = subap_to_local_grid(subB, gridStepLocal);

if useNormalizedHeatmap
    zCaseMax = max(abs([subA.z(:); subB.z(:)]));
    if zCaseMax <= 0 || ~isfinite(zCaseMax), zCaseMax = 1; end
    ZA_plot = ZA / zCaseMax;
    ZB_plot = ZB / zCaseMax;
    climHeat = [-1, 1];
    cbLabel  = 'Normalized sag';
else
    ZA_plot = ZA;          % mm
    ZB_plot = ZB;          % mm
    zCaseMax = max(abs([ZA_plot(:); ZB_plot(:)]), [], 'omitnan');
    if zCaseMax <= 0 || ~isfinite(zCaseMax), zCaseMax = 1; end
    climHeat = [-zCaseMax, zCaseMax];
    cbLabel  = 'Sag (mm)';
end

[bxA_l, byA_l] = aperture_boundary_local(subA);
[bxB_l, byB_l] = aperture_boundary_local(subB);

%% ======================== Figure 1: 三维母面型 ===========================
fig1 = figure('Color','w','Position',[80 100 920 700]);
ax1  = axes(fig1);
hold(ax1,'on'); box(ax1,'on');

step3d = 1;
hSurf = surf(ax1, truth_X(1:step3d:end,1:step3d:end), ...
                   truth_Y(1:step3d:end,1:step3d:end), ...
                   zTruth_mm(1:step3d:end,1:step3d:end), ...
    'EdgeColor','none','FaceColor','interp','FaceAlpha',1.0);
maskDown = maskT(1:step3d:end,1:step3d:end);
set(hSurf,'AlphaData',double(maskDown),'AlphaDataMapping','none');

% 子孔径边界（面上 + 投影）
plot3(ax1, bxA_w, byA_w, bzA_w, 'r-','LineWidth',surfLineWidth);
plot3(ax1, bxB_w, byB_w, bzB_w, 'b-','LineWidth',surfLineWidth);
plot3(ax1, bxA_w, byA_w, zRef_mm*ones(size(bxA_w)), 'r-','LineWidth',projLineWidth);
plot3(ax1, bxB_w, byB_w, zRef_mm*ones(size(bxB_w)), 'b-','LineWidth',projLineWidth);

% 侧壁与引导线
draw_sidewall_patch(ax1, bxA_w, byA_w, bzA_w, zRef_mm, [1 0 0], wallAlpha);
draw_sidewall_patch(ax1, bxB_w, byB_w, bzB_w, zRef_mm, [0 0 1], wallAlpha);
draw_vertical_guides(ax1, bxA_w, byA_w, bzA_w, zRef_mm, [1 0 0]);
draw_vertical_guides(ax1, bxB_w, byB_w, bzB_w, zRef_mm, [0 0 1]);

axis(ax1,'tight');
daspect(ax1,[1 1 0.06]);


view(ax1, viewAz, viewEl);
xlabel(ax1,'x (mm)');
ylabel(ax1,'y (mm)');
zlabel(ax1,'z (mm)');
grid(ax1,'on');
ax1.FontSize = 12;
title(ax1,'Peaks Freeform','FontWeight','normal','FontSize',20);

cb1 = colorbar(ax1);
cb1.Label.String = 'Sag (mm)';
cb1.FontSize = 11;

legend(ax1,{'Base surface','Sub-aperture A','Sub-aperture B'}, ...
    'Location','northeast','Box','on','FontSize',11);
camlight(ax1,'headlight');
lighting(ax1,'gouraud');

%% ======================== Figure 2: 子孔径 A 热力图 =====================
fig2 = figure('Color','w','Position',[120 120 620 540]);
ax2  = axes(fig2);
plot_subap_heatmap(ax2, XA, YA, ZA_plot, bxA_l, byA_l, climHeat, nContour, ...
    'Sub-aperture A', cbLabel);

%% ======================== Figure 3: 子孔径 B 热力图 =====================
fig3 = figure('Color','w','Position',[160 140 620 540]);
ax3  = axes(fig3);
plot_subap_heatmap(ax3, XB, YB, ZB_plot, bxB_l, byB_l, climHeat, nContour, ...
    'Sub-aperture B', cbLabel);


%% ======================== 导出图形 ==============================
% 定义目标文件夹路径
save_dir = 'figures\Peaks_Freeform';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

figs_to_save = [fig1, fig2, fig3];
axes_to_read = [ax1, ax2, ax3];

% 使用下划线代替空格，避免 Overleaf 编译报错
file_names = {'Peaks_Freeform', 'Sub-aperture_A', 'Sub-aperture_B'};

for i = 1:length(figs_to_save)
    if i == 1
        % Figure 1 (三维图): 导出为 600 DPI 高清 PNG，完美自动去白边
        full_save_path = fullfile(save_dir, sprintf('%s.png', file_names{i}));
        exportgraphics(axes_to_read(i), full_save_path, ...
            'Resolution', 600, ...
            'BackgroundColor', 'w');
        fprintf('成功导出高清 PNG: %s\n', full_save_path);
    else
        % Figure 2 & 3 (二维热力图): 继续导出为矢量 PDF
        full_save_path = fullfile(save_dir, sprintf('%s.pdf', file_names{i}));
        exportgraphics(axes_to_read(i), full_save_path, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'w');
        fprintf('成功导出矢量 PDF: %s\n', full_save_path);
    end
end



%% ======================== 保存数据到 mat 文件 ============================
% 构造 truth 结构体（母面型，未平滑）
truth_peaks = struct();
truth_peaks.X    = truth_X;
truth_peaks.Y    = truth_Y;
truth_peaks.Z    = truth_Z;
truth_peaks.mask = truth_mask;

% 构造子孔径几何结构体（仅几何参数，不含测量数据）
% 三种平滑度共用同一组子孔径几何
subap_geom_A = struct();
subap_geom_A.type  = 'circle';
subap_geom_A.c     = cA(:)';
subap_geom_A.theta = 0;
subap_geom_A.Rsub  = R_sub;

subap_geom_B = struct();
subap_geom_B.type  = 'circle';
subap_geom_B.c     = cB(:)';
subap_geom_B.theta = theta_true;
subap_geom_B.Rsub  = R_sub;

% 三组子孔径几何相同（面型差异由高斯平滑在拼接脚本中动态生成）
subap_f1_A = subap_geom_A;
subap_f1_B = subap_geom_B;


% 保存
save('subaperture_data.mat', ...
    'truth_peaks', ...
    'subap_f1_A', 'subap_f1_B', ...
    'ds', 'D_full', 'R_full', 'R_sub', ...
    'overlap_target', 'rot_deg', 'theta_true', ...
    'phi_layout_deg');

fprintf('已保存 subaperture_data.mat\n');





%% =========================================================================
%% 局部函数
%% =========================================================================

function coef = fit_plane_ls(x, y, z)
    A = [x(:), y(:), ones(numel(x),1)];
    coef = A \ z(:);
end

function d = solve_circle_sep_from_overlap(R, target_eta)
    lo = 0; hi = 2*R;
    for it = 1:80
        mid = 0.5*(lo+hi);
        eta_mid = circle_overlap_ratio(mid, R);
        if eta_mid > target_eta
            lo = mid;
        else
            hi = mid;
        end
    end
    d = 0.5*(lo+hi);
end

function eta = circle_overlap_ratio(d, R)
    if d >= 2*R, eta = 0; return; end
    if d <= 0,   eta = 1; return; end
    part1 = 2*R^2*acos(d/(2*R));
    part2 = 0.5*d*sqrt(4*R^2 - d^2);
    eta = (part1 - part2) / (pi*R^2);
end

function subap = sample_circular_subap_local(Finterp, c, theta, Rsub, ds)
    xv = -Rsub:ds:Rsub;
    yv = -Rsub:ds:Rsub;
    [Xl, Yl] = meshgrid(xv, yv);
    mask = (Xl.^2 + Yl.^2 <= Rsub^2);
    xl = Xl(mask);
    yl = Yl(mask);
    ct = cos(theta); st = sin(theta);
    xw = c(1) + ct.*xl - st.*yl;
    yw = c(2) + st.*xl + ct.*yl;
    zw = Finterp(xw, yw);
    valid = isfinite(zw);
    subap.type  = 'circle';
    subap.c     = c(:);
    subap.theta = theta;
    subap.Rsub  = Rsub;
    subap.x  = xl(valid);
    subap.y  = yl(valid);
    subap.z  = zw(valid);
end

function [bx, by] = aperture_boundary_world(S)
    t = linspace(0,2*pi,361);
    xl = S.Rsub*cos(t);
    yl = S.Rsub*sin(t);
    ct = cos(S.theta); st = sin(S.theta);
    bx = S.c(1) + ct.*xl - st.*yl;
    by = S.c(2) + st.*xl + ct.*yl;
end

function [bx, by] = aperture_boundary_local(S)
    t = linspace(0,2*pi,361);
    bx = S.Rsub*cos(t);
    by = S.Rsub*sin(t);
end

function [Xg, Yg, Zg] = subap_to_local_grid(S, ds)
    r = S.Rsub;
    xv = -r:ds:r;
    yv = -r:ds:r;
    [Xg, Yg] = meshgrid(xv, yv);
    mask = (Xg.^2 + Yg.^2 <= r^2);
    F = scatteredInterpolant(S.x, S.y, S.z, 'natural','none');
    Zg = F(Xg, Yg);
    Zg(~mask) = NaN;
end

function draw_sidewall_patch(ax, bx, by, bzTop, zBottom, faceColor, alphaVal)
    bx = bx(:); by = by(:); bzTop = bzTop(:);
    n = numel(bx);
    % 批量构建 patch 以加速（避免逐段循环）
    step = max(1, floor(n/60));  % 取约 60 段
    idx = 1:step:n;
    if idx(end) ~= n, idx(end+1) = n; end
    nSeg = numel(idx)-1;
    Xp = zeros(4, nSeg);
    Yp = zeros(4, nSeg);
    Zp = zeros(4, nSeg);
    for k = 1:nSeg
        i1 = idx(k); i2 = idx(k+1);
        Xp(:,k) = [bx(i1); bx(i2); bx(i2); bx(i1)];
        Yp(:,k) = [by(i1); by(i2); by(i2); by(i1)];
        Zp(:,k) = [zBottom; zBottom; bzTop(i2); bzTop(i1)];
    end
    patch(ax, Xp, Yp, Zp, faceColor, ...
        'FaceAlpha', alphaVal, 'EdgeColor','none');
end

function draw_vertical_guides(ax, bx, by, bzTop, zBottom, colorRGB)
    idx = round(linspace(1, numel(bx)-1, 4));
    for ii = 1:numel(idx)
        i = idx(ii);
        plot3(ax, [bx(i) bx(i)], [by(i) by(i)], [zBottom bzTop(i)], ...
            '--','Color',colorRGB,'LineWidth',1.2);
    end
end

function plot_subap_heatmap(ax, X, Y, Z, bx, by, clim, nContour, ttl, cbLabel)
    hold(ax,'on'); box(ax,'on');
    axis(ax,'equal','tight');

    valid = isfinite(Z);
    hImg = imagesc(ax, X(1,:), Y(:,1), Z);
    set(hImg,'AlphaData',double(valid));
    set(ax,'YDir','normal');
    caxis(ax, clim); %#ok<CAXIS>

    levels = linspace(clim(1), clim(2), nContour);
    [~, hc] = contour(ax, X, Y, Z, levels, 'k-','LineWidth',0.5);
    hc.Color = [0.28 0.28 0.28];

    plot(ax, bx, by, 'w-','LineWidth',1.8);

    xlabel(ax,'local x (mm)');
    ylabel(ax,'local y (mm)');
    title(ax, ttl,'FontWeight','normal','FontSize',16);
    ax.FontSize = 11;

    xlim(ax,[min(X(:)) max(X(:))]);
    ylim(ax,[min(Y(:)) max(Y(:))]);

    cb = colorbar(ax);
    cb.Label.String = cbLabel;
    cb.FontSize = 11;
end
