%% build_overlap_cases_6dof_and_save_three_mats.m
% 说明：
% 1) 图 1 仍采用二维标称子孔径布局示意，仅用于论文展示；
% 2) 真正传递到下一层计算的数据是带 6-DoF 扰动重采样后的子孔径 Z_meas；
% 3) 默认保存 3 个独立 MAT 文件，每个文件对应一个 overlap case；
% 4) 默认每个 case 使用不同的 6-DoF 随机误差；若需三组 case 共用同一组误差，请将
%    useSameErrorForAllCases = true。
%
% 输出：
%   case_O1_high_overlap_6dof.mat
%   case_O2_medium_overlap_6dof.mat
%   case_O3_low_overlap_6dof.mat
%   Fig_overlap_1x3.pdf / png
%   单独 3 张二维示意图 pdf / png

clear; clc; close all;

%% --------------------------- 输入与输出 ---------------------------------
matFile = 'subaperture_data.mat';
outDir  = fullfile('figures', 'overlap_cases_6dof');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

etaList    = [0.70, 0.50, 0.30];
caseNames  = {'O1_high_overlap', 'O2_medium_overlap', 'O3_low_overlap'};
caseTitles = {'(a) High overlap', '(b) Medium overlap', '(c) Low overlap'};

exportSingleFigure = true;

%% --------------------------- 面型预处理参数 -----------------------------
useSmoothedTruth = true;
sigmaTruth       = 5;
removePiston     = true;

%% --------------------------- 6-DoF 误差参数 -----------------------------
% 长度单位: mm
% 角度单位: degree（内部自动转成 rad）
noise_level = struct();
noise_level.t_xy = 0.50;   % x/y 平移标准差
noise_level.t_z  = 0.005;  % z 向活塞标准差
noise_level.r_xy = 0.020;  % rx/ry 标准差
noise_level.r_z  = 0.50;   % rz 标准差

% 是否三个 case 共用同一组误差
useSameErrorForAllCases = false;

% 是否固定随机种子
useFixedSeed = true;
rngSeed      = 20260324;
if useFixedSeed
    rng(rngSeed, 'twister');
else
    rng('shuffle');
end

%% --------------------------- 验证条件 -----------------------------------
validationRule = struct();
validationRule.minCoverage = 0.85;   % 有效采样覆盖率阈值
validationRule.overlapTol  = 5e-4;   % 名义 overlap 与目标值允许偏差
validationRule.saveOnlyIfPass = false;

%% --------------------------- 读取母面型 ---------------------------------
S  = load(matFile);
T0 = normalize_truth_surface(pick_base_truth_surface(S));

if isfield(S, 'ds')
    ds = double(S.ds);
else
    ds = infer_grid_spacing(T0);
end

if useSmoothedTruth
    Zuse = masked_gaussian_smooth(T0.Z, T0.mask, sigmaTruth);
else
    Zuse = T0.Z;
end

if removePiston
    Zuse = Zuse - mean(Zuse(T0.mask & isfinite(Zuse)), 'omitnan');
end

Tfixed      = T0;
Tfixed.Z    = Zuse;
Tfixed.ds   = ds;

%% --------------------------- 参考几何读取 -------------------------------
[subapA_ref, subapB_ref, hasRef] = try_get_reference_subaps(S);

if hasRef
    Rsub   = get_rsub(subapA_ref, []);
    thetaA = get_theta(subapA_ref);
    thetaB = get_theta(subapB_ref);

    dirVec = subapB_ref.c(:) - subapA_ref.c(:);
    if norm(dirVec) < eps
        dirVec = [1; 0];
    end
    dirVec = dirVec / norm(dirVec);
else
    if isfield(S, 'R_sub')
        Rsub = double(S.R_sub);
    else
        error('未找到参考子孔径，且数据中缺少 R_sub。');
    end

    thetaA = 0.0;
    if isfield(S, 'theta_true')
        thetaB = double(S.theta_true);
    else
        thetaB = 0.0;
    end
    dirVec = [1; 0];
end

fprintf('R_sub   = %.4f mm\n', Rsub);
fprintf('theta_A = %.4f deg\n', rad2deg(thetaA));
fprintf('theta_B = %.4f deg\n', rad2deg(thetaB));
fprintf('dirVec  = [%.6f, %.6f]\n\n', dirVec(1), dirVec(2));

%% --------------------------- 误差生成策略 -------------------------------
if useSameErrorForAllCases
    shared_err6_A = draw_random_6dof(noise_level);
    shared_err6_B = draw_random_6dof(noise_level);
else
    shared_err6_A = [];
    shared_err6_B = [];
end

%% --------------------------- 构造三个 case ------------------------------
cases = cell(1, numel(etaList));

for k = 1:numel(etaList)
    eta = etaList(k);
    d   = solve_separation_from_overlap_ratio(Rsub, eta);

    cA_nom = (-0.5 * d) * dirVec;
    cB_nom = (+0.5 * d) * dirVec;

    if useSameErrorForAllCases
        err6_A = shared_err6_A;
        err6_B = shared_err6_B;
    else
        err6_A = draw_random_6dof(noise_level);
        err6_B = draw_random_6dof(noise_level);
    end

    subapA = sample_subap_with_6dof_error(Tfixed, cA_nom, thetaA, Rsub, ds, err6_A);
    subapB = sample_subap_with_6dof_error(Tfixed, cB_nom, thetaB, Rsub, ds, err6_B);

    caseK = struct();
    caseK.name              = caseNames{k};
    caseK.title             = caseTitles{k};
    caseK.truth             = Tfixed;
    caseK.subapA            = subapA;
    caseK.subapB            = subapB;
    caseK.etaTarget         = eta;
    caseK.etaNominal        = circle_overlap_ratio_equal_radius(Rsub, d);
    caseK.separationNominal = d;
    caseK.err6_A            = err6_A;
    caseK.err6_B            = err6_B;
    caseK.meta = struct();
    caseK.meta.useSmoothedTruth       = useSmoothedTruth;
    caseK.meta.sigmaTruth             = sigmaTruth;
    caseK.meta.removePiston           = removePiston;
    caseK.meta.noise_level            = noise_level;
    caseK.meta.useSameErrorForAllCases= useSameErrorForAllCases;
    caseK.meta.ds                     = ds;
    caseK.meta.description = ['2D aperture layout is schematic only; ', ...
                              'sub-aperture data passed downstream are generated ', ...
                              'from 6-DoF perturbed resampling.'];

    caseK.validation = validate_case(caseK, validationRule);
    cases{k} = caseK;

    fprintf('============================================================\n');
    fprintf('%s\n', caseK.name);
    fprintf('eta_target  = %.4f\n', caseK.etaTarget);
    fprintf('eta_nominal = %.4f\n', caseK.etaNominal);
    fprintf('d_nominal   = %.4f mm\n', caseK.separationNominal);
    fprintf('A: tx=%+.4f mm, ty=%+.4f mm, tz=%+.5f mm, rx=%+.5f deg, ry=%+.5f deg, rz=%+.5f deg\n', ...
        err6_A(1), err6_A(2), err6_A(3), rad2deg(err6_A(4)), rad2deg(err6_A(5)), rad2deg(err6_A(6)));
    fprintf('B: tx=%+.4f mm, ty=%+.4f mm, tz=%+.5f mm, rx=%+.5f deg, ry=%+.5f deg, rz=%+.5f deg\n', ...
        err6_B(1), err6_B(2), err6_B(3), rad2deg(err6_B(4)), rad2deg(err6_B(5)), rad2deg(err6_B(6)));
    fprintf('coverage A  = %.4f\n', caseK.validation.coverageA);
    fprintf('coverage B  = %.4f\n', caseK.validation.coverageB);
    fprintf('validation  = %d\n', caseK.validation.pass);
    fprintf('============================================================\n\n');
end

%% --------------------------- 保存三个 MAT 文件 --------------------------
for k = 1:numel(cases)
    caseData   = cases{k};
    truth      = caseData.truth;         %#ok<NASGU>
    subap_A    = caseData.subapA;        %#ok<NASGU>
    subap_B    = caseData.subapB;        %#ok<NASGU>
    meta       = caseData.meta;          %#ok<NASGU>
    validation = caseData.validation;    %#ok<NASGU>

    etaTarget         = caseData.etaTarget;          %#ok<NASGU>
    etaNominal        = caseData.etaNominal;         %#ok<NASGU>
    separationNominal = caseData.separationNominal;  %#ok<NASGU>
    err6_A            = caseData.err6_A;             %#ok<NASGU>
    err6_B            = caseData.err6_B;             %#ok<NASGU>
    caseName          = caseData.name;               %#ok<NASGU>
    caseTitle         = caseData.title;              %#ok<NASGU>

    saveFile = fullfile(outDir, sprintf('case_%s_6dof.mat', caseData.name));

    if validationRule.saveOnlyIfPass
        if validation.pass
            save(saveFile, ...
                'truth', 'subap_A', 'subap_B', ...
                'meta', 'validation', ...
                'etaTarget', 'etaNominal', 'separationNominal', ...
                'err6_A', 'err6_B', 'caseName', 'caseTitle');
            fprintf('已保存: %s\n', saveFile);
        else
            warning('未通过验证，跳过保存: %s', saveFile);
        end
    else
        save(saveFile, ...
            'truth', 'subap_A', 'subap_B', ...
            'meta', 'validation', ...
            'etaTarget', 'etaNominal', 'separationNominal', ...
            'err6_A', 'err6_B', 'caseName', 'caseTitle');
        fprintf('已保存: %s\n', saveFile);
    end
end

%% --------------------------- 二维 1x3 示意图 ---------------------------
Zplot = Tfixed.Z;
mask  = Tfixed.mask & isfinite(Zplot);

zMin  = min(Zplot(mask), [], 'omitnan');
zMax  = max(Zplot(mask), [], 'omitnan');
Znorm = (Zplot - zMin) ./ max(zMax - zMin, eps);
Znorm(~mask) = NaN;

fig = figure('Color', 'w', 'Position', [80 80 1600 500]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

axList = gobjects(1, numel(cases));

for k = 1:numel(cases)
    ax = nexttile(tl, k);
    axList(k) = ax;

    imagesc(ax, Tfixed.X(1,:), Tfixed.Y(:,1), Znorm, 'AlphaData', double(mask));
    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    hold(ax, 'on');
    box(ax, 'on');

    colormap(ax, turbo);

    bA = aperture_boundary_world_from_subap(cases{k}.subapA);
    bB = aperture_boundary_world_from_subap(cases{k}.subapB);

    plot(ax, bA(:,1), bA(:,2), '-', 'Color', [1.00, 0.00, 0.00], 'LineWidth', 1.7);
    plot(ax, bB(:,1), bB(:,2), '-', 'Color', [0.00, 0.10, 1.00], 'LineWidth', 1.7);

    title(ax, sprintf('%s\n\\eta = %.0f%%,  d = %.2f mm', ...
        cases{k}.title, 100*cases{k}.etaNominal, cases{k}.separationNominal), ...
        'FontWeight', 'normal', 'FontSize', 12);

    xlabel(ax, 'x (mm)', 'FontSize', 11);
    ylabel(ax, 'y (mm)', 'FontSize', 11);

    ax.FontName  = 'Times New Roman';
    ax.FontSize  = 11;
    ax.LineWidth = 1.0;
    ax.TickDir   = 'in';
end

cb = colorbar(axList(end), 'eastoutside');
ylabel(cb, 'normalized sag', 'FontSize', 11);
cb.FontName  = 'Times New Roman';
cb.FontSize  = 10;
cb.LineWidth = 0.8;

sgtitle(tl, 'Representative surface with different overlap ratios', ...
    'FontWeight', 'normal', 'FontSize', 14, 'FontName', 'Times New Roman');

exportgraphics(fig, fullfile(outDir, 'Fig_overlap_1x3.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outDir, 'Fig_overlap_1x3.png'), 'Resolution', 400);

%% --------------------------- 导出单图 -----------------------------------
if exportSingleFigure
    for k = 1:numel(cases)
        fig1 = figure('Color', 'w', 'Position', [120 120 620 560]);
        ax1  = axes(fig1);

        imagesc(ax1, Tfixed.X(1,:), Tfixed.Y(:,1), Znorm, 'AlphaData', double(mask));
        axis(ax1, 'image');
        set(ax1, 'YDir', 'normal');
        hold(ax1, 'on');
        box(ax1, 'on');

        colormap(ax1, turbo);

        bA = aperture_boundary_world_from_subap(cases{k}.subapA);
        bB = aperture_boundary_world_from_subap(cases{k}.subapB);

        plot(ax1, bA(:,1), bA(:,2), '-', 'Color', [1.00, 0.00, 0.00], 'LineWidth', 1.8);
        plot(ax1, bB(:,1), bB(:,2), '-', 'Color', [0.00, 0.10, 1.00], 'LineWidth', 1.8);

        title(ax1, sprintf('%s\n\\eta = %.0f%%,  d = %.2f mm', ...
            cases{k}.title, 100*cases{k}.etaNominal, cases{k}.separationNominal), ...
            'FontWeight', 'normal', 'FontSize', 12);

        xlabel(ax1, 'x (mm)', 'FontSize', 11);
        ylabel(ax1, 'y (mm)', 'FontSize', 11);

        ax1.FontName  = 'Times New Roman';
        ax1.FontSize  = 11;
        ax1.LineWidth = 1.0;
        ax1.TickDir   = 'in';

        cb1 = colorbar(ax1, 'eastoutside');
        ylabel(cb1, 'normalized sag', 'FontSize', 11);
        cb1.FontName  = 'Times New Roman';
        cb1.FontSize  = 10;
        cb1.LineWidth = 0.8;

        exportgraphics(fig1, fullfile(outDir, sprintf('%s_layout.pdf', cases{k}.name)), 'ContentType', 'vector');
        exportgraphics(fig1, fullfile(outDir, sprintf('%s_layout.png', cases{k}.name)), 'Resolution', 400);
        close(fig1);
    end
end

fprintf('\n全部完成。\n');

%% ============================== 局部函数 ================================

function err6 = draw_random_6dof(noise_level)
    err6 = [ ...
        noise_level.t_xy * randn(1,1), ...
        noise_level.t_xy * randn(1,1), ...
        noise_level.t_z  * randn(1,1), ...
        deg2rad(noise_level.r_xy * randn(1,1)), ...
        deg2rad(noise_level.r_xy * randn(1,1)), ...
        deg2rad(noise_level.r_z  * randn(1,1))  ...
    ];
end

function subap = sample_subap_with_6dof_error(truth, c_nom, theta_nom, Rsub, ds, err6)
    % 近轴近似下的 6-DoF 扰动重采样
    % c_nom, theta_nom 仅作为标称显示与名义布局；
    % 真正传递到下游的数据为带 6-DoF 扰动后的 Z_meas 与 mask。

    tx = err6(1); ty = err6(2); tz = err6(3);
    rx = err6(4); ry = err6(5); rz = err6(6);

    c_nom = c_nom(:);

    xv = -Rsub:ds:Rsub;
    yv = -Rsub:ds:Rsub;
    [Xl, Yl] = meshgrid(xv, yv);
    mask_local = (Xl.^2 + Yl.^2 <= Rsub^2);

    R_nom = [cos(theta_nom), -sin(theta_nom), 0; ...
             sin(theta_nom),  cos(theta_nom), 0; ...
             0,               0,              1];

    Rx = [1, 0, 0; ...
          0, cos(rx), -sin(rx); ...
          0, sin(rx),  cos(rx)];

    Ry = [ cos(ry), 0, sin(ry); ...
           0,       1, 0; ...
          -sin(ry), 0, cos(ry)];

    Rz = [cos(rz), -sin(rz), 0; ...
          sin(rz),  cos(rz), 0; ...
          0,        0,       1];

    R_err = Rz * Ry * Rx;
    R_w   = R_err * R_nom;
    T_w   = [c_nom(1) + tx; c_nom(2) + ty; tz];

    P0x = R_w(1,1) .* Xl + R_w(1,2) .* Yl + T_w(1);
    P0y = R_w(2,1) .* Xl + R_w(2,2) .* Yl + T_w(2);
    P0z = R_w(3,1) .* Xl + R_w(3,2) .* Yl + T_w(3);

    Zsurf = interp2(truth.X, truth.Y, truth.Z, P0x, P0y, 'linear', NaN);
    Msurf = interp2(truth.X, truth.Y, double(truth.mask), P0x, P0y, 'nearest', 0) > 0.5;

    rayDir = R_w(:,3);
    vz = rayDir(3);
    if abs(vz) < 1e-12
        error('vz 过小，当前 6-DoF 姿态不适合该近轴模型。');
    end

    Z_meas = (Zsurf - P0z) ./ vz;
    validMask = mask_local & Msurf & isfinite(Zsurf) & isfinite(Z_meas);

    Z_meas(~validMask) = NaN;

    subap = struct();
    subap.type      = 'circle';
    subap.c         = c_nom(:).';                    % 标称二维显示中心
    subap.theta     = theta_nom;                     % 标称二维显示角度
    subap.c_nom     = c_nom(:).';
    subap.theta_nom = theta_nom;
    subap.c_meas_xy = [c_nom(1)+tx, c_nom(2)+ty];   % 仅平面投影近似
    subap.theta_meas_xy = wrap_to_pi_local(theta_nom + rz);

    subap.Rsub   = Rsub;
    subap.ds     = ds;
    subap.err6   = err6(:).';
    subap.R_w    = R_w;
    subap.T_w    = T_w;
    subap.rayDir = rayDir;

    subap.X_local    = Xl;
    subap.Y_local    = Yl;
    subap.Z_meas     = Z_meas;
    subap.mask_local = mask_local;
    subap.mask       = validMask;
end

function V = validate_case(C, rule)
    A = C.subapA;
    B = C.subapB;

    V = struct();

    V.sizeMatchA = isequal(size(A.X_local), size(A.Y_local), size(A.Z_meas), size(A.mask), size(A.mask_local));
    V.sizeMatchB = isequal(size(B.X_local), size(B.Y_local), size(B.Z_meas), size(B.mask), size(B.mask_local));

    V.coverageA = nnz(A.mask) / max(nnz(A.mask_local), 1);
    V.coverageB = nnz(B.mask) / max(nnz(B.mask_local), 1);

    V.passCoverageA = (V.coverageA >= rule.minCoverage);
    V.passCoverageB = (V.coverageB >= rule.minCoverage);

    dNom = norm(B.c(:) - A.c(:));
    V.nominalSeparationCheck = dNom;
    V.nominalOverlapCheck    = circle_overlap_ratio_equal_radius(A.Rsub, dNom);
    V.passOverlap            = abs(V.nominalOverlapCheck - C.etaTarget) <= rule.overlapTol;

    V.passFiniteErrA = all(isfinite(A.err6));
    V.passFiniteErrB = all(isfinite(B.err6));

    V.passFiniteZ_A = all(isfinite(A.Z_meas(A.mask)));
    V.passFiniteZ_B = all(isfinite(B.Z_meas(B.mask)));

    V.pass = all([ ...
        V.sizeMatchA, V.sizeMatchB, ...
        V.passCoverageA, V.passCoverageB, ...
        V.passOverlap, ...
        V.passFiniteErrA, V.passFiniteErrB, ...
        V.passFiniteZ_A, V.passFiniteZ_B]);

    V.rule = rule;
end

function T = pick_base_truth_surface(S)
    fn = fieldnames(S);

    priorityList = {'truth_peaks', 'truth_f1', 'truth_F1', 'truth'};
    for i = 1:numel(priorityList)
        if isfield(S, priorityList{i})
            T = S.(priorityList{i});
            return;
        end
    end

    idx = find(startsWith(fn, 'truth_'), 1, 'first');
    if ~isempty(idx)
        T = S.(fn{idx});
        return;
    end

    error('未找到可用母面型字段（如 truth_peaks / truth_f1 / truth_*）。');
end

function T = normalize_truth_surface(Tin)
    T = struct();

    if isfield(Tin, 'X')
        T.X = double(Tin.X);
    elseif isfield(Tin, 'x')
        T.X = double(Tin.x);
    else
        error('truth surface 缺少 X/x 字段。');
    end

    if isfield(Tin, 'Y')
        T.Y = double(Tin.Y);
    elseif isfield(Tin, 'y')
        T.Y = double(Tin.y);
    else
        error('truth surface 缺少 Y/y 字段。');
    end

    if isfield(Tin, 'Z')
        T.Z = double(Tin.Z);
    elseif isfield(Tin, 'z')
        T.Z = double(Tin.z);
    else
        error('truth surface 缺少 Z/z 字段。');
    end

    if isfield(Tin, 'mask')
        T.mask = logical(Tin.mask);
    else
        T.mask = isfinite(T.Z);
    end
end

function ds = infer_grid_spacing(T)
    dx = [];
    dy = [];

    if size(T.X, 2) >= 2
        dx = median(abs(diff(T.X(1,:))), 'omitnan');
    end
    if size(T.Y, 1) >= 2
        dy = median(abs(diff(T.Y(:,1))), 'omitnan');
    end

    cand = [dx, dy];
    cand = cand(isfinite(cand) & cand > 0);

    if isempty(cand)
        error('无法从 truth surface 自动推断网格间隔 ds。');
    end

    ds = cand(1);
end

function Zs = masked_gaussian_smooth(Z, mask, sigmaPix)
    if nargin < 3 || sigmaPix <= 0
        Zs = Z;
        return;
    end

    mask = logical(mask) & isfinite(Z);
    Z0 = Z;
    Z0(~mask) = 0;

    fsz = 2 * ceil(3 * sigmaPix) + 1;

    num = imgaussfilt(Z0, sigmaPix, 'FilterSize', fsz, 'Padding', 'replicate');
    den = imgaussfilt(double(mask), sigmaPix, 'FilterSize', fsz, 'Padding', 'replicate');

    Zs = num ./ max(den, eps);
    Zs(~mask) = NaN;
end

function [subapA, subapB, ok] = try_get_reference_subaps(S)
    ok = false;
    subapA = [];
    subapB = [];

    candA = {'subap_A', 'subap_f1_A', 'subap_F1_A'};
    candB = {'subap_B', 'subap_f1_B', 'subap_F1_B'};

    for i = 1:numel(candA)
        if isfield(S, candA{i})
            subapA = S.(candA{i});
            break;
        end
    end

    for i = 1:numel(candB)
        if isfield(S, candB{i})
            subapB = S.(candB{i});
            break;
        end
    end

    if ~isempty(subapA) && ~isempty(subapB)
        if iscolumn(subapA.c), subapA.c = subapA.c.'; end
        if iscolumn(subapB.c), subapB.c = subapB.c.'; end
        ok = true;
    end
end

function th = get_theta(Sap)
    if isfield(Sap, 'theta')
        th = double(Sap.theta);
    else
        th = 0.0;
    end
end

function Rsub = get_rsub(Sap, defaultRsub)
    if isfield(Sap, 'Rsub')
        Rsub = double(Sap.Rsub);
    else
        if isempty(defaultRsub)
            error('参考子孔径中缺少 Rsub，且未提供默认值。');
        end
        Rsub = double(defaultRsub);
    end
end

function d = solve_separation_from_overlap_ratio(R, eta)
    if eta <= 0
        d = 2 * R;
        return;
    end
    if eta >= 1
        d = 0;
        return;
    end

    fun = @(x) circle_overlap_ratio_equal_radius(R, x) - eta;
    d = fzero(fun, [0, 2*R]);
end

function eta = circle_overlap_ratio_equal_radius(R, d)
    d = max(min(d, 2*R), 0);

    if d >= 2*R
        eta = 0;
        return;
    end

    if d <= 0
        eta = 1;
        return;
    end

    Aint = 2 * R^2 * acos(d / (2*R)) - 0.5 * d * sqrt(max(4*R^2 - d^2, 0));
    eta  = Aint / (pi * R^2);
end

function Bxy = aperture_boundary_world_from_subap(S)
    if ~isfield(S, 'theta')
        S.theta = 0;
    end

    switch lower(S.type)
        case 'circle'
            t  = linspace(0, 2*pi, 361);
            xl = S.Rsub * cos(t);
            yl = S.Rsub * sin(t);

        case 'rect'
            xl = [-S.Lx/2,  S.Lx/2,  S.Lx/2, -S.Lx/2, -S.Lx/2];
            yl = [-S.Ly/2, -S.Ly/2,  S.Ly/2,  S.Ly/2, -S.Ly/2];

        otherwise
            error('未知子孔径类型: %s', S.type);
    end

    ct = cos(S.theta);
    st = sin(S.theta);

    xw = S.c(1) + ct .* xl - st .* yl;
    yw = S.c(2) + st .* xl + ct .* yl;

    Bxy = [xw(:), yw(:)];
end

function ang = wrap_to_pi_local(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end