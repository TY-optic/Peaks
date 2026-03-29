%% sag_based_6dof_bruteforce_stitching.m
% Sag-based 6-DoF brute-force registration for overlap_cases_6dof


clear; clc; close all;

%% ========================= 路径与输出 ====================================
caseDirCandidates = { ...
    fullfile(pwd, 'overlap_cases_6dof'), ...
    fullfile(pwd, 'figures', 'overlap_cases_6dof')};

caseDir = '';
for i = 1:numel(caseDirCandidates)
    if exist(caseDirCandidates{i}, 'dir')
        caseDir = caseDirCandidates{i};
        break;
    end
end
if isempty(caseDir)
    error('未找到 overlap_cases_6dof 文件夹。');
end

outDir = fullfile(caseDir, 'sag_based_6dof_bruteforce');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

saveSingleCasePDF = true;
saveSingleCasePNG = true;

errScale    = 1e3;     % mm -> um
errUnitName = '\mum';

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');

%% ========================= 搜索参数 ======================================
cfg = struct();

cfg.search = struct();
cfg.search.nPerDim        = 5;        % 每维 5 点 (5^6 = 15625)
cfg.search.levels         = 6;        % 分层数
cfg.search.shrink         = 0.45;     % 每层缩小比例
cfg.search.minOverlapFrac = 0.12;     % 最小有效重叠比例
cfg.search.minOverlapPts  = 150;      % 最小重叠点数
cfg.search.maxSearchPts   = 3500;     % 搜索阶段点数上限
cfg.search.maxEvalPts     = inf;      % 最终评价点数上限

% 搜索跨度（以名义位姿为中心，只需覆盖扰动范围）
cfg.search.initSpan = struct();
cfg.search.initSpan.tx_mm  = 0.50;
cfg.search.initSpan.ty_mm  = 0.50;
cfg.search.initSpan.tz_mm  = 0.030;
cfg.search.initSpan.rx_deg = 0.15;
cfg.search.initSpan.ry_deg = 0.15;
cfg.search.initSpan.rz_deg = 1.50;

cfg.refine = struct();
cfg.refine.useFminsearch = true;
cfg.refine.MaxIter       = 2000;
cfg.refine.MaxFunEvals   = 10000;
cfg.refine.TolX          = 1e-12;
cfg.refine.TolFun        = 1e-15;

cfg.stats = struct();
cfg.stats.innerErodeRadius = 2;       % 内区腐蚀半径，单位像素

cfg.plot = struct();
cfg.plot.cLimPrctile = 99.5;
cfg.plot.colormap    = blue_white_red(256);

fprintf('读取目录: %s\n', caseDir);

%% ========================= 读取 case 文件 ================================
caseFiles = dir(fullfile(caseDir, 'case_*_6dof.mat'));
if isempty(caseFiles)
    error('在 %s 下未找到 case_*_6dof.mat 文件。', caseDir);
end

[~, idxSort] = sort({caseFiles.name});
caseFiles = caseFiles(idxSort);

nCase = numel(caseFiles);
fprintf('发现 %d 个 case 文件。\n', nCase);

allRes      = cell(1, nCase);
summaryCell = cell(1, nCase);

%% ========================= 主循环 =======================================
for k = 1:nCase
    fprintf('\n============================================================\n');
    fprintf('Case %d/%d : %s\n', k, nCase, caseFiles(k).name);
    fprintf('============================================================\n');

    S = load(fullfile(caseDir, caseFiles(k).name));

    T  = normalize_truth(S.truth);
    SA = normalize_subap(S.subap_A);
    SB = normalize_subap(S.subap_B);

    A = build_measured_subap(SA, cfg.search.maxSearchPts, cfg.search.maxEvalPts);
    B = build_measured_subap(SB, cfg.search.maxSearchPts, cfg.search.maxEvalPts);

    % 计算真实相对位姿
    truthRel = true_relative_pose(SA, SB);

    % 计算名义相对位姿（用于搜索初始中心）
    nominalRel = nominal_relative_pose(SA, SB);

    fprintf('A 有效点数 = %d, B 有效点数 = %d\n', A.n, B.n);
    fprintf('A 局部坐标范围: x=[%.3f, %.3f], y=[%.3f, %.3f]\n', ...
        min(A.xLocal), max(A.xLocal), min(A.yLocal), max(A.yLocal));
    fprintf('B 局部坐标范围: x=[%.3f, %.3f], y=[%.3f, %.3f]\n', ...
        min(B.xLocal), max(B.xLocal), min(B.yLocal), max(B.yLocal));

    fprintf('真实相对位姿:\n');
    fprintf('  t_true = [%+.6f, %+.6f, %+.6f] mm\n', ...
        truthRel.t(1), truthRel.t(2), truthRel.t(3));
    fprintf('  eul_true = [rx=%+.6f, ry=%+.6f, rz=%+.6f] deg\n', ...
        rad2deg(truthRel.rx), rad2deg(truthRel.ry), rad2deg(truthRel.rz));

    fprintf('名义相对位姿（搜索中心）:\n');
    fprintf('  t_nom = [%+.6f, %+.6f, %+.6f] mm\n', ...
        nominalRel.t(1), nominalRel.t(2), nominalRel.t(3));
    fprintf('  eul_nom = [rx=%+.6f, ry=%+.6f, rz=%+.6f] deg\n', ...
        rad2deg(nominalRel.rx), rad2deg(nominalRel.ry), rad2deg(nominalRel.rz));

    fprintf('名义与真实的差异:\n');
    fprintf('  dt = [%.6f, %.6f, %.6f] mm\n', ...
        abs(nominalRel.t(1)-truthRel.t(1)), ...
        abs(nominalRel.t(2)-truthRel.t(2)), ...
        abs(nominalRel.t(3)-truthRel.t(3)));
    fprintf('  dr = [%.6f, %.6f, %.6f] deg\n', ...
        abs(rad2deg(nominalRel.rx-truthRel.rx)), ...
        abs(rad2deg(nominalRel.ry-truthRel.ry)), ...
        abs(rad2deg(nominalRel.rz-truthRel.rz)));

    % 执行配准
    tic;
    reg = sag_bruteforce_6dof_register(A, B, nominalRel, cfg);
    solveTime = toc;

    % 评估位姿误差
    poseErr = evaluate_pose_error(reg, truthRel);

    fprintf('估计相对位姿:\n');
    fprintf('  t_est = [%+.6f, %+.6f, %+.6f] mm\n', ...
        reg.p(1), reg.p(2), reg.p(3));
    fprintf('  eul_est = [rx=%+.6f, ry=%+.6f, rz=%+.6f] deg\n', ...
        rad2deg(reg.p(4)), rad2deg(reg.p(5)), rad2deg(reg.p(6)));
    fprintf('位姿误差:\n');
    fprintf('  e_t      = %.6f mm\n', poseErr.et_mm);
    fprintf('  e_R      = %.6f deg\n', poseErr.eR_deg);
    fprintf('  |dtx|    = %.6f mm\n', abs(poseErr.dt_mm(1)));
    fprintf('  |dty|    = %.6f mm\n', abs(poseErr.dt_mm(2)));
    fprintf('  |dtz|    = %.6f mm\n', abs(poseErr.dt_mm(3)));
    fprintf('  |drx|    = %.6f deg\n', abs(rad2deg(poseErr.deul_rad(1))));
    fprintf('  |dry|    = %.6f deg\n', abs(rad2deg(poseErr.deul_rad(2))));
    fprintf('  |drz|    = %.6f deg\n', abs(rad2deg(poseErr.deul_rad(3))));
    fprintf('重叠区:\n');
    fprintf('  overlap ratio = %.6f\n', reg.overlapRatio);
    fprintf('  overlap RMSE  = %.6e mm\n', reg.overlapRMSE);

    % 拼接与评价
    [fusedMap, errMap, unionMask, stats, boundaryA, boundaryBestB, ...
        ZAmap, ZBmap, overlapMaskMap, innerMask] = ...
        fuse_and_evaluate_6dof(T, A, B, reg, cfg);

    fprintf('拼接误差:\n');
    fprintf('  fused RMSE        = %.6e mm (%.4f um)\n', ...
        stats.fusedRMSE, errScale*stats.fusedRMSE);
    fprintf('  fused MAE         = %.6e mm (%.4f um)\n', ...
        stats.fusedMAE, errScale*stats.fusedMAE);
    fprintf('  fused PV          = %.6e mm (%.4f um)\n', ...
        stats.fusedPV, errScale*stats.fusedPV);
    fprintf('  fused STD         = %.6e mm (%.4f um)\n', ...
        stats.fusedSTD, errScale*stats.fusedSTD);
    fprintf('  overlapOnly RMSE  = %.6e mm (%.4f um)\n', ...
        stats.overlapOnlyRMSE, errScale*stats.overlapOnlyRMSE);
    fprintf('  inner RMSE        = %.6e mm (%.4f um)\n', ...
        stats.innerRMSE, errScale*stats.innerRMSE);
    fprintf('  robust PV 99.5%%   = %.6e mm (%.4f um)\n', ...
        stats.robustPV995, errScale*stats.robustPV995);
    fprintf('耗时: %.2f s\n', solveTime);

    % 保存结果
    R = struct();
    R.name            = get_case_name(S, caseFiles(k).name);
    R.title           = get_case_title(S, caseFiles(k).name);
    R.truth           = T;
    R.subapA          = SA;
    R.subapB          = SB;
    R.A               = A;
    R.B               = B;
    R.truthRel        = truthRel;
    R.nominalRel      = nominalRel;
    R.reg             = reg;
    R.poseErr         = poseErr;
    R.fusedMap        = fusedMap;
    R.errMap          = errMap;
    R.unionMask       = unionMask;
    R.overlapMaskMap  = overlapMaskMap;
    R.innerMask       = innerMask;
    R.stats           = stats;
    R.boundaryA       = boundaryA;
    R.boundaryBestB   = boundaryBestB;
    R.ZAmap           = ZAmap;
    R.ZBmap           = ZBmap;
    R.solveTime       = solveTime;

    allRes{k} = R;

    summaryCell{k} = struct( ...
        'Case',            string(R.name), ...
        'Title',           string(R.title), ...
        'tx_est_mm',       reg.p(1), ...
        'ty_est_mm',       reg.p(2), ...
        'tz_est_mm',       reg.p(3), ...
        'rx_est_deg',      rad2deg(reg.p(4)), ...
        'ry_est_deg',      rad2deg(reg.p(5)), ...
        'rz_est_deg',      rad2deg(reg.p(6)), ...
        'tx_true_mm',      truthRel.t(1), ...
        'ty_true_mm',      truthRel.t(2), ...
        'tz_true_mm',      truthRel.t(3), ...
        'rx_true_deg',     rad2deg(truthRel.rx), ...
        'ry_true_deg',     rad2deg(truthRel.ry), ...
        'rz_true_deg',     rad2deg(truthRel.rz), ...
        'tx_nom_mm',       nominalRel.t(1), ...
        'ty_nom_mm',       nominalRel.t(2), ...
        'tz_nom_mm',       nominalRel.t(3), ...
        'rx_nom_deg',      rad2deg(nominalRel.rx), ...
        'ry_nom_deg',      rad2deg(nominalRel.ry), ...
        'rz_nom_deg',      rad2deg(nominalRel.rz), ...
        'e_t_mm',          poseErr.et_mm, ...
        'e_R_deg',         poseErr.eR_deg, ...
        'overlap_ratio',   reg.overlapRatio, ...
        'overlap_RMSE_mm', reg.overlapRMSE, ...
        'fused_RMSE_mm',   stats.fusedRMSE, ...
        'fused_MAE_mm',    stats.fusedMAE, ...
        'fused_PV_mm',     stats.fusedPV, ...
        'fused_STD_mm',    stats.fusedSTD, ...
        'overlapOnly_RMSE_mm', stats.overlapOnlyRMSE, ...
        'inner_RMSE_mm',   stats.innerRMSE, ...
        'robustPV995_mm',  stats.robustPV995, ...
        'solve_time_s',    solveTime);
end

%% ========================= 公共色标范围 =================================
allAbsErr = [];
for k = 1:nCase
    E = errScale * allRes{k}.errMap;
    M = allRes{k}.unionMask;
    allAbsErr = [allAbsErr; abs(E(M))]; 
end
allAbsErr = allAbsErr(isfinite(allAbsErr));

if isempty(allAbsErr)
    cLim = 1;
else
    cLim = prctile(allAbsErr, cfg.plot.cLimPrctile);
    if ~isfinite(cLim) || cLim <= 0
        cLim = max(allAbsErr);
    end
    if ~isfinite(cLim) || cLim <= 0
        cLim = 1;
    end
end

%% ========================= 1xN 总图 =====================================
fig = figure('Color', 'w', 'Position', [60 80 520*nCase 500]);
tl  = tiledlayout(fig, 1, nCase, 'TileSpacing', 'compact', 'Padding', 'compact');
axList = gobjects(1, nCase);

for k = 1:nCase
    R = allRes{k};
    ax = nexttile(tl, k);
    axList(k) = ax;

    Eplot = errScale * R.errMap;
    imagesc(ax, R.truth.X(1,:), R.truth.Y(:,1), Eplot, ...
        'AlphaData', double(R.unionMask));

    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    box(ax, 'on');
    hold(ax, 'on');

    colormap(ax, cfg.plot.colormap);
    clim(ax, [-cLim, cLim]);

    plot(ax, R.boundaryA(:,1),     R.boundaryA(:,2),     'k-',  'LineWidth', 2);
    plot(ax, R.boundaryBestB(:,1), R.boundaryBestB(:,2), 'k--', 'LineWidth', 2);

    xlabel(ax, 'x (mm)');
    ylabel(ax, 'y (mm)');
    title(ax, sprintf('%s\nRMSE=%.2f %s, PV=%.2f %s\ne_t=%.4f mm, e_R=%.4f^o', ...
        R.title, ...
        errScale*R.stats.fusedRMSE, errUnitName, ...
        errScale*R.stats.fusedPV,   errUnitName, ...
        R.poseErr.et_mm, R.poseErr.eR_deg), ...
        'FontWeight', 'normal', 'FontSize', 10);
    

    ax.FontSize  = 11;
    ax.LineWidth = 0.9;
    ax.TickDir   = 'in';
end

cb = colorbar(axList(end), 'eastoutside');
cb.Label.String = sprintf('Stitching error (%s)', errUnitName);
cb.FontSize = 10;

sgtitle(tl, 'Sag-based 6-DoF brute-force registration', ...
    'FontWeight', 'normal', 'FontSize', 13);

exportgraphics(fig, fullfile(outDir, 'sag_based_6dof_error_maps_1xN.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(outDir, 'sag_based_6dof_error_maps_1xN.png'), ...
    'Resolution', 300);

%% ========================= 单独导出 =====================================
if saveSingleCasePDF || saveSingleCasePNG
    for k = 1:nCase
        R = allRes{k};

        fig1 = figure('Color', 'w', 'Position', [120 120 620 540]);
        ax1  = axes(fig1);

        Eplot = errScale * R.errMap;
        imagesc(ax1, R.truth.X(1,:), R.truth.Y(:,1), Eplot, ...
            'AlphaData', double(R.unionMask));

        axis(ax1, 'image');
        set(ax1, 'YDir', 'normal');
        box(ax1, 'on');
        hold(ax1, 'on');

        colormap(ax1, cfg.plot.colormap);
        clim(ax1, [-cLim, cLim]);

        %plot(ax1, R.boundaryA(:,1),     R.boundaryA(:,2),     'k-',  'LineWidth', 2);
        %plot(ax1, R.boundaryBestB(:,1), R.boundaryBestB(:,2), 'k--', 'LineWidth', 2);

        maskA = isfinite(R.ZAmap) & R.truth.mask;
        maskB = isfinite(R.ZBmap) & R.truth.mask;

        contour(ax1, R.truth.X(1,:), R.truth.Y(:,1), double(maskA), [0.1 0.1], ...
             'LineColor', [0 0 0], 'LineStyle', '-',  'LineWidth', 2.0);

        contour(ax1, R.truth.X(1,:), R.truth.Y(:,1), double(maskB), [0.1 0.1], ...
         'LineColor', [0 0 0], 'LineStyle', '--', 'LineWidth', 2.0);


        cb1 = colorbar(ax1);
        cb1.Label.String = sprintf('Error (%s)', errUnitName);

        xlabel(ax1, 'x (mm)');
        ylabel(ax1, 'y (mm)');
        % title(ax1, sprintf(['%s\n' ...
        %     'RMSE=%.2f %s, PV=%.2f %s\n' ...
        %     'e_t=%.4f mm, e_R=%.4f^o\n' ...
        %     't=[%.4f, %.4f, %.5f] mm\n' ...
        %     'r=[%.4f, %.4f, %.4f]^o'], ...
        %     R.title, ...
        %     errScale*R.stats.fusedRMSE, errUnitName, ...
        %     errScale*R.stats.fusedPV,   errUnitName, ...
        %     R.poseErr.et_mm, R.poseErr.eR_deg, ...
        %     R.reg.p(1), R.reg.p(2), R.reg.p(3), ...
        %     rad2deg(R.reg.p(4)), rad2deg(R.reg.p(5)), rad2deg(R.reg.p(6))), ...
        %     'FontWeight', 'normal', 'FontSize', 10);


        title(ax1, sprintf(['%s\n' ...
            'RMSE=%.2f %s, PV=%.2f %s'], ...
            R.title, ...
            errScale*R.stats.fusedRMSE, errUnitName, ...
            errScale*R.stats.fusedPV,   errUnitName), ...
            'FontWeight', 'normal', 'FontSize', 10);

        ax1.FontSize  = 11;
        ax1.LineWidth = 0.9;
        ax1.TickDir   = 'in';

        if saveSingleCasePDF
            exportgraphics(fig1, fullfile(outDir, ...
                sprintf('sag6dof_error_%s.pdf', R.name)), 'ContentType', 'vector');
        end
        if saveSingleCasePNG
            exportgraphics(fig1, fullfile(outDir, ...
                sprintf('sag6dof_error_%s.png', R.name)), 'Resolution', 300);
        end

        close(fig1);
    end
end

%% ========================= 统计表与保存 ==================================
Tstat = struct2table([summaryCell{:}]);

Tstat.overlap_RMSE_um      = errScale * Tstat.overlap_RMSE_mm;
Tstat.fused_RMSE_um        = errScale * Tstat.fused_RMSE_mm;
Tstat.fused_MAE_um         = errScale * Tstat.fused_MAE_mm;
Tstat.fused_PV_um          = errScale * Tstat.fused_PV_mm;
Tstat.fused_STD_um         = errScale * Tstat.fused_STD_mm;
Tstat.overlapOnly_RMSE_um  = errScale * Tstat.overlapOnly_RMSE_mm;
Tstat.inner_RMSE_um        = errScale * Tstat.inner_RMSE_mm;
Tstat.robustPV995_um       = errScale * Tstat.robustPV995_mm;

writetable(Tstat, fullfile(outDir, 'sag_based_6dof_stats.csv'));
save(fullfile(outDir, 'sag_based_6dof_results.mat'), ...
    'allRes', 'Tstat', 'cfg', 'cLim', 'errScale', 'errUnitName');

fprintf('\n===================== 结果汇总 =====================\n');
disp(Tstat(:, {'Case', 'e_t_mm', 'e_R_deg', 'overlap_RMSE_um', ...
    'fused_RMSE_um', 'overlapOnly_RMSE_um', 'inner_RMSE_um', 'robustPV995_um'}));
fprintf('结果已保存到: %s\n', outDir);

%% ========================================================================
%%                               局部函数
%% ========================================================================

%% ---- 数据预处理 --------------------------------------------------------

function T = normalize_truth(Tin)
    T = struct();

    if isfield(Tin, 'X'),     T.X = double(Tin.X);
    elseif isfield(Tin, 'x'), T.X = double(Tin.x);
    else, error('truth 缺少 X/x 字段。');
    end

    if isfield(Tin, 'Y'),     T.Y = double(Tin.Y);
    elseif isfield(Tin, 'y'), T.Y = double(Tin.y);
    else, error('truth 缺少 Y/y 字段。');
    end

    if isfield(Tin, 'Z'),     T.Z = double(Tin.Z);
    elseif isfield(Tin, 'z'), T.Z = double(Tin.z);
    else, error('truth 缺少 Z/z 字段。');
    end

    if isfield(Tin, 'mask')
        T.mask = logical(Tin.mask);
    else
        T.mask = isfinite(T.Z);
    end
end

function S = normalize_subap(Sin)
    S = Sin;

    if isfield(S, 'c')     && iscolumn(S.c),     S.c     = S.c.';     end
    if isfield(S, 'c_nom') && iscolumn(S.c_nom), S.c_nom = S.c_nom.'; end
    if isfield(S, 'err6')  && iscolumn(S.err6),  S.err6  = S.err6.';  end

    S.X_local = double(S.X_local);
    S.Y_local = double(S.Y_local);
    S.Z_meas  = double(S.Z_meas);

    if isfield(S, 'mask')
        S.mask = logical(S.mask);
    else
        S.mask = isfinite(S.Z_meas);
    end

    if isfield(S, 'mask_local')
        S.mask_local = logical(S.mask_local);
    else
        S.mask_local = (S.X_local.^2 + S.Y_local.^2 <= S.Rsub^2);
    end

    if isfield(S, 'R_w')
        S.R_w = double(S.R_w);
    else
        error('subap 缺少 R_w 字段。');
    end

    if isfield(S, 'T_w')
        S.T_w = double(S.T_w(:));
    else
        error('subap 缺少 T_w 字段。');
    end

    if ~isfield(S, 'theta'), S.theta = 0; end
    S.theta = double(S.theta);

    if ~isfield(S, 'type'), S.type = 'circle'; end

    % 确保有名义姿态字段（用于计算名义相对位姿）
    if ~isfield(S, 'R_w_nom')
        S.R_w_nom = S.R_w;  % 如果没有名义值，用实际值
    else
        S.R_w_nom = double(S.R_w_nom);
    end
    if ~isfield(S, 'T_w_nom')
        S.T_w_nom = S.T_w;
    else
        S.T_w_nom = double(S.T_w_nom(:));
    end
end

%% ---- 构建测量子孔径结构 ------------------------------------------------

function M = build_measured_subap(Sap, maxSearchPts, maxEvalPts)
    valid = Sap.mask & isfinite(Sap.Z_meas);
    xL = Sap.X_local(valid);
    yL = Sap.Y_local(valid);
    zL = Sap.Z_meas(valid);

    if isempty(xL)
        error('子孔径有效点为空。');
    end

    % 世界坐标
    P  = Sap.R_w * [xL(:).'; yL(:).'; zL(:).'] + Sap.T_w(:);
    xW = P(1,:).';
    yW = P(2,:).';
    zW = P(3,:).';

    % 构建局部网格插值器（确保单调递增）
    xVec = double(Sap.X_local(1,:));
    yVec = double(Sap.Y_local(:,1));
    Zg   = double(Sap.Z_meas);
    Zg(~Sap.mask) = NaN;

    % 确保 yVec 严格单调递增
    if numel(yVec) > 1 && yVec(2) < yVec(1)
        yVec = flipud(yVec);
        Zg   = flipud(Zg);
    end
    % 确保 xVec 严格单调递增
    if numel(xVec) > 1 && xVec(2) < xVec(1)
        xVec = fliplr(xVec);
        Zg   = fliplr(Zg);
    end

    % 检查严格单调性
    assert(all(diff(yVec) > 0), 'yVec 不是严格单调递增');
    assert(all(diff(xVec) > 0), 'xVec 不是严格单调递增');

    Fgrid = griddedInterpolant({yVec, xVec}, Zg, 'linear', 'none');

    % 子采样索引
    idxSearch = pick_subsample_indices(numel(xL), maxSearchPts);
    idxEval   = pick_subsample_indices(numel(xL), maxEvalPts);

    M = struct();
    M.pose          = Sap;
    M.xLocal        = xL(:);
    M.yLocal        = yL(:);
    M.zLocal        = zL(:);
    M.xWorld        = xW(:);
    M.yWorld        = yW(:);
    M.zWorld        = zW(:);
    M.n             = numel(M.xLocal);
    M.Fgrid         = Fgrid;
    M.idxSearch     = idxSearch(:);
    M.idxEval       = idxEval(:);
    M.boundaryWorld = measured_boundary_world(Sap);
end

function idx = pick_subsample_indices(n, maxN)
    if isinf(maxN) || maxN >= n || maxN <= 0
        idx = (1:n).';
        return;
    end
    idx = round(linspace(1, n, maxN));
    idx = unique(max(min(idx, n), 1)).';
end

%% ---- 相对位姿计算 ------------------------------------------------------

function truthRel = true_relative_pose(SA, SB)
    % 真实相对位姿：B 相对于 A
    % PA_local = R_true * PB_local + t_true
    RA = SA.R_w;
    RB = SB.R_w;
    TA = SA.T_w(:);
    TB = SB.T_w(:);

    R_true = RA.' * RB;
    t_true = RA.' * (TB - TA);

    [rx, ry, rz] = rotm_to_eul_zyx(R_true);

    truthRel = struct();
    truthRel.R  = R_true;
    truthRel.t  = t_true(:).';
    truthRel.rx = rx;
    truthRel.ry = ry;
    truthRel.rz = rz;
end

function nomRel = nominal_relative_pose(SA, SB)
    % 名义相对位姿（不含扰动）
    RA = SA.R_w_nom;
    RB = SB.R_w_nom;
    TA = SA.T_w_nom(:);
    TB = SB.T_w_nom(:);

    R_nom = RA.' * RB;
    t_nom = RA.' * (TB - TA);

    [rx, ry, rz] = rotm_to_eul_zyx(R_nom);

    nomRel = struct();
    nomRel.R  = R_nom;
    nomRel.t  = t_nom(:).';
    nomRel.rx = rx;
    nomRel.ry = ry;
    nomRel.rz = rz;
end

%% ---- 6-DoF 强行搜索配准 ------------------------------------------------

function reg = sag_bruteforce_6dof_register(A, B, nominalRel, cfg)
    % 使用名义相对位姿作为搜索初始中心
    center = [nominalRel.t(1), nominalRel.t(2), nominalRel.t(3), ...
              nominalRel.rx,   nominalRel.ry,   nominalRel.rz];

    span = [ ...
        cfg.search.initSpan.tx_mm, ...
        cfg.search.initSpan.ty_mm, ...
        cfg.search.initSpan.tz_mm, ...
        deg2rad(cfg.search.initSpan.rx_deg), ...
        deg2rad(cfg.search.initSpan.ry_deg), ...
        deg2rad(cfg.search.initSpan.rz_deg)];

    bestP     = center;
    bestScore = inf;
    bestRMSE  = inf;
    bestOut   = empty_obj_out();

    fprintf('开始 6-DoF 强行搜索...\n');
    fprintf('搜索中心: t=[%+.4f, %+.4f, %+.4f] mm, r=[%+.4f, %+.4f, %+.4f] deg\n', ...
        center(1), center(2), center(3), ...
        rad2deg(center(4)), rad2deg(center(5)), rad2deg(center(6)));
    fprintf('初始跨度: tx=%.4f mm, ty=%.4f mm, tz=%.5f mm, rx=%.4f deg, ry=%.4f deg, rz=%.4f deg\n', ...
        span(1), span(2), span(3), ...
        rad2deg(span(4)), rad2deg(span(5)), rad2deg(span(6)));

    nPD = cfg.search.nPerDim;
    totalEvals = nPD^6 * cfg.search.levels;
    fprintf('每层 %d^6 = %d 个候选, 共 %d 层, 总计 %d 次评估\n', ...
        nPD, nPD^6, cfg.search.levels, totalEvals);

    for lev = 1:cfg.search.levels
        v1 = center(1) + linspace(-span(1), span(1), nPD);
        v2 = center(2) + linspace(-span(2), span(2), nPD);
        v3 = center(3) + linspace(-span(3), span(3), nPD);
        v4 = center(4) + linspace(-span(4), span(4), nPD);
        v5 = center(5) + linspace(-span(5), span(5), nPD);
        v6 = center(6) + linspace(-span(6), span(6), nPD);

        levelBestScore = inf;
        levelBestP     = center;
        levelBestRMSE  = inf;
        levelBestOut   = empty_obj_out();

        nEval = 0;
        for i1 = 1:numel(v1)
            for i2 = 1:numel(v2)
                for i3 = 1:numel(v3)
                    for i4 = 1:numel(v4)
                        for i5 = 1:numel(v5)
                            for i6 = 1:numel(v6)
                                p = [v1(i1), v2(i2), v3(i3), ...
                                     v4(i4), v5(i5), v6(i6)];
                                [score, rmse, out] = sag_obj_6dof( ...
                                    p, A, B, cfg, 'search');
                                nEval = nEval + 1;
                                if score < levelBestScore
                                    levelBestScore = score;
                                    levelBestP     = p;
                                    levelBestRMSE  = rmse;
                                    levelBestOut   = out;
                                end
                            end
                        end
                    end
                end
            end
        end

        center = levelBestP;
        span   = span * cfg.search.shrink;

        fprintf('  Level %d/%d (%d evals):\n', lev, cfg.search.levels, nEval);
        fprintf('    t = [%+.6f, %+.6f, %+.6f] mm\n', ...
            center(1), center(2), center(3));
        fprintf('    r = [%+.6f, %+.6f, %+.6f] deg\n', ...
            rad2deg(center(4)), rad2deg(center(5)), rad2deg(center(6)));
        fprintf('    overlap RMSE = %.6e mm, ratio = %.4f, N = %d\n', ...
            levelBestRMSE, levelBestOut.overlapRatio, levelBestOut.overlapN);

        bestP     = center;
        bestScore = levelBestScore;
        bestRMSE  = levelBestRMSE;
        bestOut   = levelBestOut;
    end

    % fminsearch 精炼
    if cfg.refine.useFminsearch
        fprintf('开始 fminsearch 精炼 (从搜索最优解出发)...\n');
        objFn = @(p) sag_obj_6dof_scalar(p, A, B, cfg);
        options = optimset('Display', 'iter', ...
            'MaxIter',    cfg.refine.MaxIter, ...
            'MaxFunEvals', cfg.refine.MaxFunEvals, ...
            'TolX',       cfg.refine.TolX, ...
            'TolFun',     cfg.refine.TolFun);

        [pOpt, fOpt] = fminsearch(objFn, bestP, options);
        [scoreOpt, rmseOpt, outOpt] = sag_obj_6dof(pOpt, A, B, cfg, 'eval');

        fprintf('fminsearch 结果:\n');
        fprintf('  t = [%+.6f, %+.6f, %+.6f] mm\n', pOpt(1), pOpt(2), pOpt(3));
        fprintf('  r = [%+.6f, %+.6f, %+.6f] deg\n', ...
            rad2deg(pOpt(4)), rad2deg(pOpt(5)), rad2deg(pOpt(6)));
        fprintf('  score = %.6e, RMSE = %.6e mm, ratio = %.4f\n', ...
            scoreOpt, rmseOpt, outOpt.overlapRatio);

        if scoreOpt < bestScore
            bestP     = pOpt;
            bestScore = scoreOpt;
            bestRMSE  = rmseOpt;
            bestOut   = outOpt;
            fprintf('  -> fminsearch 改善了结果\n');
        else
            fprintf('  -> fminsearch 未改善，保留搜索结果\n');
        end
    end

    R = eul_zyx_to_rotm(bestP(4), bestP(5), bestP(6));

    reg = struct();
    reg.p            = bestP;
    reg.R            = R;
    reg.score        = bestScore;
    reg.overlapRMSE  = bestRMSE;
    reg.overlapRatio = bestOut.overlapRatio;
    reg.overlapN     = bestOut.overlapN;
    reg.residual     = bestOut.residual;
    reg.xAov         = bestOut.xAov;
    reg.yAov         = bestOut.yAov;
    reg.zAref        = bestOut.zAref;
    reg.zAfromB      = bestOut.zAfromB;
end

%% ---- 目标函数 ----------------------------------------------------------

function f = sag_obj_6dof_scalar(p, A, B, cfg)
    [f, ~, ~] = sag_obj_6dof(p, A, B, cfg, 'eval');
end

function [score, overlapRMSE, out] = sag_obj_6dof(p, A, B, cfg, mode)
    % 选择点集
    switch lower(mode)
        case 'search'
            idx = B.idxSearch;
        case 'eval'
            idx = B.idxEval;
        otherwise
            idx = B.idxSearch;
    end

    tx = p(1); ty = p(2); tz = p(3);
    rx = p(4); ry = p(5); rz = p(6);

    R = eul_zyx_to_rotm(rx, ry, rz);

    % B 的局部坐标点
    xB = B.xLocal(idx);
    yB = B.yLocal(idx);
    zB = B.zLocal(idx);

    % 变换到 A 的局部坐标系
    PB = [xB(:).'; yB(:).'; zB(:).'];
    PA = R * PB + [tx; ty; tz];

    xAq     = PA(1,:).';   % A 局部 x
    yAq     = PA(2,:).';   % A 局部 y
    zAfromB = PA(3,:).';   % 变换后的 z（应与 A 的 sag 一致）

    % 在 A 的局部 sag 网格上插值
    zAref = A.Fgrid(yAq, xAq);

    % 有效重叠判断
    ov = isfinite(zAref) & isfinite(zAfromB);

    nOv  = nnz(ov);
    nRef = numel(idx);
    overlapRatio = nOv / max(nRef, 1);

    nMin = max(cfg.search.minOverlapPts, ...
               round(cfg.search.minOverlapFrac * nRef));
    if nOv < nMin
        score       = 1e12 + (nMin - nOv)^2;
        overlapRMSE = inf;
        out = empty_obj_out();
        out.overlapRatio = overlapRatio;
        out.overlapN     = nOv;
        return;
    end

    % 残差 = 变换后 B 的 z - A 的 sag
    res = zAfromB(ov) - zAref(ov);
    overlapRMSE = sqrt(mean(res.^2));

    % 惩罚低重叠率
    penalty = 25 * max(0, cfg.search.minOverlapFrac - overlapRatio)^2;
    score   = overlapRMSE * (1 + penalty);

    out = struct();
    out.overlapRatio = overlapRatio;
    out.overlapN     = nOv;
    out.residual     = res;
    out.xAov         = xAq(ov);
    out.yAov         = yAq(ov);
    out.zAref        = zAref(ov);
    out.zAfromB      = zAfromB(ov);
end

function out = empty_obj_out()
    out = struct();
    out.overlapRatio = 0;
    out.overlapN     = 0;
    out.residual     = [];
    out.xAov         = [];
    out.yAov         = [];
    out.zAref        = [];
    out.zAfromB      = [];
end

%% ---- 位姿误差评估 ------------------------------------------------------

function poseErr = evaluate_pose_error(reg, truthRel)
    pEst  = reg.p(:).';
    REst  = reg.R;
    RTrue = truthRel.R;

    % 平移误差
    dt = pEst(1:3) - truthRel.t(:).';
    et = norm(dt);

    % 旋转误差（测地线距离）
    dR = REst * RTrue.';
    eR = rotm_geodesic_deg(dR);

    % 欧拉角逐分量误差
    [rxE, ryE, rzE] = rotm_to_eul_zyx(REst);
    deul = [rxE - truthRel.rx, ryE - truthRel.ry, rzE - truthRel.rz];
    deul = wrap_to_pi_vec(deul);

    poseErr = struct();
    poseErr.dt_mm    = dt(:).';
    poseErr.et_mm    = et;
    poseErr.eR_deg   = eR;
    poseErr.deul_rad = deul(:).';
end

%% ---- 拼接与评价 --------------------------------------------------------

function [fusedMap, errMap, unionMask, stats, boundaryA, boundaryBestB, ...
    ZAmap, ZBmap, overlapMaskMap, innerMask] = ...
    fuse_and_evaluate_6dof(T, A, B, reg, cfg)

    % A 的世界坐标点（直接使用）
    xAw = A.xWorld;
    yAw = A.yWorld;
    zAw = A.zWorld;

    % B 的点：先用估计位姿变换到 A 局部，再用 A 姿态变换到世界
    p = reg.p(:).';
    R = reg.R;

    PB = [B.xLocal(:).'; B.yLocal(:).'; B.zLocal(:).'];
    PA = R * PB + [p(1); p(2); p(3)];

    RA = A.pose.R_w;
    TA = A.pose.T_w(:);

    PW = RA * PA + TA;
    xBw = PW(1,:).';
    yBw = PW(2,:).';
    zBw = PW(3,:).';

    % 栅格化到 truth 网格
    FAw = scatteredInterpolant(xAw, yAw, zAw, 'natural', 'none');
    FBw = scatteredInterpolant(xBw, yBw, zBw, 'natural', 'none');

    ZAmap = FAw(T.X, T.Y);
    ZBmap = FBw(T.X, T.Y);

    ZAmap(~T.mask) = NaN;
    ZBmap(~T.mask) = NaN;

    vA = isfinite(ZAmap) & T.mask;
    vB = isfinite(ZBmap) & T.mask;

    unionMask   = vA | vB;
    overlapMask = vA & vB;
    onlyA       = vA & ~vB;
    onlyB       = ~vA & vB;

    % 融合
    fusedMap = nan(size(T.Z));
    fusedMap(overlapMask) = 0.5 * (ZAmap(overlapMask) + ZBmap(overlapMask));
    fusedMap(onlyA)       = ZAmap(onlyA);
    fusedMap(onlyB)       = ZBmap(onlyB);

    % 误差 = 融合 - 真值
    errMap = nan(size(T.Z));
    errMap(unionMask) = fusedMap(unionMask) - T.Z(unionMask);

    % 统计
    ev = errMap(unionMask);
    ev = ev(isfinite(ev));

    stats = struct();
    stats.fusedRMSE = safe_rmse(ev);
    stats.fusedMAE  = safe_mae(ev);
    stats.fusedPV   = safe_pv(ev);
    stats.fusedSTD  = safe_std(ev);

    ovErr = errMap(overlapMask);
    ovErr = ovErr(isfinite(ovErr));
    stats.overlapOnlyRMSE = safe_rmse(ovErr);
    stats.overlapOnlyMAE  = safe_mae(ovErr);
    stats.overlapOnlyPV   = safe_pv(ovErr);

    innerMask = erode_mask_disk(unionMask, cfg.stats.innerErodeRadius);
    innerErr  = errMap(innerMask);
    innerErr  = innerErr(isfinite(innerErr));
    stats.innerRMSE = safe_rmse(innerErr);
    stats.innerMAE  = safe_mae(innerErr);
    stats.innerPV   = safe_pv(innerErr);

    if isempty(ev)
        stats.robustPV995 = NaN;
    else
        q1 = prctile(ev, 0.5);
        q2 = prctile(ev, 99.5);
        stats.robustPV995 = q2 - q1;
    end

    stats.unionPixels   = nnz(unionMask);
    stats.overlapPixels = nnz(overlapMask);
    stats.innerPixels   = nnz(innerMask);

    boundaryA     = A.boundaryWorld;
    boundaryBestB = transformed_boundary_world(B.pose, A.pose, reg.p);
    overlapMaskMap = overlapMask;
end

%% ---- 边界计算 ----------------------------------------------------------

function Bxy = measured_boundary_world(Sap)
    t  = linspace(0, 2*pi, 361);
    xb = Sap.Rsub * cos(t);
    yb = Sap.Rsub * sin(t);

    zb = interp2(Sap.X_local, Sap.Y_local, Sap.Z_meas, ...
                 xb, yb, 'linear', NaN);
    bad = isnan(zb);
    if any(bad)
        zb(bad) = interp2(Sap.X_local, Sap.Y_local, Sap.Z_meas, ...
                          xb(bad), yb(bad), 'nearest', 0);
    end

    P = Sap.R_w * [xb; yb; zb] + Sap.T_w(:);
    Bxy = [P(1,:).', P(2,:).'];
end

function Bxy = transformed_boundary_world(Bpose, Apose, p)
    t  = linspace(0, 2*pi, 361);
    xb = Bpose.Rsub * cos(t);
    yb = Bpose.Rsub * sin(t);

    zb = interp2(Bpose.X_local, Bpose.Y_local, Bpose.Z_meas, ...
                 xb, yb, 'linear', NaN);
    bad = isnan(zb);
    if any(bad)
        zb(bad) = interp2(Bpose.X_local, Bpose.Y_local, Bpose.Z_meas, ...
                          xb(bad), yb(bad), 'nearest', 0);
    end

    PB = [xb; yb; zb];
    R  = eul_zyx_to_rotm(p(4), p(5), p(6));
    PA = R * PB + [p(1); p(2); p(3)];

    PW = Apose.R_w * PA + Apose.T_w(:);
    Bxy = [PW(1,:).', PW(2,:).'];
end

%% ---- 旋转矩阵工具 ------------------------------------------------------

function R = eul_zyx_to_rotm(rx, ry, rz)
    % ZYX 内旋欧拉角 -> 旋转矩阵
    % R = Rz * Ry * Rx
    cx = cos(rx); sx = sin(rx);
    cy = cos(ry); sy = sin(ry);
    cz = cos(rz); sz = sin(rz);

    R = [ cz*cy,  cz*sy*sx - sz*cx,  cz*sy*cx + sz*sx; ...
          sz*cy,  sz*sy*sx + cz*cx,  sz*sy*cx - cz*sx; ...
         -sy,     cy*sx,             cy*cx            ];
end

function [rx, ry, rz] = rotm_to_eul_zyx(R)
    % 旋转矩阵 -> ZYX 内旋欧拉角
    s = -R(3,1);
    s = min(max(s, -1), 1);
    ry = asin(s);

    cy = cos(ry);
    if abs(cy) > 1e-12
        rx = atan2(R(3,2), R(3,3));
        rz = atan2(R(2,1), R(1,1));
    else
        rx = 0;
        rz = atan2(-R(1,2), R(2,2));
    end
end

function angDeg = rotm_geodesic_deg(R)
    c = 0.5 * (trace(R) - 1);
    c = min(max(c, -1), 1);
    angDeg = rad2deg(acos(c));
end

function v = wrap_to_pi_vec(v)
    v = mod(v + pi, 2*pi) - pi;
end

%% ---- 统计工具 ----------------------------------------------------------

function v = safe_rmse(x)
    if isempty(x), v = NaN; else, v = sqrt(mean(x.^2)); end
end

function v = safe_mae(x)
    if isempty(x), v = NaN; else, v = mean(abs(x)); end
end

function v = safe_pv(x)
    if isempty(x), v = NaN; else, v = max(x) - min(x); end
end

function v = safe_std(x)
    if isempty(x), v = NaN; else, v = std(x, 1); end
end

function maskOut = erode_mask_disk(maskIn, radius)
    if radius <= 0
        maskOut = maskIn;
        return;
    end
    [xx, yy] = meshgrid(-radius:radius, -radius:radius);
    ker = (xx.^2 + yy.^2 <= radius^2);
    cnt = conv2(double(maskIn), double(ker), 'same');
    maskOut = (cnt == nnz(ker));
end

%% ---- 颜色映射 ----------------------------------------------------------

function cmap = blue_white_red(m)
    if nargin < 1, m = 256; end
    n1 = floor(m/2);
    n2 = m - n1;

    c1 = [linspace(0.05, 1.00, n1)', ...
          linspace(0.25, 1.00, n1)', ...
          linspace(0.65, 1.00, n1)'];

    c2 = [linspace(1.00, 0.70, n2)', ...
          linspace(1.00, 0.05, n2)', ...
          linspace(1.00, 0.05, n2)'];

    cmap = [c1; c2];
    cmap = min(max(cmap, 0), 1);
end

%% ---- 辅助命名 ----------------------------------------------------------

function name = get_case_name(S, defaultName)
    if isfield(S, 'caseName')
        name = S.caseName;
    else
        [~, name] = fileparts(defaultName);
    end
end

function titleStr = get_case_title(S, defaultName)
    if isfield(S, 'caseTitle')
        titleStr = S.caseTitle;
    elseif isfield(S, 'etaTarget')
        titleStr = sprintf('%s, \\eta = %.0f%%', defaultName, 100*S.etaTarget);
    else
        titleStr = defaultName;
    end
end
