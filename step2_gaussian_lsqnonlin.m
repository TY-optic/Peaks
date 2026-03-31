%% step2_gaussian_lsqnonlin_refine_2d6dof.m
% 基于 Gaussian-curvature step2 初值的 lsqnonlin 三维位姿精配准
% 评价方式与 sag_based_6dof_bruteforce_stitching.m 统一：
%   1) 真值相对位姿：由 R_w / T_w 计算
%   2) 位姿误差：evaluate_pose_error
%   3) 最终拼接误差：fuse_and_evaluate_6dof (truth 网格 + natural 插值)
%
% 输入：
%   results\step2_gaussian_curvature_2d_6dof\step2_gaussian_all_cases.mat
%
% 输出：
%   results\step3_gaussian_lsqnonlin_refine_2d6dof\step3_result_<case>.mat
%   results\step3_gaussian_lsqnonlin_refine_2d6dof\step3_gaussian_all_cases.mat
%   results\step3_gaussian_lsqnonlin_refine_2d6dof\step3_gaussian_summary.csv

clear; clc;

%% ========================= 路径与参数 ===================================
step2File = fullfile('results', 'step2_gaussian_curvature_2d_6dof', ...
    'step2_gaussian_all_cases.mat');

if ~exist(step2File, 'file')
    error('未找到 step2 输出文件: %s', step2File);
end

if exist('lsqnonlin', 'file') ~= 2
    error('当前环境未检测到 lsqnonlin，请确认已安装 Optimization Toolbox。');
end

outDir = fullfile('results', 'step3_gaussian_lsqnonlin_refine_2d6dof');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

params = struct();

params.caseMode = 'all';          % 'all' / 'O1' / 'O2' / 'O3' / 完整 caseName

% 候选点
params.minCandidatePts    = 800;
params.useStep2OverlapOnly = true;

% A 面连续模型
params.normalSmoothSigmaPix = 1.0;
params.maskInterpThresh     = 0.55;

% 优化权重
params.lambdaMaskPenalty = 2.0e-2;
params.lambdaRxRy        = 2.0e-3;
params.lambdaTz          = 2.0e-3;

% 优化边界（相对 step2 初值）
params.bound.rxDeg = 1.0;
params.bound.ryDeg = 1.0;
params.bound.rzDeg = 2.0;
params.bound.txMM  = 1.5;
params.bound.tyMM  = 1.5;
params.bound.tzMM  = 0.2;

% lsqnonlin 选项
params.opt = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'SpecifyObjectiveGradient', false, ...
    'MaxIterations', 80, ...
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12);

% 与 sag 统一的评价配置
cfgEval = struct();
cfgEval.search = struct();
cfgEval.search.minOverlapFrac = 0.12;
cfgEval.search.minOverlapPts  = 150;
cfgEval.stats = struct();
cfgEval.stats.innerErodeRadius = 2;

% 输出控制
params.savePerCaseMat = true;
params.saveSummaryCsv = true;

errScale = 1e3;   % mm -> um

%% ========================= 读取 step2 结果 ===============================
S = load(step2File);
if ~isfield(S, 'allResults')
    error('step2 文件中缺少 allResults。');
end

allResults = S.allResults;
allResults = filter_cases(allResults, params.caseMode);

fprintf('Using step2 file: %s\n', step2File);
fprintf('Detected %d case(s).\n', numel(allResults));

allStep3 = cell(1, numel(allResults));
summaryCell = cell(1, numel(allResults));

%% ========================= 主循环 =======================================
for kCase = 1:numel(allResults)
    R2 = allResults{kCase};

    fprintf('\n============================================================\n');
    fprintf('Running step3-lsqnonlin case: %s\n', R2.caseName);
    fprintf('============================================================\n');

    % -------- 数据规范化 --------
    T  = normalize_truth(R2.truth);
    SA = normalize_subap_from_step2(R2.subapA);
    SB = normalize_subap_from_step2(R2.subapB);

    % -------- 构建与 sag 统一的测量子孔径结构 --------
    A = build_measured_subap_from_step2(SA, R2.localA, inf, inf);
    B = build_measured_subap_from_step2(SB, R2.localB, inf, inf);

    % -------- 真值相对位姿（统一口径）--------
    truthRel = true_relative_pose(SA, SB);

    % -------- A 面连续模型 --------
    modelA = build_surface_model_A(R2, params);

    % -------- step2 初值 -> 内部参数 pose0 = [rx ry rz tx ty tz] --------
    p0 = build_initial_pose_from_step2(R2.fine);

    % -------- 从 step2 overlap 构造固定 B 候选点 --------
    [PBcand, metaCand] = build_candidate_points_from_step2(R2, p0, params);

    if size(PBcand, 1) < params.minCandidatePts
        error('Case %s 候选点数不足: %d < %d', ...
            R2.caseName, size(PBcand,1), params.minCandidatePts);
    end

    % -------- 初值残差 --------
    stat0 = evaluate_pose_against_surface(PBcand, p0, modelA, params);

    % -------- 优化边界 --------
    [lb, ub] = build_pose_bounds(p0, params);

    % -------- lsqnonlin --------
    fun = @(p) residual_pose_point2surface(p, PBcand, modelA, params);

    tic;
    [pEst, ~, residual, exitflag, output] = lsqnonlin(fun, p0, lb, ub, params.opt); %#ok<ASGLU>
    solveTime = toc;

    % -------- 优化后残差 --------
    stat1 = evaluate_pose_against_surface(PBcand, pEst, modelA, params);

    % -------- 转换为 sag 兼容 reg 结构 --------
    reg = gaussian_pose_to_reg(pEst);

    % -------- 统一的重叠区误差评价 --------
    reg = attach_overlap_metrics_to_reg(reg, A, B, cfgEval);

    % -------- 统一的位姿误差 --------
    poseErr = evaluate_pose_error(reg, truthRel);

    % -------- 统一的拼接与误差评价 --------
    [fusedMap, errMap, unionMask, stats, boundaryA, boundaryBestB, ...
        ZAmap, ZBmap, overlapMaskMap, innerMask] = ...
        fuse_and_evaluate_6dof(T, A, B, reg, cfgEval);

    % -------- 结果打包 --------
    R = struct();
    R.name           = get_case_name_from_step2(R2);
    R.title          = get_case_title_from_step2(R2);
    R.caseName       = R2.caseName;
    R.caseTitle      = R.title;
    if isfield(R2, 'etaTarget'),  R.etaTarget  = R2.etaTarget;  end
    if isfield(R2, 'etaNominal'), R.etaNominal = R2.etaNominal; end

    R.truth          = T;
    R.subapA         = SA;
    R.subapB         = SB;
    R.localA         = R2.localA;
    R.localB         = R2.localB;
    R.A              = A;
    R.B              = B;

    R.step2          = R2.fine;
    R.step2Truth2D   = R2.truth2d;

    R.truthRel       = truthRel;
    R.reg            = reg;
    R.poseErr        = poseErr;

    R.initialPoseInternal   = p0;
    R.estimatedPoseInternal = pEst;

    [R_init, t_init] = posevec_to_rt_internal(p0);
    [R_est,  t_est]  = posevec_to_rt_internal(pEst);
    R.initialTransformInternal   = struct('R', R_init, 't', t_init);
    R.estimatedTransformInternal = struct('R', R_est,  't', t_est);

    R.initialResidualStat = stat0;
    R.finalResidualStat   = stat1;

    R.surfaceModelA    = rmfield(modelA, {'Zfill','dZdx','dZdy','mask'});
    R.candidateMeta    = metaCand;
    R.numCandidatePts  = size(PBcand,1);

    R.fusedMap         = fusedMap;
    R.errMap           = errMap;
    R.unionMask        = unionMask;
    R.overlapMaskMap   = overlapMaskMap;
    R.innerMask        = innerMask;
    R.stats            = stats;
    R.boundaryA        = boundaryA;
    R.boundaryBestB    = boundaryBestB;
    R.ZAmap            = ZAmap;
    R.ZBmap            = ZBmap;
    R.solveTime        = solveTime;
    R.lsqnonlin        = struct( ...
        'exitflag', exitflag, ...
        'output', output, ...
        'resnorm', sum(residual.^2));

    allStep3{kCase} = R;

    summaryCell{kCase} = struct( ...
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
        'initPlane_RMSE_mm', stat0.rmsePlane, ...
        'finalPlane_RMSE_mm', stat1.rmsePlane, ...
        'initDz_RMSE_mm',  stat0.rmseZ, ...
        'finalDz_RMSE_mm', stat1.rmseZ, ...
        'solve_time_s',    solveTime, ...
        'exitflag',        exitflag);

    if params.savePerCaseMat
        step3Result = R; %#ok<NASGU>
        save(fullfile(outDir, sprintf('step3_result_%s.mat', R.caseName)), 'step3Result');
    end

    fprintf('Initial valid points        = %d\n', stat0.numValid);
    fprintf('Final   valid points        = %d\n', stat1.numValid);
    fprintf('Initial RMSE plane residual = %.6e mm\n', stat0.rmsePlane);
    fprintf('Final   RMSE plane residual = %.6e mm\n', stat1.rmsePlane);
    fprintf('Initial RMSE dz             = %.6e mm\n', stat0.rmseZ);
    fprintf('Final   RMSE dz             = %.6e mm\n', stat1.rmseZ);

    fprintf('统一位姿误差:\n');
    fprintf('  e_t      = %.6f mm\n', poseErr.et_mm);
    fprintf('  e_R      = %.6f deg\n', poseErr.eR_deg);
    fprintf('  |dtx|    = %.6f mm\n', abs(poseErr.dt_mm(1)));
    fprintf('  |dty|    = %.6f mm\n', abs(poseErr.dt_mm(2)));
    fprintf('  |dtz|    = %.6f mm\n', abs(poseErr.dt_mm(3)));
    fprintf('  |drx|    = %.6f deg\n', abs(rad2deg(poseErr.deul_rad(1))));
    fprintf('  |dry|    = %.6f deg\n', abs(rad2deg(poseErr.deul_rad(2))));
    fprintf('  |drz|    = %.6f deg\n', abs(rad2deg(poseErr.deul_rad(3))));

    fprintf('统一重叠区评价:\n');
    fprintf('  overlap ratio = %.6f\n', reg.overlapRatio);
    fprintf('  overlap RMSE  = %.6e mm\n', reg.overlapRMSE);

    fprintf('统一拼接误差:\n');
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
end

%% ========================= 汇总保存 =====================================
Tstat = struct2table([summaryCell{:}]);

Tstat.overlap_RMSE_um      = errScale * Tstat.overlap_RMSE_mm;
Tstat.fused_RMSE_um        = errScale * Tstat.fused_RMSE_mm;
Tstat.fused_MAE_um         = errScale * Tstat.fused_MAE_mm;
Tstat.fused_PV_um          = errScale * Tstat.fused_PV_mm;
Tstat.fused_STD_um         = errScale * Tstat.fused_STD_mm;
Tstat.overlapOnly_RMSE_um  = errScale * Tstat.overlapOnly_RMSE_mm;
Tstat.inner_RMSE_um        = errScale * Tstat.inner_RMSE_mm;
Tstat.robustPV995_um       = errScale * Tstat.robustPV995_mm;

summary = Tstat; %#ok<NASGU>
save(fullfile(outDir, 'step3_gaussian_all_cases.mat'), ...
    'allStep3', 'Tstat', 'summary', 'params', 'cfgEval', 'errScale');

if params.saveSummaryCsv
    writetable(Tstat, fullfile(outDir, 'step3_gaussian_summary.csv'));
end

fprintf('\n===================== 结果汇总 =====================\n');
disp(Tstat(:, {'Case', 'e_t_mm', 'e_R_deg', 'overlap_RMSE_um', ...
    'fused_RMSE_um', 'overlapOnly_RMSE_um', 'inner_RMSE_um', 'robustPV995_um'}));
fprintf('结果已保存到: %s\n', outDir);

%% ========================================================================
%%                               局部函数
%% ========================================================================

function cases = filter_cases(cases, caseMode)
    if strcmpi(caseMode, 'all')
        return;
    end

    keep = false(1, numel(cases));
    for i = 1:numel(cases)
        nm = cases{i}.caseName;
        keep(i) = strcmpi(nm, caseMode) || startsWith(nm, [caseMode '_'], 'IgnoreCase', true);
    end

    if ~any(keep)
        error('未找到指定 case: %s', caseMode);
    end

    cases = cases(keep);
end

function name = get_case_name_from_step2(R2)
    if isfield(R2, 'caseName') && ~isempty(R2.caseName)
        name = char(string(R2.caseName));
    else
        name = 'unknown_case';
    end
end

function ttl = get_case_title_from_step2(R2)
    if isfield(R2, 'caseTitle') && ~isempty(R2.caseTitle)
        ttl = char(string(R2.caseTitle));
    elseif isfield(R2, 'caseName') && ~isempty(R2.caseName)
        ttl = char(string(R2.caseName));
    else
        ttl = 'unknown title';
    end
end

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

    T.mask = T.mask & isfinite(T.Z);
end

function S = normalize_subap_from_step2(Sin)
    S = Sin;

    if isfield(S, 'c') && iscolumn(S.c), S.c = S.c.'; end
    if isfield(S, 'c_nom') && iscolumn(S.c_nom), S.c_nom = S.c_nom.'; end
    if isfield(S, 'err6') && iscolumn(S.err6), S.err6 = S.err6.'; end

    if isfield(S, 'X_local')
        S.X_local = double(S.X_local);
    else
        error('subap 缺少 X_local 字段。');
    end
    if isfield(S, 'Y_local')
        S.Y_local = double(S.Y_local);
    else
        error('subap 缺少 Y_local 字段。');
    end
    if isfield(S, 'Z_meas')
        S.Z_meas = double(S.Z_meas);
    else
        error('subap 缺少 Z_meas 字段。');
    end

    if isfield(S, 'mask')
        S.mask = logical(S.mask);
    else
        S.mask = isfinite(S.Z_meas);
    end

    if isfield(S, 'mask_local')
        S.mask_local = logical(S.mask_local);
    else
        if isfield(S, 'Rsub')
            S.mask_local = (S.X_local.^2 + S.Y_local.^2 <= S.Rsub^2);
        else
            S.mask_local = isfinite(S.Z_meas);
        end
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

    if ~isfield(S, 'R_w_nom')
        S.R_w_nom = S.R_w;
    else
        S.R_w_nom = double(S.R_w_nom);
    end

    if ~isfield(S, 'T_w_nom')
        S.T_w_nom = S.T_w;
    else
        S.T_w_nom = double(S.T_w_nom(:));
    end
end

function M = build_measured_subap_from_step2(Sap, local, maxSearchPts, maxEvalPts)
    X = double(local.X);
    Y = double(local.Y);
    Z = double(local.Z);
    mask = logical(local.mask) & isfinite(Z);

    valid = mask;
    xL = X(valid);
    yL = Y(valid);
    zL = Z(valid);

    if isempty(xL)
        error('子孔径有效点为空。');
    end

    P = Sap.R_w * [xL(:).'; yL(:).'; zL(:).'] + Sap.T_w(:);
    xW = P(1,:).';
    yW = P(2,:).';
    zW = P(3,:).';

    xVec = double(X(1,:));
    yVec = double(Y(:,1));
    Zg   = double(Z);
    Zg(~mask) = NaN;

    if numel(yVec) > 1 && yVec(2) < yVec(1)
        yVec = flipud(yVec);
        Zg   = flipud(Zg);
    end
    if numel(xVec) > 1 && xVec(2) < xVec(1)
        xVec = fliplr(xVec);
        Zg   = fliplr(Zg);
    end

    assert(all(diff(yVec) > 0), 'yVec 不是严格单调递增。');
    assert(all(diff(xVec) > 0), 'xVec 不是严格单调递增。');

    Fgrid = griddedInterpolant({yVec, xVec}, Zg, 'linear', 'none');

    idxSearch = pick_subsample_indices(numel(xL), maxSearchPts);
    idxEval   = pick_subsample_indices(numel(xL), maxEvalPts);

    boundaryLocal = extract_mask_boundary(X, Y, mask);
    boundaryWorld = local_boundary_to_world(boundaryLocal, Sap.R_w, Sap.T_w);

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
    M.boundaryLocal = boundaryLocal;
    M.boundaryWorld = boundaryWorld;
end

function idx = pick_subsample_indices(n, maxN)
    if isinf(maxN) || maxN >= n || maxN <= 0
        idx = (1:n).';
        return;
    end
    idx = round(linspace(1, n, maxN));
    idx = unique(max(min(idx, n), 1)).';
end

function modelA = build_surface_model_A(R2, params)
    XA = double(R2.localA.X);
    YA = double(R2.localA.Y);
    ZA = double(R2.localA.Z);
    MA = logical(R2.localA.mask) & isfinite(ZA);

    if isfield(R2, 'sagA_smooth') && ~isempty(R2.sagA_smooth)
        Zsmooth = double(R2.sagA_smooth);
    else
        Zsmooth = masked_gaussian_filter(ZA, MA, params.normalSmoothSigmaPix);
    end

    Zfill = fill_invalid_by_nearest(XA, YA, Zsmooth, MA);
    ds = infer_local_ds(XA, YA);

    [dZdx, dZdy] = gradient(Zfill, ds, ds);

    modelA = struct();
    modelA.X     = XA;
    modelA.Y     = YA;
    modelA.xVec  = XA(1,:);
    modelA.yVec  = YA(:,1);
    modelA.Zfill = Zfill;
    modelA.dZdx  = dZdx;
    modelA.dZdy  = dZdy;
    modelA.mask  = MA;
end

function ds = infer_local_ds(X, Y)
    dx = [];
    dy = [];
    if size(X,2) >= 2
        dx = median(abs(diff(X(1,:))), 'omitnan');
    end
    if size(Y,1) >= 2
        dy = median(abs(diff(Y(:,1))), 'omitnan');
    end
    cand = [dx, dy];
    cand = cand(isfinite(cand) & cand > 0);
    if isempty(cand)
        error('无法自动推断局部网格间隔 ds。');
    end
    ds = cand(1);
end

function p0 = build_initial_pose_from_step2(fine)
    p0 = zeros(6,1);
    p0(1) = 0;
    p0(2) = 0;
    p0(3) = fine.thetaRad;
    p0(4) = fine.tx;
    p0(5) = fine.ty;
    p0(6) = 0;
end

function [PBcand, meta] = build_candidate_points_from_step2(R2, p0, params)
    XA = double(R2.localA.X);
    YA = double(R2.localA.Y);
    ZA = double(R2.localA.Z);
    MA = logical(R2.localA.mask) & isfinite(ZA);

    XB = double(R2.localB.X);
    YB = double(R2.localB.Y);
    ZB = double(R2.localB.Z);
    MB = logical(R2.localB.mask) & isfinite(ZB);

    if params.useStep2OverlapOnly
        if isfield(R2, 'fine') && isfield(R2.fine, 'overlapMask') && ~isempty(R2.fine.overlapMask)
            overlapMaskA = logical(R2.fine.overlapMask) & MA;
        elseif isfield(R2, 'overlapMaskEst') && ~isempty(R2.overlapMaskEst)
            overlapMaskA = logical(R2.overlapMaskEst) & MA;
        else
            overlapMaskA = MA;
        end
    else
        overlapMaskA = MA;
    end

    xA = XA(overlapMaskA);
    yA = YA(overlapMaskA);

    [R0, t0] = posevec_to_rt_internal(p0);
    Rxy = R0(1:2,1:2);
    txy = t0(1:2);

    qA = [xA.'; yA.'];
    qB = Rxy.' * (qA - txy);

    xB = qB(1,:).';
    yB = qB(2,:).';

    zB = interp2(XB(1,:), YB(:,1), ZB, xB, yB, 'linear', NaN);
    mB = interp2(XB(1,:), YB(:,1), double(MB), xB, yB, 'linear', 0) >= 0.95;

    valid = isfinite(zB) & mB;
    PBcand = [xB(valid), yB(valid), zB(valid)];

    meta = struct();
    meta.numRawAOverlap = nnz(overlapMaskA);
    meta.numCandidateFromB = size(PBcand,1);
end

function [lb, ub] = build_pose_bounds(p0, params)
    lb = p0;
    ub = p0;

    lb(1) = -deg2rad(params.bound.rxDeg);
    ub(1) =  deg2rad(params.bound.rxDeg);

    lb(2) = -deg2rad(params.bound.ryDeg);
    ub(2) =  deg2rad(params.bound.ryDeg);

    lb(3) = p0(3) - deg2rad(params.bound.rzDeg);
    ub(3) = p0(3) + deg2rad(params.bound.rzDeg);

    lb(4) = p0(4) - params.bound.txMM;
    ub(4) = p0(4) + params.bound.txMM;

    lb(5) = p0(5) - params.bound.tyMM;
    ub(5) = p0(5) + params.bound.tyMM;

    lb(6) = -params.bound.tzMM;
    ub(6) =  params.bound.tzMM;
end

function r = residual_pose_point2surface(p, PB, modelA, params)
    [R, t] = posevec_to_rt_internal(p);
    PA = apply_transform(PB, R, t);

    x = PA(:,1);
    y = PA(:,2);
    z = PA(:,3);

    zA = interp2(modelA.xVec, modelA.yVec, modelA.Zfill, x, y, 'linear', NaN);
    gx = interp2(modelA.xVec, modelA.yVec, modelA.dZdx,  x, y, 'linear', NaN);
    gy = interp2(modelA.xVec, modelA.yVec, modelA.dZdy,  x, y, 'linear', NaN);
    mv = interp2(modelA.xVec, modelA.yVec, double(modelA.mask), x, y, 'linear', 0);

    valid = isfinite(zA) & isfinite(gx) & isfinite(gy) & (mv > 0);

    N = size(PB,1);
    rPlane = zeros(N,1);

    if any(valid)
        denom = sqrt(1 + gx(valid).^2 + gy(valid).^2);
        rPlane(valid) = (z(valid) - zA(valid)) ./ denom;
    end

    mv = min(max(mv, 0), 1);
    rMask = params.lambdaMaskPenalty * (1 - mv);

    rReg = [ ...
        params.lambdaRxRy * p(1); ...
        params.lambdaRxRy * p(2); ...
        params.lambdaTz   * p(6)];

    r = [rPlane; rMask; rReg];
end

function stat = evaluate_pose_against_surface(PB, p, modelA, params)
    [R, t] = posevec_to_rt_internal(p);
    PA = apply_transform(PB, R, t);

    x = PA(:,1);
    y = PA(:,2);
    z = PA(:,3);

    zA = interp2(modelA.xVec, modelA.yVec, modelA.Zfill, x, y, 'linear', NaN);
    gx = interp2(modelA.xVec, modelA.yVec, modelA.dZdx,  x, y, 'linear', NaN);
    gy = interp2(modelA.xVec, modelA.yVec, modelA.dZdy,  x, y, 'linear', NaN);
    mv = interp2(modelA.xVec, modelA.yVec, double(modelA.mask), x, y, 'linear', 0);

    valid = isfinite(zA) & isfinite(gx) & isfinite(gy) & (mv >= params.maskInterpThresh);

    dz = z(valid) - zA(valid);
    denom = sqrt(1 + gx(valid).^2 + gy(valid).^2);
    rPlane = dz ./ denom;

    stat = struct();
    stat.numValid = nnz(valid);

    if isempty(dz)
        stat.rmsePlane = NaN;
        stat.maePlane  = NaN;
        stat.rmseZ     = NaN;
        stat.maeZ      = NaN;
        stat.pvZ       = NaN;
        stat.stdZ      = NaN;
    else
        stat.rmsePlane = sqrt(mean(rPlane.^2));
        stat.maePlane  = mean(abs(rPlane));
        stat.rmseZ     = sqrt(mean(dz.^2));
        stat.maeZ      = mean(abs(dz));
        stat.pvZ       = max(dz) - min(dz);
        stat.stdZ      = std(dz, 0);
    end
end

function reg = gaussian_pose_to_reg(p)
    [R, t] = posevec_to_rt_internal(p);

    reg = struct();
    reg.p = [t(1), t(2), t(3), p(1), p(2), p(3)];
    reg.R = R;
    reg.score = NaN;
    reg.overlapRMSE = NaN;
    reg.overlapRatio = NaN;
    reg.overlapN = NaN;
    reg.residual = [];
    reg.xAov = [];
    reg.yAov = [];
    reg.zAref = [];
    reg.zAfromB = [];
end

function reg = attach_overlap_metrics_to_reg(reg, A, B, cfg)
    idx = B.idxEval;

    tx = reg.p(1); ty = reg.p(2); tz = reg.p(3);
    R  = reg.R;

    xB = B.xLocal(idx);
    yB = B.yLocal(idx);
    zB = B.zLocal(idx);

    PB = [xB(:).'; yB(:).'; zB(:).'];
    PA = R * PB + [tx; ty; tz];

    xAq = PA(1,:).';
    yAq = PA(2,:).';
    zAfromB = PA(3,:).';

    zAref = BypassFgridEval(A.Fgrid, yAq, xAq);
    valid = isfinite(zAref) & isfinite(zAfromB);

    overlapN = nnz(valid);
    overlapBase = min(A.n, B.n);
    overlapRatio = overlapN / max(overlapBase, 1);

    if overlapN < max(cfg.search.minOverlapPts, 1) || overlapRatio < cfg.search.minOverlapFrac
        overlapRMSE = inf;
        residual = [];
        xAov = [];
        yAov = [];
        zArefv = [];
        zAfromBv = [];
    else
        residual = zAfromB(valid) - zAref(valid);
        overlapRMSE = safe_rmse(residual);
        xAov = xAq(valid);
        yAov = yAq(valid);
        zArefv = zAref(valid);
        zAfromBv = zAfromB(valid);
    end

    reg.overlapRMSE  = overlapRMSE;
    reg.overlapRatio = overlapRatio;
    reg.overlapN     = overlapN;
    reg.residual     = residual;
    reg.xAov         = xAov;
    reg.yAov         = yAov;
    reg.zAref        = zArefv;
    reg.zAfromB      = zAfromBv;
end

function v = BypassFgridEval(Fgrid, yq, xq)
    v = Fgrid(yq, xq);
end

function truthRel = true_relative_pose(SA, SB)
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

function poseErr = evaluate_pose_error(reg, truthRel)
    dt = reg.p(1:3) - truthRel.t(:).';
    deul = [ ...
        wrap_to_pi_local(reg.p(4) - truthRel.rx), ...
        wrap_to_pi_local(reg.p(5) - truthRel.ry), ...
        wrap_to_pi_local(reg.p(6) - truthRel.rz)];

    poseErr = struct();
    poseErr.dt_mm    = dt(:).';
    poseErr.deul_rad = deul(:).';
    poseErr.et_mm    = norm(dt);
    poseErr.eR_deg   = rotm_geodesic_deg(reg.R, truthRel.R);
end

function deg = rotm_geodesic_deg(R1, R2)
    Rerr = R1 * R2.';
    val = (trace(Rerr) - 1) / 2;
    val = min(max(val, -1), 1);
    deg = rad2deg(acos(val));
end

function [fusedMap, errMap, unionMask, stats, boundaryA, boundaryBestB, ...
    ZAmap, ZBmap, overlapMaskMap, innerMask] = ...
    fuse_and_evaluate_6dof(T, A, B, reg, cfg)

    xAw = A.xWorld;
    yAw = A.yWorld;
    zAw = A.zWorld;

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

    fusedMap = nan(size(T.Z));
    fusedMap(overlapMask) = 0.5 * (ZAmap(overlapMask) + ZBmap(overlapMask));
    fusedMap(onlyA)       = ZAmap(onlyA);
    fusedMap(onlyB)       = ZBmap(onlyB);

    errMap = nan(size(T.Z));
    errMap(unionMask) = fusedMap(unionMask) - T.Z(unionMask);

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

    stats.robustPV995 = safe_robust_pv995(ev);

    overlapMaskMap = overlapMask;

    boundaryA = A.boundaryWorld;
    boundaryBestB = boundary_local_to_world_via_reg(B.boundaryLocal, A.pose, reg);
end

function boundaryWorld = boundary_local_to_world_via_reg(boundaryLocalB, poseA, reg)
    if isempty(boundaryLocalB)
        boundaryWorld = zeros(0,2);
        return;
    end

    p = reg.p(:);
    R = reg.R;

    PB = [boundaryLocalB(:,1).'; boundaryLocalB(:,2).'; zeros(1, size(boundaryLocalB,1))];
    PA = R * PB + p(1:3);

    PW = poseA.R_w * PA + poseA.T_w(:);
    boundaryWorld = [PW(1,:).', PW(2,:).'];
end

function boundaryWorld = local_boundary_to_world(boundaryLocal, Rw, Tw)
    if isempty(boundaryLocal)
        boundaryWorld = zeros(0,2);
        return;
    end

    P = Rw * [boundaryLocal(:,1).'; boundaryLocal(:,2).'; zeros(1, size(boundaryLocal,1))] + Tw(:);
    boundaryWorld = [P(1,:).', P(2,:).'];
end

function b = extract_mask_boundary(X, Y, mask)
    mask = logical(mask);
    if ~any(mask(:))
        b = zeros(0,2);
        return;
    end

    xVec = X(1,:);
    yVec = Y(:,1);

    C = contourc(xVec, yVec, double(mask), [0.5 0.5]);
    if isempty(C)
        b = zeros(0,2);
        return;
    end

    segs = parse_contourc_segments(C);
    if isempty(segs)
        b = zeros(0,2);
        return;
    end

    segLen = zeros(numel(segs), 1);
    for i = 1:numel(segs)
        xy = segs{i};
        if size(xy,1) < 2
            segLen(i) = 0;
        else
            dxy = diff(xy, 1, 1);
            segLen(i) = sum(sqrt(sum(dxy.^2, 2)));
        end
    end

    [~, idx] = max(segLen);
    b = segs{idx};

    if ~isempty(b) && norm(b(1,:) - b(end,:)) > 0
        b(end+1,:) = b(1,:);
    end
end

function segs = parse_contourc_segments(C)
    segs = {};
    k = 1;
    while k < size(C,2)
        npt = C(2,k);
        if npt <= 0 || (k + npt) > size(C,2)
            break;
        end
        xy = C(:, k+1:k+npt).';
        segs{end+1} = xy; %#ok<AGROW>
        k = k + npt + 1;
    end
end

function m = erode_mask_disk(mask, radiusPix)
    mask = logical(mask);
    if nargin < 2 || radiusPix <= 0
        m = mask;
        return;
    end

    rad = max(1, round(radiusPix));
    [Xk, Yk] = meshgrid(-rad:rad, -rad:rad);
    K = double((Xk.^2 + Yk.^2) <= rad^2);

    cnt = conv2(double(mask), K, 'same');
    m = (cnt >= sum(K(:)) - 1e-12) & mask;
end

function v = safe_rmse(x)
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = sqrt(mean(x.^2));
    end
end

function v = safe_mae(x)
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = mean(abs(x));
    end
end

function v = safe_pv(x)
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = max(x) - min(x);
    end
end

function v = safe_std(x)
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        v = std(x, 1);
    end
end

function v = safe_robust_pv995(x)
    x = x(isfinite(x));
    if isempty(x)
        v = NaN;
    else
        q = prctile(x, [0.25, 99.75]);
        v = q(2) - q(1);
    end
end

function [R, t] = posevec_to_rt_internal(p)
    rx = p(1);
    ry = p(2);
    rz = p(3);

    Rx = [1, 0, 0; ...
          0, cos(rx), -sin(rx); ...
          0, sin(rx),  cos(rx)];

    Ry = [ cos(ry), 0, sin(ry); ...
           0,       1, 0; ...
          -sin(ry), 0, cos(ry)];

    Rz = [cos(rz), -sin(rz), 0; ...
          sin(rz),  cos(rz), 0; ...
          0,        0,       1];

    R = Rz * Ry * Rx;
    t = p(4:6);
end

function P2 = apply_transform(P1, R, t)
    P2 = (R * P1.' + t);
    P2 = P2.';
end

function R = eul_zyx_to_rotm(rx, ry, rz)
    Rx = [1, 0, 0; ...
          0, cos(rx), -sin(rx); ...
          0, sin(rx),  cos(rx)];

    Ry = [ cos(ry), 0, sin(ry); ...
           0,       1, 0; ...
          -sin(ry), 0, cos(ry)];

    Rz = [cos(rz), -sin(rz), 0; ...
          sin(rz),  cos(rz), 0; ...
          0,        0,       1];

    R = Rz * Ry * Rx;
end

function [rx, ry, rz] = rotm_to_eul_zyx(R)
    sy = -R(3,1);
    sy = min(max(sy, -1), 1);
    ry = asin(sy);

    cy = cos(ry);
    if abs(cy) > 1e-12
        rx = atan2(R(3,2), R(3,3));
        rz = atan2(R(2,1), R(1,1));
    else
        rx = atan2(-R(2,3), R(2,2));
        rz = 0;
    end
end

function Zs = masked_gaussian_filter(Z, mask, sigmaPix)
    if nargin < 3 || sigmaPix <= 0
        Zs = Z;
        return;
    end

    mask = logical(mask) & isfinite(Z);
    Z0 = Z;
    Z0(~mask) = 0;

    G = gaussian_kernel_2d(sigmaPix);
    num = conv2(Z0, G, 'same');
    den = conv2(double(mask), G, 'same');

    Zs = num ./ max(den, eps);
    Zs(~mask) = NaN;
end

function G = gaussian_kernel_2d(sigmaPix)
    rad = max(1, ceil(3 * sigmaPix));
    [X, Y] = meshgrid(-rad:rad, -rad:rad);
    G = exp(-(X.^2 + Y.^2) / (2 * sigmaPix^2));
    G = G / sum(G(:));
end

function Zfill = fill_invalid_by_nearest(X, Y, Z, mask)
    valid = mask & isfinite(Z);

    if ~any(valid(:))
        Zfill = zeros(size(Z));
        return;
    end

    xv = X(valid);
    yv = Y(valid);
    zv = Z(valid);

    F = scatteredInterpolant(xv(:), yv(:), zv(:), 'nearest', 'nearest');
    Zfill = F(X, Y);
end

function ang = wrap_to_pi_local(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end