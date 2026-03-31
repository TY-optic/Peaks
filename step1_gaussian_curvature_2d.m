%% step1_gaussian_curvature_2d_6dof.m
% 基于 Gaussian curvature 的二维粗配准
% 仅支持读取 figures\overlap_cases_6dof 下的三个独立 case MAT 文件
%
% 输入：
%   figures\overlap_cases_6dof\case_*_6dof.mat
%
% 约定（与 build_overlap_cases_6dof_and_save_three_mats.m 保持一致）：
%   每个 MAT 文件中至少包含：
%       truth, subap_A, subap_B, caseName, caseTitle, etaTarget, etaNominal
%   其中 subap_A / subap_B 主要字段：
%       X_local, Y_local, Z_meas, mask, mask_local,
%       c_meas_xy, theta_meas_xy, c, theta, Rsub, ds
%
% 输出：
%   results\step2_gaussian_curvature_2d_6dof\step2_result_<case>.mat
%   results\step2_gaussian_curvature_2d_6dof\step2_gaussian_all_cases.mat
%   results\step2_gaussian_curvature_2d_6dof\step2_gaussian_summary.csv
%
% 说明：
%   1) 不再回退到单面型模式；
%   2) 不生成任何绘图；
%   3) 直接使用各 case 中的 Z_meas 数据做二维高斯曲率配准；
%   4) 真值采用 measured xy 投影：c_meas_xy / theta_meas_xy；
%   5) 仅输出二维条件估计值与相对真值误差。

clear; clc;

%% ========================= 路径与参数 ===================================
inDir = fullfile('figures', 'overlap_cases_6dof');
outDir = fullfile('results', 'step2_gaussian_curvature_2d_6dof');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

params = struct();

% 预处理
params.sagSmoothSigmaPix = 1.25;   % sag 平滑
params.curvMarginPix     = 3;      % 曲率评价内缩
params.clipSigma         = 6.0;    % 曲率标准化裁剪

% 角度粗搜索
params.thetaSearchDeg = -4.0 : 0.05 : 4.0;

% 粗搜索最小重叠像素
params.minOverlapPix = 1200;

% 细化搜索
params.doFineRefine     = true;
params.thetaFineSpanDeg = 0.12;
params.thetaFineNum     = 13;
params.shiftFineSpanPx  = 1.2;
params.shiftFineNum     = 13;
params.maskInterpThresh = 0.85;

% 输出控制
params.savePerCaseMat = true;
params.saveSummaryCsv = true;

%% ========================= 读取 case 文件 ================================
matFiles = dir(fullfile(inDir, 'case_*_6dof.mat'));
if isempty(matFiles)
    error('未找到输入文件：%s', fullfile(inDir, 'case_*_6dof.mat'));
end

matFiles = sort_case_files(matFiles);

fprintf('Input directory: %s\n', inDir);
fprintf('Detected %d case file(s).\n\n', numel(matFiles));

allResults = cell(1, numel(matFiles));

%% ========================= 主循环 =======================================
for kCase = 1:numel(matFiles)
    matPath = fullfile(matFiles(kCase).folder, matFiles(kCase).name);
    S = load(matPath);

    C = load_case_from_file(S, matPath);

    fprintf('============================================================\n');
    fprintf('Running case: %s\n', C.name);
    fprintf('File       : %s\n', matFiles(kCase).name);
    fprintf('============================================================\n');

    subA = normalize_subap_meas_struct(C.subapA);
    subB = normalize_subap_meas_struct(C.subapB);

    % -------- 检查网格一致性 --------
    assert(isequal(size(subA.X), size(subA.Y), size(subA.Z), size(subA.mask)), ...
        'A 子孔径网格尺寸不一致。');
    assert(isequal(size(subB.X), size(subB.Y), size(subB.Z), size(subB.mask)), ...
        'B 子孔径网格尺寸不一致。');

    dsA = subA.ds;
    dsB = subB.ds;
    if abs(dsA - dsB) > 1e-12
        warning('A/B 网格间隔不完全一致：dsA=%.12g, dsB=%.12g，后续按各自网格处理。', dsA, dsB);
    end

    % -------- 构造 Gaussian curvature 图 --------
    [KA, maskKA, ZA_smooth, auxA] = build_gaussian_curvature_map(subA, params);
    [KB, maskKB, ZB_smooth, auxB] = build_gaussian_curvature_map(subB, params);

    % -------- 真值二维变换（B -> A 的局部坐标系）--------
    gt2d = compute_ground_truth_transform_2d(subA, subB);

    % -------- 粗搜索：旋转 + 整数像素平移 --------
    coarse = coarse_search_gaussian_se2( ...
        KA, maskKA, KB, maskKB, ...
        subA.X, subA.Y, subB.X, subB.Y, params);

    % -------- 细化：亚像素 / 亚角度 --------
    if params.doFineRefine
        fine = refine_search_gaussian_se2( ...
            KA, maskKA, KB, maskKB, ...
            subA.X, subA.Y, subB.X, subB.Y, ...
            coarse, params);
    else
        fine = coarse_to_fine_passthrough(coarse);
    end

    % -------- 以最佳变换将 B 重采样到 A 坐标系 --------
    [KB_reg, maskKB_reg] = warp_B_to_A_grid( ...
        KB, maskKB, subB.X, subB.Y, subA.X, subA.Y, ...
        fine.thetaRad, fine.tx, fine.ty);

    overlapMaskEst = maskKA & maskKB_reg & isfinite(KA) & isfinite(KB_reg);

    % -------- 误差统计 --------
    poseErr = struct();
    poseErr.tx        = fine.tx - gt2d.tx;
    poseErr.ty        = fine.ty - gt2d.ty;
    poseErr.thetaDeg  = rad2deg(wrap_to_pi_local(fine.thetaRad - gt2d.thetaRad));
    poseErr.transNorm = hypot(poseErr.tx, poseErr.ty);

    caseResult = struct();
    caseResult.caseName = C.name;
    caseResult.caseTitle = C.title;
    caseResult.sourceMatFile = matPath;
    caseResult.etaTarget  = C.etaTarget;
    caseResult.etaNominal = C.etaNominal;

    caseResult.truth = C.truth;
    caseResult.subapA = C.subapA;
    caseResult.subapB = C.subapB;

    caseResult.localA = subA;
    caseResult.localB = subB;

    caseResult.sagA_smooth = ZA_smooth;
    caseResult.sagB_smooth = ZB_smooth;

    caseResult.curvatureA = KA;
    caseResult.curvatureB = KB;
    caseResult.curvatureA_mask = maskKA;
    caseResult.curvatureB_mask = maskKB;

    caseResult.coarse = coarse;
    caseResult.fine   = fine;

    caseResult.truth2d   = gt2d;
    caseResult.poseError = poseErr;

    caseResult.registeredB.curvature = KB_reg;
    caseResult.registeredB.mask      = maskKB_reg;
    caseResult.overlapMaskEst        = overlapMaskEst;

    caseResult.auxA = auxA;
    caseResult.auxB = auxB;

    if isfield(C, 'validation')
        caseResult.validation = C.validation;
    end
    if isfield(C, 'meta')
        caseResult.meta = C.meta;
    end

    % -------- 保存单案例结果 --------
    if params.savePerCaseMat
        save(fullfile(outDir, sprintf('step2_result_%s.mat', C.name)), 'caseResult');
    end

    % -------- 打印 --------
    fprintf('Truth2D  : tx = %+9.4f mm, ty = %+9.4f mm, theta = %+8.4f deg\n', ...
        gt2d.tx, gt2d.ty, rad2deg(gt2d.thetaRad));
    fprintf('Estimate : tx = %+9.4f mm, ty = %+9.4f mm, theta = %+8.4f deg\n', ...
        fine.tx, fine.ty, rad2deg(fine.thetaRad));
    fprintf('Error    : dtx= %+9.4f mm, dty= %+9.4f mm, dth   = %+8.4f deg\n', ...
        poseErr.tx, poseErr.ty, poseErr.thetaDeg);
    fprintf('TransErr : %.6f mm\n', poseErr.transNorm);
    fprintf('BestNCC  : %.6f\n', fine.score);
    fprintf('Overlap  : %d pixels\n\n', nnz(overlapMaskEst));

    allResults{kCase} = caseResult;
end

%% ========================= 汇总保存 =====================================
summary = build_summary_table(allResults);

save(fullfile(outDir, 'step2_gaussian_all_cases.mat'), ...
    'allResults', 'summary', 'params');

if params.saveSummaryCsv
    writetable(summary, fullfile(outDir, 'step2_gaussian_summary.csv'));
end

disp(summary);
fprintf('\n已保存总结果: %s\n', fullfile(outDir, 'step2_gaussian_all_cases.mat'));
if params.saveSummaryCsv
    fprintf('已保存汇总表: %s\n', fullfile(outDir, 'step2_gaussian_summary.csv'));
end

%% ========================================================================
%% 局部函数
%% ========================================================================

function matFiles = sort_case_files(matFiles)
    if isempty(matFiles)
        return;
    end

    names = {matFiles.name};
    key = inf(size(names));

    for i = 1:numel(names)
        tok = regexp(names{i}, 'case_(O\d+)_', 'tokens', 'once');
        if ~isempty(tok)
            key(i) = parse_overlap_order(tok{1});
        end
    end

    [~, idx] = sortrows([key(:), (1:numel(names)).']);
    matFiles = matFiles(idx);
end

function k = parse_overlap_order(tag)
    tag = upper(string(tag));
    switch tag
        case "O1"
            k = 1;
        case "O2"
            k = 2;
        case "O3"
            k = 3;
        otherwise
            k = inf;
    end
end

function C = load_case_from_file(S, matPath)
    reqFields = {'truth', 'subap_A', 'subap_B'};
    for i = 1:numel(reqFields)
        if ~isfield(S, reqFields{i})
            error('文件缺少字段 %s ：%s', reqFields{i}, matPath);
        end
    end

    C = struct();
    C.truth  = S.truth;
    C.subapA = S.subap_A;
    C.subapB = S.subap_B;

    if isfield(S, 'caseName')
        C.name = char(string(S.caseName));
    else
        [~, nm] = fileparts(matPath);
        C.name = nm;
    end

    if isfield(S, 'caseTitle')
        C.title = char(string(S.caseTitle));
    else
        C.title = C.name;
    end

    if isfield(S, 'etaTarget')
        C.etaTarget = double(S.etaTarget);
    else
        C.etaTarget = NaN;
    end

    if isfield(S, 'etaNominal')
        C.etaNominal = double(S.etaNominal);
    else
        C.etaNominal = NaN;
    end

    if isfield(S, 'validation')
        C.validation = S.validation;
    end

    if isfield(S, 'meta')
        C.meta = S.meta;
    end
end

function sub = normalize_subap_meas_struct(Sin)
    sub = struct();

    % 局部坐标
    if isfield(Sin, 'X_local')
        sub.X = double(Sin.X_local);
    elseif isfield(Sin, 'X')
        sub.X = double(Sin.X);
    else
        error('subap 缺少 X_local/X 字段。');
    end

    if isfield(Sin, 'Y_local')
        sub.Y = double(Sin.Y_local);
    elseif isfield(Sin, 'Y')
        sub.Y = double(Sin.Y);
    else
        error('subap 缺少 Y_local/Y 字段。');
    end

    % 测量 sag
    if isfield(Sin, 'Z_meas')
        sub.Z = double(Sin.Z_meas);
    elseif isfield(Sin, 'Z')
        sub.Z = double(Sin.Z);
    else
        error('subap 缺少 Z_meas/Z 字段。');
    end

    % mask
    if isfield(Sin, 'mask')
        sub.mask = logical(Sin.mask);
    else
        sub.mask = isfinite(sub.Z);
    end

    if isfield(Sin, 'mask_local')
        sub.maskLocal = logical(Sin.mask_local);
        sub.mask = sub.mask & sub.maskLocal;
    else
        sub.maskLocal = true(size(sub.mask));
    end

    sub.mask = sub.mask & isfinite(sub.Z);

    % 2D 真值优先采用 measured xy 投影
    if isfield(Sin, 'c_meas_xy')
        sub.c2d = double(Sin.c_meas_xy(:)).';
        sub.centerSource = 'c_meas_xy';
    elseif isfield(Sin, 'c')
        sub.c2d = double(Sin.c(:)).';
        sub.centerSource = 'c';
    else
        error('subap 缺少 c_meas_xy/c 字段。');
    end

    if isfield(Sin, 'theta_meas_xy')
        sub.theta2d = double(Sin.theta_meas_xy);
        sub.thetaSource = 'theta_meas_xy';
    elseif isfield(Sin, 'theta')
        sub.theta2d = double(Sin.theta);
        sub.thetaSource = 'theta';
    else
        sub.theta2d = 0;
        sub.thetaSource = 'default_zero';
    end

    if isfield(Sin, 'Rsub')
        sub.Rsub = double(Sin.Rsub);
    else
        sub.Rsub = NaN;
    end

    if isfield(Sin, 'ds')
        sub.ds = double(Sin.ds);
    else
        sub.ds = infer_grid_spacing_from_local_grid(sub.X, sub.Y);
    end

    sub.raw = Sin;
end

function ds = infer_grid_spacing_from_local_grid(X, Y)
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
        error('无法从局部网格自动推断 ds。');
    end

    ds = cand(1);
end

function [Knorm, evalMask, Zsmooth, aux] = build_gaussian_curvature_map(sub, params)
    Zsmooth = masked_gaussian_filter(sub.Z, sub.mask, params.sagSmoothSigmaPix);
    Zfill   = fill_invalid_by_nearest(sub.X, sub.Y, Zsmooth, sub.mask);

    dx = sub.ds;
    dy = sub.ds;

    [Zx, Zy]   = gradient(Zfill, dx, dy);
    [Zxx, Zxy] = gradient(Zx, dx, dy);
    [Zyx, Zyy] = gradient(Zy, dx, dy);
    Zxy = 0.5 * (Zxy + Zyx);

    K = (Zxx .* Zyy - Zxy.^2) ./ (1 + Zx.^2 + Zy.^2).^2;

    evalMask = erode_mask_square(sub.mask, params.curvMarginPix);
    evalMask = evalMask & isfinite(K);

    if nnz(evalMask) < 16
        warning('有效曲率像素过少，退化为原 mask。');
        evalMask = sub.mask & isfinite(K);
    end

    kval = K(evalMask);
    mu = median(kval, 'omitnan');
    sig = 1.4826 * median(abs(kval - mu), 'omitnan');

    if ~(isfinite(sig) && sig > 0)
        sig = std(kval, 0, 'omitnan');
    end
    if ~(isfinite(sig) && sig > 0)
        sig = 1;
    end

    Knorm = nan(size(K));
    Knorm(evalMask) = (K(evalMask) - mu) / sig;
    Knorm(evalMask) = max(min(Knorm(evalMask), params.clipSigma), -params.clipSigma);

    aux = struct();
    aux.Kraw = K;
    aux.Zx   = Zx;
    aux.Zy   = Zy;
    aux.Zxx  = Zxx;
    aux.Zxy  = Zxy;
    aux.Zyy  = Zyy;
end

function gt = compute_ground_truth_transform_2d(subA, subB)
    thA = subA.theta2d;
    thB = subB.theta2d;

    RA = [cos(thA), -sin(thA); ...
          sin(thA),  cos(thA)];

    dWorld = (subB.c2d(:) - subA.c2d(:));
    tLocalA = RA.' * dWorld;

    gt = struct();
    gt.tx = tLocalA(1);
    gt.ty = tLocalA(2);
    gt.thetaRad = wrap_to_pi_local(thB - thA);

    gt.centerA = subA.c2d;
    gt.centerB = subB.c2d;
    gt.thetaA  = thA;
    gt.thetaB  = thB;
    gt.centerSourceA = subA.centerSource;
    gt.centerSourceB = subB.centerSource;
    gt.thetaSourceA  = subA.thetaSource;
    gt.thetaSourceB  = subB.thetaSource;
end

function coarse = coarse_search_gaussian_se2(A, MA, B, MB, XA, YA, XB, YB, params)
    xVecA = XA(1,:);
    yVecA = YA(:,1);

    bestScore = -inf;
    bestThetaDeg = 0;
    bestRowShift = 0;
    bestColShift = 0;
    thetaScores = nan(size(params.thetaSearchDeg));

    for i = 1:numel(params.thetaSearchDeg)
        thDeg = params.thetaSearchDeg(i);
        thRad = deg2rad(thDeg);

        [Br, MBr] = warp_B_to_A_grid(B, MB, XB, YB, XA, YA, thRad, 0, 0);

        [scoreMap, overlapMap] = masked_ncc_full_integer(A, MA, Br, MBr, params.minOverlapPix);

        if all(~isfinite(scoreMap(:)))
            thetaScores(i) = NaN;
            continue;
        end

        [curBestScore, idx] = max(scoreMap(:));
        thetaScores(i) = curBestScore;

        if curBestScore > bestScore
            [ir, ic] = ind2sub(size(scoreMap), idx);

            m = size(A,1);
            n = size(A,2);

            rowShift = ir - m;
            colShift = ic - n;

            bestScore = curBestScore;
            bestThetaDeg = thDeg;
            bestRowShift = rowShift;
            bestColShift = colShift;
            bestOverlapMap = overlapMap; %#ok<NASGU>
        end
    end

    coarse = struct();
    coarse.score = bestScore;
    coarse.thetaDeg = bestThetaDeg;
    coarse.thetaRad = deg2rad(bestThetaDeg);
    coarse.rowShiftPx = bestRowShift;
    coarse.colShiftPx = bestColShift;
    coarse.ty = bestRowShift * median(diff(yVecA));
    coarse.tx = bestColShift * median(diff(xVecA));
    coarse.thetaScores = thetaScores;
    coarse.thetaSearchDeg = params.thetaSearchDeg;
end

function fine = refine_search_gaussian_se2(A, MA, B, MB, XA, YA, XB, YB, coarse, params)
    dsA = median(diff(XA(1,:)));

    theta0 = coarse.thetaRad;
    tx0 = coarse.tx;
    ty0 = coarse.ty;

    thetaVec = linspace(theta0 - deg2rad(params.thetaFineSpanDeg), ...
                        theta0 + deg2rad(params.thetaFineSpanDeg), ...
                        params.thetaFineNum);

    txVec = linspace(tx0 - params.shiftFineSpanPx*dsA, ...
                     tx0 + params.shiftFineSpanPx*dsA, ...
                     params.shiftFineNum);

    tyVec = linspace(ty0 - params.shiftFineSpanPx*dsA, ...
                     ty0 + params.shiftFineSpanPx*dsA, ...
                     params.shiftFineNum);

    bestScore = -inf;
    bestTheta = theta0;
    bestTx = tx0;
    bestTy = ty0;
    bestOverlapMask = false(size(A));

    for it = 1:numel(thetaVec)
        th = thetaVec(it);
        for ix = 1:numel(txVec)
            tx = txVec(ix);
            for iy = 1:numel(tyVec)
                ty = tyVec(iy);

                [score, overlapMask] = masked_ncc_direct( ...
                    A, MA, B, MB, ...
                    XA, YA, XB, YB, ...
                    th, tx, ty, ...
                    params.minOverlapPix, params.maskInterpThresh);

                if score > bestScore
                    bestScore = score;
                    bestTheta = th;
                    bestTx = tx;
                    bestTy = ty;
                    bestOverlapMask = overlapMask;
                end
            end
        end
    end

    fine = struct();
    fine.score = bestScore;
    fine.thetaRad = bestTheta;
    fine.thetaDeg = rad2deg(bestTheta);
    fine.tx = bestTx;
    fine.ty = bestTy;
    fine.overlapMask = bestOverlapMask;
end

function fine = coarse_to_fine_passthrough(coarse)
    fine = struct();
    fine.score = coarse.score;
    fine.thetaRad = coarse.thetaRad;
    fine.thetaDeg = coarse.thetaDeg;
    fine.tx = coarse.tx;
    fine.ty = coarse.ty;
    fine.overlapMask = [];
end

function [Bw, Mw] = warp_B_to_A_grid(B, MB, XB, YB, XA, YA, theta, tx, ty)
    xVecB = XB(1,:);
    yVecB = YB(:,1);

    c = cos(theta);
    s = sin(theta);

    Xb = c .* (XA - tx) + s .* (YA - ty);
    Yb = -s .* (XA - tx) + c .* (YA - ty);

    Bw = interp2(xVecB, yVecB, B,  Xb, Yb, 'linear', NaN);
    Mw = interp2(xVecB, yVecB, double(MB), Xb, Yb, 'linear', 0) >= 0.999;
end

function [scoreMap, overlapMap] = masked_ncc_full_integer(A, MA, B, MB, minOverlap)
    A0 = A;
    B0 = B;

    A0(~MA | ~isfinite(A0)) = 0;
    B0(~MB | ~isfinite(B0)) = 0;

    MA = double(MA);
    MB = double(MB);

    crossAB = conv2(A0, rot90(B0, 2), 'full');
    sumA    = conv2(A0, rot90(MB, 2), 'full');
    sumB    = conv2(MA, rot90(B0, 2), 'full');
    sumA2   = conv2(A0.^2, rot90(MB, 2), 'full');
    sumB2   = conv2(MA, rot90(B0.^2, 2), 'full');
    N       = conv2(MA, rot90(MB, 2), 'full');

    denomN = max(N, 1);
    num = crossAB - (sumA .* sumB) ./ denomN;
    denA = sumA2 - (sumA.^2) ./ denomN;
    denB = sumB2 - (sumB.^2) ./ denomN;

    den = sqrt(max(denA, 0) .* max(denB, 0));
    scoreMap = num ./ den;

    invalid = (N < minOverlap) | ~(den > 0);
    scoreMap(invalid) = -inf;

    overlapMap = N;
end

function [score, overlapMask] = masked_ncc_direct(A, MA, B, MB, XA, YA, XB, YB, theta, tx, ty, minOverlap, maskInterpThresh)
    xVecB = XB(1,:);
    yVecB = YB(:,1);

    c = cos(theta);
    s = sin(theta);

    Xb = c .* (XA - tx) + s .* (YA - ty);
    Yb = -s .* (XA - tx) + c .* (YA - ty);

    Bw = interp2(xVecB, yVecB, B, Xb, Yb, 'linear', NaN);
    Mw = interp2(xVecB, yVecB, double(MB), Xb, Yb, 'linear', 0) >= maskInterpThresh;

    overlapMask = MA & Mw & isfinite(A) & isfinite(Bw);
    if nnz(overlapMask) < minOverlap
        score = -inf;
        return;
    end

    a = A(overlapMask);
    b = Bw(overlapMask);

    a = a - mean(a, 'omitnan');
    b = b - mean(b, 'omitnan');

    da = sqrt(sum(a.^2));
    db = sqrt(sum(b.^2));

    if da <= 0 || db <= 0
        score = -inf;
    else
        score = sum(a .* b) / (da * db);
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

function maskErode = erode_mask_square(mask, marginPix)
    if nargin < 2 || marginPix <= 0
        maskErode = logical(mask);
        return;
    end

    ksz = 2 * marginPix + 1;
    K = ones(ksz, ksz);
    cnt = conv2(double(mask), K, 'same');
    maskErode = (cnt == numel(K));
end

function T = build_summary_table(allResults)
    n = numel(allResults);

    caseName   = strings(n,1);
    caseTitle  = strings(n,1);
    etaTarget  = nan(n,1);
    etaNominal = nan(n,1);

    txTrue = nan(n,1);
    tyTrue = nan(n,1);
    thTrue = nan(n,1);

    txEst = nan(n,1);
    tyEst = nan(n,1);
    thEst = nan(n,1);

    dtx    = nan(n,1);
    dty    = nan(n,1);
    dth    = nan(n,1);
    dtrans = nan(n,1);
    score  = nan(n,1);
    overlapPix = nan(n,1);

    centerSourceA = strings(n,1);
    centerSourceB = strings(n,1);
    thetaSourceA  = strings(n,1);
    thetaSourceB  = strings(n,1);

    for i = 1:n
        R = allResults{i};

        caseName(i)   = string(R.caseName);
        caseTitle(i)  = string(R.caseTitle);
        etaTarget(i)  = R.etaTarget;
        etaNominal(i) = R.etaNominal;

        txTrue(i) = R.truth2d.tx;
        tyTrue(i) = R.truth2d.ty;
        thTrue(i) = rad2deg(R.truth2d.thetaRad);

        txEst(i) = R.fine.tx;
        tyEst(i) = R.fine.ty;
        thEst(i) = R.fine.thetaDeg;

        dtx(i)    = R.poseError.tx;
        dty(i)    = R.poseError.ty;
        dth(i)    = R.poseError.thetaDeg;
        dtrans(i) = R.poseError.transNorm;
        score(i)  = R.fine.score;
        overlapPix(i) = nnz(R.overlapMaskEst);

        centerSourceA(i) = string(R.truth2d.centerSourceA);
        centerSourceB(i) = string(R.truth2d.centerSourceB);
        thetaSourceA(i)  = string(R.truth2d.thetaSourceA);
        thetaSourceB(i)  = string(R.truth2d.thetaSourceB);
    end

    T = table(caseName, caseTitle, etaTarget, etaNominal, ...
              txTrue, tyTrue, thTrue, ...
              txEst, tyEst, thEst, ...
              dtx, dty, dth, dtrans, score, overlapPix, ...
              centerSourceA, centerSourceB, thetaSourceA, thetaSourceB);
end

function ang = wrap_to_pi_local(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end