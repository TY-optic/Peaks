function [subA, subB, truth_data] = generate_noisy_peaks_subaps(sigma_n_um, Lc_mm, rng_seed)
% GENERATE_NOISY_PEAKS_SUBAPS 生成带有指定空间相关噪声的光学子孔径对
% 
% 输入参数:
%   sigma_n_um - 噪声的均方根(RMS)幅值，单位: 微米 (um)
%   Lc_mm      - 噪声的高斯相关长度，单位: 毫米 (mm)
%   rng_seed   - 随机数种子(可选)。用于保证蒙特卡洛实验的单次可复现性。
%
% 输出参数:
%   subA       - 子孔径 A 的数据结构体 (包含局部坐标 x, y, z 等)
%   subB       - 子孔径 B 的数据结构体
%   truth_data - 包含真实面形、噪声场和真实相对位姿的结构体，用于误差计算

    %% 1. 初始化与随机种子设置
    if nargin > 2 && ~isempty(rng_seed)
        rng(rng_seed);
    else
        rng('shuffle');
    end

    %% 2. 全局固定参数 (保持与原实验一致)
    ds = 0.5;                 % 网格分辨率 mm
    D_full = 100;             % 母面形直径 mm
    R_full = D_full / 2;
    R_sub  = 30;              % 子孔径半径 mm
    overlap_target = 0.50;    % 目标重叠比例 50%
    rot_deg = 2.0;            % B 相对于 A 的真实旋转角
    theta_true = deg2rad(rot_deg);
    phi_layout_deg = 25;      % 布局方向
    phi_layout = deg2rad(phi_layout_deg);
    target_pv_f1 = 2.5;       % 目标 PV 值 mm
    
    sigma_n = sigma_n_um * 1e-3; % 转换为 mm

    %% 3. 构造母面型
    xv = -R_full : ds : R_full;
    yv = -R_full : ds : R_full;
    [Xm, Ym] = meshgrid(xv, yv);
    Rm = hypot(Xm, Ym);
    mask_full = (Rm <= R_full);

    Xbar = 6 * Xm / D_full;
    Ybar = 6 * Ym / D_full;
    Z_base = peaks(Xbar, Ybar);
    Z_base(~mask_full) = NaN;

    % 去平面项并归一化
    coef0 = fit_plane_ls(Xm(mask_full), Ym(mask_full), Z_base(mask_full));
    Z_base = Z_base - (coef0(1)*Xm + coef0(2)*Ym + coef0(3));
    vmax = max(abs(Z_base(mask_full)));
    if vmax <= 0, vmax = 1; end
    Z_base(mask_full) = Z_base(mask_full) / vmax;

    % 缩放到目标 PV
    Z_clean = Z_base;
    coef1 = fit_plane_ls(Xm(mask_full), Ym(mask_full), Z_clean(mask_full));
    Z_clean(mask_full) = Z_clean(mask_full) - ...
        (coef1(1)*Xm(mask_full) + coef1(2)*Ym(mask_full) + coef1(3));
    
    pv0 = max(Z_clean(mask_full)) - min(Z_clean(mask_full));
    if pv0 > 0
        Z_clean(mask_full) = Z_clean(mask_full) * (target_pv_f1 / pv0);
    end

    %% 4. 生成特定频率与幅值的相关噪声
    noise_raw = randn(size(Xm));
    noise_raw(~mask_full) = 0;

    sigma_px = Lc_mm / ds;
    ker_half = max(3, ceil(4*sigma_px));
    [xk, yk] = meshgrid(-ker_half:ker_half, -ker_half:ker_half);
    gk = exp(-(xk.^2 + yk.^2) / (2*sigma_px^2));
    gk = gk / sum(gk(:));

    noise_corr = conv2(noise_raw, gk, 'same');
    noise_corr(~mask_full) = NaN;

    % 去均值并强制缩放到指定的 RMS (sigma_n)
    tmp = noise_corr(mask_full);
    tmp = tmp - mean(tmp, 'omitnan');
    rms0 = sqrt(mean(tmp.^2, 'omitnan'));
    if rms0 <= 0
        error('相关噪声 RMS 为零，无法归一化。');
    end

    Noise = nan(size(noise_corr));
    Noise(mask_full) = sigma_n * tmp / rms0;

    %% 5. 叠加噪声生成含噪面形
    Z_noisy = Z_clean;
    Z_noisy(mask_full) = Z_clean(mask_full) + Noise(mask_full);

    %% 6. 计算布局并提取子孔径
    sep = solve_circle_sep_from_overlap(R_sub, overlap_target);
    u_layout = [cos(phi_layout); sin(phi_layout)];
    cA = -0.5 * sep * u_layout;
    cB =  0.5 * sep * u_layout;

    mask_valid = mask_full & isfinite(Z_noisy);
    F_noisy = scatteredInterpolant(Xm(mask_valid), Ym(mask_valid), ...
                                   Z_noisy(mask_valid), 'natural', 'none');

    subA = sample_circular_subap_local(F_noisy, cA, 0,          R_sub, ds);
    subB = sample_circular_subap_local(F_noisy, cB, theta_true, R_sub, ds);

    %% 7. 打包 Ground Truth 数据供外部计算误差
    truth_data = struct();
    truth_data.Z_clean = Z_clean;
    truth_data.Noise   = Noise;
    truth_data.mask    = mask_full;
    % 真实的相对位姿 (B 相对于 A)
    truth_data.true_tx = cB(1) - cA(1);
    truth_data.true_ty = cB(2) - cA(2);
    truth_data.true_theta = theta_true;
    truth_data.pv_clean = max(Z_clean(mask_full)) - min(Z_clean(mask_full));
    truth_data.rmse_noise_actual = sqrt(mean(Noise(mask_full).^2, 'omitnan'));

end

%% =========================================================================
% 以下为局部辅助函数 (与原脚本完全一致)
%% =========================================================================
function coef = fit_plane_ls(x, y, z)
    A = [x(:), y(:), ones(numel(x),1)];
    coef = A \ z(:);
end

function d = solve_circle_sep_from_overlap(R, target_eta)
    lo = 0;
    hi = 2*R;
    for it = 1:80
        mid = 0.5 * (lo + hi);
        eta_mid = circle_overlap_ratio(mid, R);
        if eta_mid > target_eta
            lo = mid;
        else
            hi = mid;
        end
    end
    d = 0.5 * (lo + hi);
end

function eta = circle_overlap_ratio(d, R)
    if d >= 2*R
        eta = 0;
        return;
    end
    if d <= 0
        eta = 1;
        return;
    end
    part1 = 2 * R^2 * acos(d/(2*R));
    part2 = 0.5 * d * sqrt(4*R^2 - d^2);
    eta = (part1 - part2) / (pi * R^2);
end

function subap = sample_circular_subap_local(Finterp, c, theta, Rsub, ds)
    xv = -Rsub:ds:Rsub;
    yv = -Rsub:ds:Rsub;
    [Xl, Yl] = meshgrid(xv, yv);
    mask = (Xl.^2 + Yl.^2 <= Rsub^2);

    xl = Xl(mask);
    yl = Yl(mask);

    ct = cos(theta);
    st = sin(theta);

    xw = c(1) + ct.*xl - st.*yl;
    yw = c(2) + st.*xl + ct.*yl;
    zw = Finterp(xw, yw);

    valid = isfinite(zw);

    subap = struct();
    subap.type  = 'circle';
    subap.c     = c(:);
    subap.theta = theta;
    subap.Rsub  = Rsub;
    subap.x     = xl(valid);
    subap.y     = yl(valid);
    subap.z     = zw(valid);
end