function [P, f] = lomb_scargle_custom(t, y, f_vec)
% LOMB_SCARGLE_CUSTOM 自定义 Lomb-Scargle 周期图计算
%   替代 MATLAB 自带 plomb 函数，解除"严格单调递增"的限制。
%
% 输入:
%   t     - 自变量向量 (对于 GNSS-R 是 sin(E))，无需排序，可包含重复值 [1 x N]
%   y     - 信号向量 (去趋势后的 SNR)，必须与 t 长度相同 [1 x N]
%   f_vec - 需要计算的频率向量 (对应 2*H/lambda) [M x 1]
%
% 输出:
%   P     - 归一化功率谱 (Normalized Power) [M x 1]
%   f     - 对应的频率向量 [M x 1]

    % --- 1. 数据维度标准化 ---
    t = t(:)';       % 强制转为行向量 (1 x N)
    y = y(:)';       % 强制转为行向量 (1 x N)
    f_vec = f_vec(:);% 强制转为列向量 (M x 1)
    
    N = length(y);
    
    % --- 2. 预处理：去均值与计算方差 ---
    % 这一步对应公式中的 (y - y_bar) 和 sigma^2
    y = y - mean(y);
    sigma2 = var(y);
    
    % 保护：如果信号是一条直线（方差为0），直接返回0
    if sigma2 < 1e-12
        P = zeros(size(f_vec));
        f = f_vec;
        return;
    end
    
    % --- 3. 核心计算（完全矩阵化，无循环）---
    % 角频率 omega
    omega = 2 * pi * f_vec; % (M x 1)
    
    % 构建 (2 * omega * t) 矩阵 -> 维度 (M x N)
    % 利用 MATLAB 的隐式广播 (Implicit Expansion) 或矩阵乘法
    two_wt = 2 * omega * t; 
    
    % 计算时间偏移量 tau
    % 公式: tan(2w*tau) = sum(sin(2wt)) / sum(cos(2wt))
    % sum(..., 2) 表示沿着行求和 (把所有时间点加起来)
    S2 = sum(sin(two_wt), 2); % (M x 1)
    C2 = sum(cos(two_wt), 2); % (M x 1)
    
    % 计算 tau = (1/2w) * atan2(S2, C2)
    tau = (1 ./ (2 * omega)) .* atan2(S2, C2); % (M x 1)
    
    % --- 4. 计算功率谱主公式 ---
    % 核心项: omega * (t - tau)
    % t 是 (1xN), tau 是 (Mx1), 我们需要生成 (MxN) 矩阵
    % 使用 bsxfun 或者 MATLAB R2016b+ 的自动广播
    wt_tau = bsxfun(@minus, omega * t, omega .* tau);
    
    % 计算正弦和余弦矩阵
    cos_wt = cos(wt_tau); % (M x N)
    sin_wt = sin(wt_tau); % (M x N)
    
    % 计算分子 (加权和的平方)
    % sum(y .* cos_wt, 2) 相当于矩阵乘法 (cos_wt * y')
    % y 是 (1xN), cos_wt 是 (MxN)
    % 这一步是在做投影
    num_cos = (cos_wt * y').^2; % (M x 1)
    num_sin = (sin_wt * y').^2; % (M x 1)
    
    % 计算分母 (权重的平方和)
    den_cos = sum(cos_wt.^2, 2); % (M x 1)
    den_sin = sum(sin_wt.^2, 2); % (M x 1)
    
    % 防止分母极小值导致除零 (数值稳定性)
    den_cos(den_cos < 1e-10) = 1.0;
    den_sin(den_sin < 1e-10) = 1.0;
    
    % --- 5. 组合结果 ---
    % 标准 Lomb-Scargle 归一化公式
    P = (1 / (2 * sigma2)) * (num_cos ./ den_cos + num_sin ./ den_sin);
    
    f = f_vec;
end