function [E_app] = eleCorr(E_input)
% 输入参数: E_input: 真几何高度角 (度，Degrees)
% 输出参数: E_app:   表观折射高度角 (度)

    T_C = 20;        % 当前气温 (摄氏度)
    P_hPa = 1013;    % 当前大气压 (hPa)
    RH_percent = 50; % 相对湿度 (0 到 100)

    % 1. 计算水汽压 e (hPa) 
    es = 6.112 * exp((17.67 * T_C) / (T_C + 243.5)); 
    e = (RH_percent / 100.0) * es;                 
    
    % 2. 计算绝对温度 (Kelvin)
    T_K = T_C + 273.15;
    
    % 3. 计算地表折射率 N_s
    N_dry = 77.6 * (P_hPa / T_K);
    N_wet = 3.73e5 * (e / (T_K^2));
    N_s = N_dry + N_wet;
    
    % 4. 计算折射弯曲角 (注意括号位置！)
    % 算出的是角分 (arcminutes)
    tau_arcmin = 1.02 ./ tand( E_input + (10.3 ./ (E_input + 5.11)) );
    
    % 转成度 (degrees)，并根据实际 N_s 缩放 (标准状态下 N_s 约等于 315)
    dE_deg = tau_arcmin .* (N_s / 315.0) ./ 60.0;
    
    % 5. 得到接收机实际“看到”的高度角
    E_app = E_input + dE_deg;
end