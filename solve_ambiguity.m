function [N_fixed, corrected_frac_phase, phase_bias] = solve_ambiguity(water_level, sinE, station_alt, lambda)

% solve_ambiguity computes the fixed integer ambiguity by removing systematic phase bias
%
% inputs:
% water_level: vector of water levels from tide gauge (m) [Nx1]
% sinE: vector of sine of satellite elevation angles [Nx1]
% station_alt: station antenna altitude (m)
% lambda: signal wavelength (m)
% outputs:
% N_fixed: Bias-corrected Integer Ambiguity (integer)
% corrected_frac_phase: Residual phase after removing bias (should be near 0)
% phase_bias: Estimated systematic bias (cycles)

    % 1. Calculate Geometry (Reflector Height)
    % H_ref = Antenna Height - Water Level
    H_ref = station_alt - water_level;
    
    % 2. Calculate Theoretical Total Phase (Float)
    % path_delay = 2 * H * sin(E)
    path_delay = 2 * H_ref .* sinE;
    total_phase_float = path_delay / lambda;
    
    % 3. Estimate Systematic Phase Bias
    % Calculate the raw fractional residuals [-0.5, 0.5]
    % This represents how far the Gauge-derived phase is from the nearest integer
    raw_residuals = total_phase_float - round(total_phase_float);
    
    % The median of these residuals is the "Systematic Bias"
    % (e.g., hardware delay, gauge zero-point offset)
    phase_bias = median(raw_residuals, 'omitnan');
    
    % 4. Correct and Fix Ambiguity
    % Remove the bias to shift the phase to the "integer center"
    % This prevents the "sawtooth" jumping when residuals are near 0.5
    total_phase_corrected = total_phase_float - phase_bias;
    
    % Safe Rounding (Fixing)
    N_fixed = round(total_phase_corrected);
    
    % 5. Calculate Corrected Residual Phase
    % This is useful for diagnostics (should now be centered at 0, not 0.3 or -0.4)
    corrected_frac_phase = total_phase_corrected - N_fixed;

end