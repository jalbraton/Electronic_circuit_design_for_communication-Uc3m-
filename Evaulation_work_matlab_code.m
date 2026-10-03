%% DCEC - Evaluation Work Block 1: 1st Order Switched-Capacitor Filter
% Author: Jalbatron
% Date: October 2026
% Description: Matlab script to analyze the frequency response, pole/zero 
% locations, and time-domain transient response of a 1st order SC filter.

clear; clc; close all;

%% 1. Circuit Specifications & Parameters
fclk = 100e3;                         % Clock frequency: 100 kHz
T = 1/fclk;                           % Clock period: 10 us

% Given component values (in Farads)
C1 = 2e-9;                            % 2 nF
C2 = 1e-9;                            % 1 nF
C3 = 1e-9;                            % 1 nF
C4 = 10e-9;                           % 10 nF

% Sum of all capacitors (needed for non-ideal op-amp modeling)
Sigma = C1 + C2 + C3 + C4;

%% 2. Transfer Function Definition H(f, A)
% Transfer function in the z-domain:
% H(z) = (C2 * z^(1/2) - C1) / ( C4*(1+1/A) - ((C3+C4) + Sigma/A)*z )
% 
% Note: z^(1/2) = exp(j * pi * f / fclk) represents the half-period delay 
% between switching phases (phi1 to phi2).
% A represents the open-loop gain of the op-amp (A = Inf for ideal op-amp).

H = @(f, A) (C2*exp(1j*pi*f/fclk) - C1) ./ ...
    ( C4*(1 + 1/A) - ((C3 + C4) + Sigma/A).*exp(1j*2*pi*f/fclk) );

%% 3. Ideal Case Analysis: Pole, Zero, and DC Gain (A = Inf)
% DC gain evaluated near DC (f -> 0, z -> 1)
H0 = abs(H(1e-6, Inf));

% Discrete pole and zero locations in the z-plane
zp = C4 / (C3 + C4);                  % zp = 10/11 ~ 0.9091 (inside unit circle -> stable)
zc = (C1 / C2)^2;                     % zc = 4 (outside unit circle)

% Equivalent continuous-time pole frequency using zp = exp(-wp * T)
fp = -log(zp) * fclk / (2*pi);

% Display calculated analytical parameters
fprintf('=== 1. Theoretical Analysis (Ideal Op-Amp) ===\n');
fprintf('DC Gain (H0)     = %.4f (%.2f dB)\n', H0, 20*log10(H0));
fprintf('Pole (zp)        = %.5f  --> Equivalent pole freq (fp) = %.1f Hz\n', zp, fp);
fprintf('Zero (zc)        = %.1f     --> Outside unit circle (sqrt(zc) = %.1f)\n\n', zc, C1/C2);

%% 4. Frequency Response Plotting (Ideal vs. Finite Gain A = 10)
f = logspace(1, log10(fclk/2), 4000); % Frequency sweep from 10 Hz to Nyquist limit (50 kHz)

Hi = H(f, Inf);                       % Ideal op-amp response (A = infinity)
Hf = H(f, 10);                        % Non-ideal op-amp response (A = 10)

figure('Name', 'Switched-Capacitor Filter Frequency Response', 'Color', 'w');

% Magnitude Plot
subplot(2, 1, 1);
semilogx(f, 20*log10(abs(Hi)), 'b', 'LineWidth', 1.6); hold on;
semilogx(f, 20*log10(abs(Hf)), 'r--', 'LineWidth', 1.4);
grid on; 
ylabel('Magnitude |H(f)| [dB]');
title('1st Order SC Filter Frequency Response (f_{clk} = 100 kHz)');
legend('Ideal Op-Amp (A = \infty)', 'Finite Gain (A = 10)', 'Location', 'southwest');

% Phase Plot
subplot(2, 1, 2);
semilogx(f, unwrap(angle(Hi))*180/pi, 'b', 'LineWidth', 1.6); hold on;
semilogx(f, unwrap(angle(Hf))*180/pi, 'r--', 'LineWidth', 1.4);
grid on; 
xlabel('Frequency [Hz]'); 
ylabel('Phase [Degrees]');

%% 5. Cut-off Frequency Calculation (-3 dB from In-Band/DC Gain)
fprintf('=== 2. Cut-off Frequency Calculation (-3 dB) ===\n');
for A_val = [Inf, 10]
    G0 = abs(H(1e-6, A_val));
    % Find frequency where magnitude drops by 3 dB (G0 / sqrt(2))
    fc = fzero(@(x) abs(H(x, A_val)) - G0/sqrt(2), [10, fclk/4]);
    
    if isinf(A_val)
        fprintf('Ideal Op-Amp (A = inf): DC Gain = %.4f (%.2f dB), fc(-3 dB) = %.1f Hz\n', G0, 20*log10(G0), fc);
        % Highlight cut-off frequency on the magnitude plot
        subplot(2, 1, 1);
        plot(fc, 20*log10(G0/sqrt(2)), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
    else
        fprintf('Finite Gain  (A = 10) : DC Gain = %.4f (%.2f dB), fc(-3 dB) = %.1f Hz\n\n', G0, 20*log10(G0), fc);
    end
end

%% 6. Magnitude and Phase at Specific Test Frequencies
fprintf('=== 3. Frequency Evaluation at 100 Hz, 1 kHz, and 10 kHz (Ideal Case) ===\n');
test_freqs = [100, 1e3, 10e3];

for k = 1:numel(test_freqs)
    fk = test_freqs(k);
    hk = H(fk, Inf);
    fprintf('f = %5.0f Hz:  |H| = %.4f (%6.2f dB), Phase = %6.1f deg\n', ...
        fk, abs(hk), 20*log10(abs(hk)), angle(hk)*180/pi);
end
fprintf('\n');

%% 7. Discrete-Time Transient Response (First 5 Output Samples)
% Difference equation derived from charge conservation:
% ((C3 + C4) + Sigma/A) * Vo[n] = C4*(1 + 1/A) * Vo[n-1] + (C1 - C2) * Vi

Vi = 1;                                % Step input voltage: 1 V DC
N_samples = 5;                         % Number of samples to calculate

fprintf('=== 4. First 5 Output Samples (Vi = 1V DC, Vo(0) = 0V) ===\n');
for A_val = [Inf, 10]
    a = (C3 + C4) + Sigma/A_val;
    b = C4 * (1 + 1/A_val);
    
    Vo = zeros(1, N_samples + 1);       % Initialize output vector (Vo(0) = 0 V)
    
    for n = 1:N_samples
        Vo(n+1) = (b * Vo(n) + (C1 - C2) * Vi) / a;
    end
    
    if isinf(A_val)
        fprintf('Ideal Op-Amp (A = inf) : Vo(1..5) = %s V\n', mat2str(Vo(2:end), 5));
    else
        fprintf('Finite Gain  (A = 10)  : Vo(1..5) = %s V\n', mat2str(Vo(2:end), 5));
    end
end