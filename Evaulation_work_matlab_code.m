%% DCEC - Trabajo de evaluacion Bloque 1: filtro de capacidades conmutadas de 1er orden
clear; clc; close all;
%% Code made by: Jalbatron, 2/10/2026

%% Datos
fclk = 100e3;                       % Frecuencia de reloj [Hz]
C1 = 2e-9; C2 = 1e-9; C3 = 1e-9; C4 = 10e-9;
Sigma = C1 + C2 + C3 + C4;

%% H(z) con z^(1/2):  H = (C2*z^(1/2) - C1) / (C4 - (C3+C4)*z)
% Ganancia de lazo abierto A finita (A = Inf -> amplificador ideal)
H = @(f,A) (C2*exp(1j*pi*f/fclk) - C1) ./ ...
    ( C4*(1+1/A) - ((C3+C4) + Sigma/A).*exp(1j*2*pi*f/fclk) );

%% Ganancia DC, polo y cero (caso ideal)
H0     = abs(H(1e-6, Inf));                 % ganancia DC (z = 1)
zp     = C4/(C3+C4);                        % polo en z
zc     = (C1/C2)^2;                         % cero en z (w = z^(1/2) = C1/C2)
fp     = -log(zp)*fclk/(2*pi);              % frecuencia equivalente del polo
fprintf('Ganancia DC = %.4f (%.2f dB)\n', H0, 20*log10(H0));
fprintf('Polo z = %.5f  (f equiv. = %.1f Hz)\n', zp, fp);
fprintf('Cero z = %.1f (w = z^(1/2) = %.1f) -> fuera del circulo unidad\n', zc, C1/C2);

%% Respuesta en frecuencia
f  = logspace(1, log10(fclk/2), 4000);      % 10 Hz ... 50 kHz
Hi = H(f, Inf);                             % ideal
Hf = H(f, 10);                              % A = 10 (apartado 6)

figure('Name','Respuesta en frecuencia','Color','w');
subplot(2,1,1);
semilogx(f, 20*log10(abs(Hi)), 'b', 'LineWidth', 1.6); hold on;
semilogx(f, 20*log10(abs(Hf)), 'r--', 'LineWidth', 1.4);
grid on; ylabel('|H| [dB]'); title('Filtro SC 1er orden, f_{clk} = 100 kHz');
legend('A = \infty','A = 10','Location','southwest');
subplot(2,1,2);
semilogx(f, unwrap(angle(Hi))*180/pi, 'b', 'LineWidth', 1.6); hold on;
semilogx(f, unwrap(angle(Hf))*180/pi, 'r--', 'LineWidth', 1.4);
grid on; xlabel('f [Hz]'); ylabel('Fase [grados]');

%% Frecuencia de corte (-3 dB respecto a la ganancia en banda = ganancia DC)
for A = [Inf 10]
    G0 = abs(H(1e-6, A));
    fc = fzero(@(x) abs(H(x,A)) - G0/sqrt(2), [10 fclk/4]);
    fprintf('A = %g:  ganancia DC = %.4f,  fc(-3 dB) = %.1f Hz\n', A, G0, fc);
    if isinf(A)
        subplot(2,1,1); plot(fc, 20*log10(G0/sqrt(2)), 'ko', 'MarkerFaceColor','y');
    end
end

%% Modulos y fases en 100 Hz, 1 kHz y 10 kHz
fs = [100 1e3 10e3];
for k = 1:numel(fs)
    h = H(fs(k), Inf);
    fprintf('f = %6.0f Hz: |H| = %.4f (%.2f dB), fase = %.1f grados\n', ...
        fs(k), abs(h), 20*log10(abs(h)), angle(h)*180/pi);
end

%% Muestras con entrada DC (Vi = 1 V, Vo(0) = 0 V)
% (C3+C4+Sigma/A)*Vo(n) = C4*(1+1/A)*Vo(n-1) + (C1-C2)*Vi
Vi = 1; N = 5;
for A = [Inf 10]
    a = (C3+C4) + Sigma/A;  b = C4*(1+1/A);
    Vo = zeros(1,N+1);
    for n = 1:N
        Vo(n+1) = (b*Vo(n) + (C1-C2)*Vi)/a;
    end
    fprintf('A = %g: Vo(1..5) = %s V\n', A, mat2str(Vo(2:end), 5));
end
