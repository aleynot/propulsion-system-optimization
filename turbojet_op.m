%% İDEAL TURBOJET MOTOR OPTİMİZASYONU
% Ödev 2 - Parametrik ve PSO Optimizasyonu
% =========================================
% Tablo 1: Hava debisi = 100 kg/s
% Tablo 2: H=9 km, M=0.9, gamma=1.4, R=287 J/kgK, cp=1005 J/kgK, LHV=42.8 MJ/kg
% Tablo 3: İdeal motor, lülede optimum genişleme, sabit gaz özellikleri
% Tablo 4: OPR: 2-50, TIT: 1500-2100 K

clc; clear; close all;

%% ========== PARAMETRİK TARAMA ==========
N_OPR   = 200;
N_TIT   = 200;
OPR_vec = linspace(2, 50, N_OPR);
TIT_vec = linspace(1500, 2100, N_TIT);

F_map    = zeros(N_OPR, N_TIT);
TSFC_map = zeros(N_OPR, N_TIT);

for i = 1:N_OPR
    for j = 1:N_TIT
        [F_map(i,j), TSFC_map(i,j)] = perf(OPR_vec(i), TIT_vec(j));
    end
end

valid_mask = F_map > 0 & TSFC_map < 1e5;

%% ========== A) MİNİMUM TSFC ==========
fprintf('========== PARAMETRİK OPTİMİZASYON ==========\n\n');

% A1 - Sınırlama yok
TSFC_tmp = TSFC_map; TSFC_tmp(~valid_mask) = inf;
[minTSFC_A1, idx] = min(TSFC_tmp(:));
[iA1, jA1] = ind2sub([N_OPR, N_TIT], idx);
[F_A1, ~]  = perf(OPR_vec(iA1), TIT_vec(jA1));
fprintf('A1 - Min TSFC (Sinirlama Yok):\n');
fprintf('  OPR = %.3f | TIT = %.1f K | TSFC = %.4f g/(kN.s) | F = %.2f kN\n\n',...
    OPR_vec(iA1), TIT_vec(jA1), minTSFC_A1, F_A1/1000);

% A2 - F >= 90 kN
TSFC_tmp = TSFC_map; TSFC_tmp(~valid_mask | F_map < 90e3) = inf;
[minTSFC_A2, idx] = min(TSFC_tmp(:));
[iA2, jA2] = ind2sub([N_OPR, N_TIT], idx);
[F_A2, ~]  = perf(OPR_vec(iA2), TIT_vec(jA2));
fprintf('A2 - Min TSFC (F >= 90 kN):\n');
fprintf('  OPR = %.3f | TIT = %.1f K | TSFC = %.4f g/(kN.s) | F = %.2f kN\n\n',...
    OPR_vec(iA2), TIT_vec(jA2), minTSFC_A2, F_A2/1000);

% A3 - F >= 100 kN
TSFC_tmp = TSFC_map; TSFC_tmp(~valid_mask | F_map < 100e3) = inf;
[minTSFC_A3, idx] = min(TSFC_tmp(:));
[iA3, jA3] = ind2sub([N_OPR, N_TIT], idx);
[F_A3, ~]  = perf(OPR_vec(iA3), TIT_vec(jA3));
fprintf('A3 - Min TSFC (F >= 100 kN):\n');
fprintf('  OPR = %.3f | TIT = %.1f K | TSFC = %.4f g/(kN.s) | F = %.2f kN\n\n',...
    OPR_vec(iA3), TIT_vec(jA3), minTSFC_A3, F_A3/1000);

%% ========== B) MAKSİMUM İTKİ ==========

% B1 - Sınırlama yok
F_tmp = F_map; F_tmp(~valid_mask) = -inf;
[maxF_B1, idx] = max(F_tmp(:));
[iB1, jB1] = ind2sub([N_OPR, N_TIT], idx);
[~, TSFC_B1] = perf(OPR_vec(iB1), TIT_vec(jB1));
fprintf('B1 - Max Itki (Sinirlama Yok):\n');
fprintf('  OPR = %.3f | TIT = %.1f K | F = %.2f kN | TSFC = %.4f g/(kN.s)\n\n',...
    OPR_vec(iB1), TIT_vec(jB1), maxF_B1/1000, TSFC_B1);

% B2 - TSFC < 25 g/(kN.s)
F_tmp = F_map; F_tmp(~valid_mask | TSFC_map >= 25) = -inf;
[maxF_B2, idx] = max(F_tmp(:));
[iB2, jB2] = ind2sub([N_OPR, N_TIT], idx);
[~, TSFC_B2] = perf(OPR_vec(iB2), TIT_vec(jB2));
fprintf('B2 - Max Itki (TSFC < 25 g/(kN.s)):\n');
fprintf('  OPR = %.3f | TIT = %.1f K | F = %.2f kN | TSFC = %.4f g/(kN.s)\n\n',...
    OPR_vec(iB2), TIT_vec(jB2), maxF_B2/1000, TSFC_B2);

% B3 - TSFC < 30 g/(kN.s)
F_tmp = F_map; F_tmp(~valid_mask | TSFC_map >= 30) = -inf;
[maxF_B3, idx] = max(F_tmp(:));
[iB3, jB3] = ind2sub([N_OPR, N_TIT], idx);
[~, TSFC_B3] = perf(OPR_vec(iB3), TIT_vec(jB3));
fprintf('B3 - Max Itki (TSFC < 30 g/(kN.s)):\n');
fprintf('  OPR = %.3f | TIT = %.1f K | F = %.2f kN | TSFC = %.4f g/(kN.s)\n\n',...
    OPR_vec(iB3), TIT_vec(jB3), maxF_B3/1000, TSFC_B3);

%% ========== PSO OPTİMİZASYONU ==========
fprintf('========== PSO OPTIMIZASYONU ==========\n\n');

lb = [2,    1500];
ub = [50,   2100];

% A1 PSO
[pos_A1, ~] = pso(@obj_A1, lb, ub);
[F_pA1, TSFC_pA1] = perf(pos_A1(1), pos_A1(2));
fprintf('PSO A1 - Min TSFC (Sinirlama Yok):\n');
fprintf('  OPR = %.4f | TIT = %.2f K | TSFC = %.4f g/(kN.s) | F = %.2f kN\n\n',...
    pos_A1(1), pos_A1(2), TSFC_pA1, F_pA1/1000);

% A2 PSO
[pos_A2, ~] = pso(@obj_A2, lb, ub);
[F_pA2, TSFC_pA2] = perf(pos_A2(1), pos_A2(2));
fprintf('PSO A2 - Min TSFC (F >= 90 kN):\n');
fprintf('  OPR = %.4f | TIT = %.2f K | TSFC = %.4f g/(kN.s) | F = %.2f kN\n\n',...
    pos_A2(1), pos_A2(2), TSFC_pA2, F_pA2/1000);

% A3 PSO
[pos_A3, ~] = pso(@obj_A3, lb, ub);
[F_pA3, TSFC_pA3] = perf(pos_A3(1), pos_A3(2));
fprintf('PSO A3 - Min TSFC (F >= 100 kN):\n');
fprintf('  OPR = %.4f | TIT = %.2f K | TSFC = %.4f g/(kN.s) | F = %.2f kN\n\n',...
    pos_A3(1), pos_A3(2), TSFC_pA3, F_pA3/1000);

% B1 PSO
[pos_B1, ~] = pso(@obj_B1, lb, ub);
[F_pB1, TSFC_pB1] = perf(pos_B1(1), pos_B1(2));
fprintf('PSO B1 - Max Itki (Sinirlama Yok):\n');
fprintf('  OPR = %.4f | TIT = %.2f K | F = %.2f kN | TSFC = %.4f g/(kN.s)\n\n',...
    pos_B1(1), pos_B1(2), F_pB1/1000, TSFC_pB1);

% B2 PSO
[pos_B2, ~] = pso(@obj_B2, lb, ub);
[F_pB2, TSFC_pB2] = perf(pos_B2(1), pos_B2(2));
fprintf('PSO B2 - Max Itki (TSFC < 25 g/(kN.s)):\n');
fprintf('  OPR = %.4f | TIT = %.2f K | F = %.2f kN | TSFC = %.4f g/(kN.s)\n\n',...
    pos_B2(1), pos_B2(2), F_pB2/1000, TSFC_pB2);

% B3 PSO
[pos_B3, ~] = pso(@obj_B3, lb, ub);
[F_pB3, TSFC_pB3] = perf(pos_B3(1), pos_B3(2));
fprintf('PSO B3 - Max Itki (TSFC < 30 g/(kN.s)):\n');
fprintf('  OPR = %.4f | TIT = %.2f K | F = %.2f kN | TSFC = %.4f g/(kN.s)\n\n',...
    pos_B3(1), pos_B3(2), F_pB3/1000, TSFC_pB3);

%% ========== KARŞILAŞTIRMA TABLOSU ==========
fprintf('=================================================================\n');
fprintf(' KARSILASTIRMA TABLOSU\n');
fprintf('=================================================================\n');
fprintf('%-4s | %-9s | %6s | %8s | %10s | %8s\n',...
    'Durum','Yontem','OPR','TIT [K]','TSFC[g/kNs]','F [kN]');
fprintf('%s\n', repmat('-',1,58));

tbl = {
    'A1','Param', OPR_vec(iA1), TIT_vec(jA1), minTSFC_A1,  F_A1/1000;
    'A2','Param', OPR_vec(iA2), TIT_vec(jA2), minTSFC_A2,  F_A2/1000;
    'A3','Param', OPR_vec(iA3), TIT_vec(jA3), minTSFC_A3,  F_A3/1000;
    'B1','Param', OPR_vec(iB1), TIT_vec(jB1), TSFC_B1,     maxF_B1/1000;
    'B2','Param', OPR_vec(iB2), TIT_vec(jB2), TSFC_B2,     maxF_B2/1000;
    'B3','Param', OPR_vec(iB3), TIT_vec(jB3), TSFC_B3,     maxF_B3/1000;
    'A1','PSO',   pos_A1(1),    pos_A1(2),    TSFC_pA1,    F_pA1/1000;
    'A2','PSO',   pos_A2(1),    pos_A2(2),    TSFC_pA2,    F_pA2/1000;
    'A3','PSO',   pos_A3(1),    pos_A3(2),    TSFC_pA3,    F_pA3/1000;
    'B1','PSO',   pos_B1(1),    pos_B1(2),    TSFC_pB1,    F_pB1/1000;
    'B2','PSO',   pos_B2(1),    pos_B2(2),    TSFC_pB2,    F_pB2/1000;
    'B3','PSO',   pos_B3(1),    pos_B3(2),    TSFC_pB3,    F_pB3/1000;
};

for k = 1:size(tbl,1)
    fprintf('%-4s | %-9s | %6.3f | %8.1f | %10.4f | %8.2f\n',...
        tbl{k,1}, tbl{k,2}, tbl{k,3}, tbl{k,4}, tbl{k,5}, tbl{k,6});
end
fprintf('=================================================================\n');

%% ========== GRAFİKLER ==========
[OPR_g, TIT_g] = meshgrid(OPR_vec, TIT_vec);
F_plt    = F_map'   / 1000;
TSFC_plt = TSFC_map';
TSFC_plt(TSFC_plt > 80) = NaN;

figure('Name','Performans Haritalari','Position',[50,50,1300,500]);

subplot(1,2,1);
contourf(OPR_g, TIT_g, F_plt, 25); colorbar;
xlabel('OPR [-]'); ylabel('TIT [K]'); title('Itki [kN]');
hold on;
plot(OPR_vec(iA1),TIT_vec(jA1),'w*','MarkerSize',14,'LineWidth',2,'DisplayName','A1');
plot(OPR_vec(iB1),TIT_vec(jB1),'ws','MarkerSize',14,'LineWidth',2,'DisplayName','B1');
legend('Location','northwest'); colormap(gca,'jet');

subplot(1,2,2);
contourf(OPR_g, TIT_g, TSFC_plt, 25); colorbar;
xlabel('OPR [-]'); ylabel('TIT [K]'); title('TSFC [g/(kN-s)]');
hold on;
plot(OPR_vec(iA1),TIT_vec(jA1),'w*','MarkerSize',14,'LineWidth',2,'DisplayName','A1');
plot(OPR_vec(iB1),TIT_vec(jB1),'ws','MarkerSize',14,'LineWidth',2,'DisplayName','B1');
legend('Location','northwest'); colormap(gca,'parula');

figure('Name','OPR Taramasi','Position',[50,600,1300,400]);
TIT_fixed = [1500, 1800, 2100];
clrs = {'b','r--','g-.'};

subplot(1,2,1); hold on;
for k = 1:3
    Fv = zeros(1,N_OPR);
    for i = 1:N_OPR
        [Fv(i), ~] = perf(OPR_vec(i), TIT_fixed(k));
    end
    plot(OPR_vec, Fv/1000, clrs{k}, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('TIT=%d K', TIT_fixed(k)));
end
xlabel('OPR'); ylabel('Itki [kN]'); title('OPR – Itki'); legend; grid on;

subplot(1,2,2); hold on;
for k = 1:3
    TSv = zeros(1,N_OPR);
    for i = 1:N_OPR
        [~, TSv(i)] = perf(OPR_vec(i), TIT_fixed(k));
    end
    TSv(TSv > 80) = NaN;
    plot(OPR_vec, TSv, clrs{k}, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('TIT=%d K', TIT_fixed(k)));
end
xlabel('OPR'); ylabel('TSFC [g/(kN-s)]'); title('OPR – TSFC'); legend; grid on;

fprintf('\nTamamlandi.\n');

%% ================================================================
%%  YEREL FONKSİYONLAR  (her biri yalnizca BİR kez tanimlanmistir)
%% ================================================================

function [F, TSFC] = perf(OPR, TIT)
%PERF  Ideal turbojet performans modeli
%  Girdi : OPR – kompresör basınç oranı, TIT – türbin giris sicakligi [K]
%  Cikti : F   – itki [N],  TSFC – özgül yakit tüketimi [g/(kN·s)]

    gamma_ = 1.4;
    cp_    = 1005;       % J/(kg·K)
    LHV_   = 42.8e6;    % J/kg
    mdot_  = 100;        % kg/s
    R_     = 287;        % J/(kg·K)
    M_     = 0.9;
    T_amb  = 229.74;     % ISA @ 9 km [K]
    P_amb  = 30742;      % ISA @ 9 km [Pa]

    % Giris stagnasyon kosullari
    T0in = T_amb * (1 + (gamma_-1)/2 * M_^2);
    P0in = P_amb * (T0in/T_amb)^(gamma_/(gamma_-1));

    % Ucus hizi
    V0 = M_ * sqrt(gamma_ * R_ * T_amb);

    % Kompresör cikisi
    T03 = T0in * OPR^((gamma_-1)/gamma_);
    P03 = P0in * OPR;

    % Yakit-hava orani
    f = cp_ * (TIT - T03) / LHV_;
    if f <= 0
        F = 0; TSFC = 1e6; return;
    end

    % Türbin (is dengesi, ideal)
    T05 = TIT - (T03 - T0in);
    if T05 <= 0
        F = 0; TSFC = 1e6; return;
    end

    % Türbin cikis basinci
    pi_t = (T05/TIT)^(gamma_/(gamma_-1));
    P05  = P03 * pi_t;

    % Nozzle – optimum genisleme (Pe = P_amb)
    if P05 <= P_amb
        F = 0; TSFC = 1e6; return;
    end
    T9 = T05 * (P_amb/P05)^((gamma_-1)/gamma_);
    if T9 <= 0 || T05 <= T9
        F = 0; TSFC = 1e6; return;
    end
    V9 = sqrt(2 * cp_ * (T05 - T9));

    % Performans
    F = mdot_ * ((1+f)*V9 - V0);
    if F <= 0
        F = 0; TSFC = 1e6; return;
    end
    TSFC = (mdot_*f) / (F/1000) * 1000;   % g/(kN·s)
end

% ----------------------------------------------------------------
function [best_pos, best_val] = pso(obj_fun, lb, ub)
%PSO  Parcacik Suru Optimizasyonu
    n_p   = 60;
    n_it  = 300;
    w     = 0.7;
    c1    = 1.5;
    c2    = 1.5;
    n_dim = length(lb);

    pos = lb + rand(n_p, n_dim) .* (ub - lb);
    vel = zeros(n_p, n_dim);

    pbest_pos = pos;
    pbest_val = zeros(1, n_p);
    for k = 1:n_p
        pbest_val(k) = obj_fun(pos(k,:));
    end

    [best_val, gi] = min(pbest_val);
    best_pos = pbest_pos(gi,:);

    for it = 1:n_it                                         %#ok<FORPF>
        r1  = rand(n_p, n_dim);
        r2  = rand(n_p, n_dim);
        vel = w*vel + c1*r1.*(pbest_pos - pos) + c2*r2.*(best_pos - pos);
        pos = max(min(pos + vel, ub), lb);

        for k = 1:n_p
            v = obj_fun(pos(k,:));
            if v < pbest_val(k)
                pbest_val(k)   = v;
                pbest_pos(k,:) = pos(k,:);
            end
        end

        [iter_best, gi] = min(pbest_val);
        if iter_best < best_val
            best_val = iter_best;
            best_pos = pbest_pos(gi,:);
        end
    end
end

% ----------------------------------------------------------------
%  Amac fonksiyonlari – her biri yalnizca bir kez tanimlanmistir
% ----------------------------------------------------------------
function v = obj_A1(x)          % Min TSFC, sinirlama yok
    [F, TS] = perf(x(1), x(2));
    if F <= 0, v = 1e6; else, v = TS; end
end

function v = obj_A2(x)          % Min TSFC, F >= 90 kN
    [F, TS] = perf(x(1), x(2));
    if F < 90e3, v = 1e6; else, v = TS; end
end

function v = obj_A3(x)          % Min TSFC, F >= 100 kN
    [F, TS] = perf(x(1), x(2));
    if F < 100e3, v = 1e6; else, v = TS; end
end

function v = obj_B1(x)          % Max itki, sinirlama yok
    [F, ~] = perf(x(1), x(2));
    if F <= 0, v = 1e9; else, v = -F; end
end

function v = obj_B2(x)          % Max itki, TSFC < 25
    [F, TS] = perf(x(1), x(2));
    if F <= 0 || TS >= 25, v = 1e9; else, v = -F; end
end

function v = obj_B3(x)          % Max itki, TSFC < 30
    [F, TS] = perf(x(1), x(2));
    if F <= 0 || TS >= 30, v = 1e9; else, v = -F; end
end