%% ÖDEV 4 - Tasarım Noktası Performans Karşılaştırması
% Turbojet, Ramjet ve Ters Brayton Çevrimli Motor

clear; clc; close all;

%% ============================================================
%  SABIT PARAMETRELER (Tablo 4 & 5)
%% ============================================================
gamma   = 1.4;          % Isı sığaları oranı
R       = 287;          % J/(kg.K)
cp      = 1005;         % J/(kg.K)
h_f     = 42.8e6;       % J/kg  (yakıt ısıl değeri)
T0_amb  = 288.15;       % K     (deniz seviyesi statik sıcaklık, ISA)
P0_amb  = 101325;       % Pa    (deniz seviyesi statik basınç, ISA)
m_dot   = 100;          % kg/s  (Soru 2 için hava debisi)

%  İdeal motor + lülede optimum genişleme => P_exit = P_amb

%% ============================================================
%  YARDIMCI FONKSİYONLAR
%% ============================================================
% Toplam sıcaklık: Tt = T_s * (1 + (gamma-1)/2 * M^2)
Tt = @(T_s, M) T_s .* (1 + (gamma-1)/2 .* M.^2);
% Toplam basınç oranı: Pt/P = (Tt/T)^(gamma/(gamma-1))
Pt_ratio = @(M) (1 + (gamma-1)/2 .* M.^2).^(gamma/(gamma-1));

%% ============================================================
%  1A) TURBOJET MOTORU
%% ============================================================
% Parametreler (Tablo 1)
pi_c    = 20;           % Kompresör basınç oranı
T4      = 2000;         % K  (Türbin giriş sıcaklığı)

% Mach aralığı: kalkış (0) → itki = 0 noktası
M_vec_TJ = linspace(0, 5, 500);

Isp_TJ  = zeros(size(M_vec_TJ));
TSFC_TJ = zeros(size(M_vec_TJ));

for k = 1:length(M_vec_TJ)
    M = M_vec_TJ(k);
    
    % Serbest akım toplam sıcaklık & basınç
    Tt0 = Tt(T0_amb, M);
    Pt0_Ps = Pt_ratio(M);            % Pt0 / Ps
    
    % Difüzör (ideal): Tt2 = Tt0, Pt2 = Pt0
    Tt2 = Tt0;
    
    % Kompresör (ideal isentropic)
    Tt3 = Tt2 * pi_c^((gamma-1)/gamma);
    
    % Yanma odası: Tt4 = T4 (sabit)
    if Tt3 >= T4
        % Kompresör çıkışı yanma odası sıcaklığını aşıyor → geçersiz
        Isp_TJ(k)  = 0;
        TSFC_TJ(k) = NaN;
        continue;
    end
    
    % Yakıt-hava oranı
    f = cp*(T4 - Tt3) / h_f;
    
    % Türbin (ideal): türbin işi = kompresör işi (shaft balance)
    delta_Tt_turbine = (Tt3 - Tt2);    % cp*deltaT_c = cp*deltaT_t
    Tt5 = T4 - delta_Tt_turbine;
    
    % Türbin basınç oranı (isentropic)
    pi_t = (Tt5/T4)^(gamma/(gamma-1));
    
    % Lüle: toplam basınç (Pt5 = Pt4 = Pt3 = Pt2 * pi_c)
    % Optimum genişleme: P_exit = P_amb (statik)
    % Tt_lule = Tt5 (ideal lüle)
    % V_exit: h0 = h_exit => cp*Tt5 = cp*T_exit + V_exit^2/2
    % Optimum => P_exit = Ps => T_exit = Ts = Tt5 / (1 + (g-1)/2 * M_exit^2)
    % Daha basit: V_exit^2 = 2*cp*(Tt5 - T0_amb * (Pt5/Ps)^(-(g-1)/g) )
    %              Pt5 / Ps = (Pt2/Ps)*pi_c*pi_t
    Pt5_Ps = Pt0_Ps * pi_c * pi_t;   % Pt5 / P_static_amb
    T_exit = Tt5 / (Pt5_Ps)^((gamma-1)/gamma);
    
    if T_exit < 0; Isp_TJ(k) = 0; TSFC_TJ(k) = NaN; continue; end
    
    V_exit = sqrt(2 * cp * (Tt5 - T_exit));
    
    % Serbest akım hızı
    a0 = sqrt(gamma * R * T0_amb);
    V0 = M * a0;
    
    % Özgül itki [N/(kg/s)]
    Isp_TJ(k) = (1+f)*V_exit - V0;   % momentum thrust / m_dot_air
    % (optimum genişleme → basınç kuvveti = 0)
    
    % Özgül yakıt tüketimi [g/(kN.s)]
    if Isp_TJ(k) > 0
        TSFC_TJ(k) = (f / Isp_TJ(k)) * 1e6;  % kg/(N.s) → g/(kN.s)
    else
        TSFC_TJ(k) = NaN;
        Isp_TJ(k)  = 0;
    end
end

% Sıfır itki Mach sayısını bul
idx_zero_TJ = find(Isp_TJ <= 0, 1, 'first');
if isempty(idx_zero_TJ)
    idx_zero_TJ = length(M_vec_TJ) + 1;
elseif idx_zero_TJ == 1
    idx_zero_TJ = 2;
end
M_vec_TJ_plot  = M_vec_TJ(1:idx_zero_TJ-1);
Isp_TJ_plot    = Isp_TJ(1:idx_zero_TJ-1);
TSFC_TJ_plot   = TSFC_TJ(1:idx_zero_TJ-1);

%% ============================================================
%  1B) RAMJET MOTORU
%% ============================================================
% Parametreler (Tablo 2)
T4_rj = 2000;  % K (yanma odası çıkış sıcaklığı)

% Ramjet kalkışta (M=0) itki üretemez → M_min ~ 0.5 dan başlatılabilir
% ama soru kalkıştan istiyor; M=0 da itki=0 alınır.
M_vec_RJ = linspace(0.01, 8, 800);

Isp_RJ  = zeros(size(M_vec_RJ));
TSFC_RJ = zeros(size(M_vec_RJ));

for k = 1:length(M_vec_RJ)
    M = M_vec_RJ(k);
    
    Tt0 = Tt(T0_amb, M);
    Pt0_Ps = Pt_ratio(M);
    
    % Difüzör (ideal): Tt2 = Tt0
    Tt2 = Tt0;
    
    if Tt2 >= T4_rj
        Isp_RJ(k) = 0; TSFC_RJ(k) = NaN; continue;
    end
    
    % Yakıt-hava oranı
    f = cp*(T4_rj - Tt2) / h_f;
    
    % Lüle (ideal, optimum genişleme)
    % Pt4 = Pt0 (ideal difüzör + ideal yanma odası basınç kaybı yok)
    Pt4_Ps = Pt0_Ps;
    T_exit = T4_rj / Pt4_Ps^((gamma-1)/gamma);
    
    if T_exit < 0; Isp_RJ(k)=0; TSFC_RJ(k)=NaN; continue; end
    
    V_exit = sqrt(2 * cp * (T4_rj - T_exit));
    
    a0 = sqrt(gamma * R * T0_amb);
    V0 = M * a0;
    
    Isp_RJ(k) = (1+f)*V_exit - V0;
    
    if Isp_RJ(k) > 0
        TSFC_RJ(k) = (f / Isp_RJ(k)) * 1e6;
    else
        TSFC_RJ(k) = NaN;
        Isp_RJ(k)  = 0;
    end
end

idx_zero_RJ = find(Isp_RJ <= 0, 1, 'first');
if isempty(idx_zero_RJ)
    idx_zero_RJ = length(M_vec_RJ) + 1;
elseif idx_zero_RJ == 1
    idx_zero_RJ = 2;
end
M_vec_RJ_plot  = M_vec_RJ(1:idx_zero_RJ-1);
Isp_RJ_plot    = Isp_RJ(1:idx_zero_RJ-1);
TSFC_RJ_plot   = TSFC_RJ(1:idx_zero_RJ-1);

%% ============================================================
%  1C) TERS BRAYTON ÇEVRİMLİ MOTOR
%% ============================================================
% Parametreler (Tablo 3)
T_pb     = 1800;  % K  - Ön yakıcı çıkış sıcaklığı
T_ab     = 2000;  % K  - Art yakıcı çıkış sıcaklığı
dTt_cool = 400;   % K  - Soğutma bölümü toplam sıcaklık düşüşü
pi_t_inv = 0.1;   % Türbin basınç oranı

% Doğru Ters Brayton topolojisi:
%  İstasyon 0: Serbest akım
%  İstasyon 2: Difüzör çıkışı  (Tt2 = Tt0, Pt2 = Pt0)
%  İstasyon 3: Ön yakıcı çıkışı (Tt3 = T_pb)
%  İstasyon 4: Türbin çıkışı    (Pt4 = Pt3 * pi_t,  isentropic)
%  İstasyon 5: Soğutma çıkışı   (Tt5 = Tt4 - dTt_cool)
%  İstasyon 6: Kompresör çıkışı (Pt6 = Pt5 / pi_t  → basıncı geri yükle)
%  İstasyon 7: Art yakıcı çıkışı (Tt7 = T_ab)
%  İstasyon 8: Lüle çıkışı      (optimum genişleme)
%
% Shaft balance: Türbin işi = Kompresör işi
%   cp*(Tt3 - Tt4) = cp*(Tt6 - Tt5)
%   → Tt6 = Tt5 + (Tt3 - Tt4)

M_vec_TB = linspace(0, 5, 500);

Isp_TB  = zeros(size(M_vec_TB));
TSFC_TB = zeros(size(M_vec_TB));

for k = 1:length(M_vec_TB)
    M = M_vec_TB(k);
    
    % İst. 0-2: Serbest akım & difüzör (ideal)
    Tt2    = Tt(T0_amb, M);
    Pt2_Ps = Pt_ratio(M);
    
    % İst. 3: Ön yakıcı
    if T_pb < Tt2
        Tt3 = Tt2; f_pb = 0;
    else
        f_pb = cp*(T_pb - Tt2) / h_f;
        Tt3  = T_pb;
    end
    Pt3_Ps = Pt2_Ps;   % ideal yanma odası
    
    % İst. 4: Türbin (isentropic genişleme, pi_t = 0.1)
    Tt4    = Tt3 * pi_t_inv^((gamma-1)/gamma);
    Pt4_Ps = Pt3_Ps * pi_t_inv;
    
    % İst. 5: Soğutma
    Tt5    = Tt4 - dTt_cool;
    Pt5_Ps = Pt4_Ps;
    
    if Tt5 <= 0
        Isp_TB(k) = 0; TSFC_TB(k) = NaN; continue;
    end
    
    % İst. 6: Kompresör (shaft balance ile türbin işine eşit)
    % Türbin işi: w_t = cp*(Tt3 - Tt4)
    % Kompresör: Tt6 = Tt5 + w_t/cp = Tt5 + (Tt3 - Tt4)
    w_t  = cp * (Tt3 - Tt4);
    Tt6  = Tt5 + w_t / cp;
    % Kompresör basınç oranı (isentropic): pi_c = (Tt6/Tt5)^(g/(g-1))
    pi_c_inv = (Tt6/Tt5)^(gamma/(gamma-1));
    Pt6_Ps   = Pt5_Ps * pi_c_inv;
    
    % İst. 7: Art yakıcı
    if T_ab < Tt6
        Tt7 = Tt6; f_ab = 0;
    else
        f_ab = cp*(T_ab - Tt6) / h_f;
        Tt7  = T_ab;
    end
    Pt7_Ps = Pt6_Ps;
    
    f_total = f_pb + f_ab;
    
    % İst. 8: Lüle (optimum genişleme: P_exit = P_amb)
    if Pt7_Ps < 1
        % Basınç atmosferin altında → genişleme yok, itki üretemez
        Isp_TB(k) = 0; TSFC_TB(k) = NaN; continue;
    end
    T_exit = Tt7 / Pt7_Ps^((gamma-1)/gamma);
    V_exit = sqrt(2 * cp * (Tt7 - T_exit));
    
    a0 = sqrt(gamma * R * T0_amb);
    V0 = M * a0;
    
    Isp_TB(k) = (1+f_total)*V_exit - V0;
    
    if Isp_TB(k) > 0
        TSFC_TB(k) = (f_total / Isp_TB(k)) * 1e6;
    else
        TSFC_TB(k) = NaN;
        Isp_TB(k)  = 0;
    end
end

idx_zero_TB = find(Isp_TB <= 0, 1, 'first');
if isempty(idx_zero_TB)
    idx_zero_TB = length(M_vec_TB) + 1;
elseif idx_zero_TB == 1
    idx_zero_TB = 2;
end
M_vec_TB_plot  = M_vec_TB(1:idx_zero_TB-1);
Isp_TB_plot    = Isp_TB(1:idx_zero_TB-1);
TSFC_TB_plot   = TSFC_TB(1:idx_zero_TB-1);

%% ============================================================
%  GRAFİK 1: Özgül İtki vs Mach
%% ============================================================
figure('Name','Özgül İtki ve TSFC Karşılaştırması','Position',[100 100 1100 850]);

subplot(2,1,1);
hold on; grid on;
plot(M_vec_TJ_plot, Isp_TJ_plot, 'b-',  'LineWidth', 2, 'DisplayName','Turbojet');
plot(M_vec_RJ_plot, Isp_RJ_plot, 'r--', 'LineWidth', 2, 'DisplayName','Ramjet');
plot(M_vec_TB_plot, Isp_TB_plot, 'g-.', 'LineWidth', 2, 'DisplayName','Ters Brayton');
xlabel('Mach Sayısı','FontSize',13);
ylabel('Özgül İtki  [N/(kg/s)]','FontSize',13);
title('Özgül İtki – Mach Sayısı İlişkisi','FontSize',14,'FontWeight','bold');
legend('Location','best','FontSize',12);
M_ends = [];
if ~isempty(M_vec_TJ_plot), M_ends(end+1) = M_vec_TJ_plot(end); end
if ~isempty(M_vec_RJ_plot),  M_ends(end+1) = M_vec_RJ_plot(end);  end
if ~isempty(M_vec_TB_plot),  M_ends(end+1) = M_vec_TB_plot(end);  end
if ~isempty(M_ends), xlim([0 max(M_ends)+0.2]); end

subplot(2,1,2);
hold on; grid on;
plot(M_vec_TJ_plot, TSFC_TJ_plot, 'b-',  'LineWidth', 2, 'DisplayName','Turbojet');
plot(M_vec_RJ_plot, TSFC_RJ_plot, 'r--', 'LineWidth', 2, 'DisplayName','Ramjet');
plot(M_vec_TB_plot, TSFC_TB_plot, 'g-.', 'LineWidth', 2, 'DisplayName','Ters Brayton');
xlabel('Mach Sayısı','FontSize',13);
ylabel('TSFC  [g/(kN \cdot s)]','FontSize',13);
title('Özgül Yakıt Tüketimi – Mach Sayısı İlişkisi','FontSize',14,'FontWeight','bold');
legend('Location','best','FontSize',12);
if ~isempty(M_ends), xlim([0 max(M_ends)+0.2]); end
ylim([0 300]);

sgtitle('Motor Performans Karşılaştırması  (Deniz Seviyesi, İdeal)','FontSize',15,'FontWeight','bold');

%% ============================================================
%  GRAFİK 2: Ters Brayton – İtki vs Soğutma Sıcaklık Düşüşü (Kalkış, M=0)
%% ============================================================
dTt_vec   = linspace(0, 800, 400);  % K
Thrust_TB = zeros(size(dTt_vec));

M0 = 0;  % Kalkış

for k = 1:length(dTt_vec)
    dT = dTt_vec(k);
    
    % İst. 0-2: Difüzör (M=0 → Tt2 = T0_amb, Pt2 = Ps)
    Tt2    = Tt(T0_amb, M0);    % = 288.15 K
    Pt2_Ps = Pt_ratio(M0);      % = 1.0
    
    % İst. 3: Ön yakıcı
    if T_pb < Tt2
        Tt3 = Tt2; f_pb = 0;
    else
        f_pb = cp*(T_pb - Tt2) / h_f;
        Tt3  = T_pb;
    end
    Pt3_Ps = Pt2_Ps;
    
    % İst. 4: Türbin (pi_t = 0.1)
    Tt4    = Tt3 * pi_t_inv^((gamma-1)/gamma);
    Pt4_Ps = Pt3_Ps * pi_t_inv;
    
    % İst. 5: Soğutma
    Tt5    = Tt4 - dT;
    Pt5_Ps = Pt4_Ps;
    
    if Tt5 <= 0
        Thrust_TB(k) = 0; continue;
    end
    
    % İst. 6: Kompresör (shaft balance)
    w_t      = cp * (Tt3 - Tt4);
    Tt6      = Tt5 + w_t / cp;
    pi_c_inv = (Tt6/Tt5)^(gamma/(gamma-1));
    Pt6_Ps   = Pt5_Ps * pi_c_inv;
    
    % İst. 7: Art yakıcı
    if T_ab < Tt6
        Tt7 = Tt6; f_ab = 0;
    else
        f_ab = cp*(T_ab - Tt6) / h_f;
        Tt7  = T_ab;
    end
    Pt7_Ps  = Pt6_Ps;
    f_total = f_pb + f_ab;
    
    % İst. 8: Lüle
    if Pt7_Ps < 1
        Thrust_TB(k) = 0; continue;
    end
    T_exit = Tt7 / Pt7_Ps^((gamma-1)/gamma);
    V_exit = sqrt(2 * cp * (Tt7 - T_exit));
    V0_val = 0;  % kalkış
    
    Isp_val = (1+f_total)*V_exit - V0_val;
    
    if Isp_val > 0
        Thrust_TB(k) = Isp_val * m_dot / 1000;  % kN
    else
        Thrust_TB(k) = 0;
    end
end

figure('Name','Ters Brayton - Kalkis Itkisi vs Sogutma DT','Position',[200 150 820 560]);
plot(dTt_vec, Thrust_TB, 'm-', 'LineWidth', 2.5);
grid on;
xlabel('\DeltaT_{sogutma}  [K]','FontSize',13);
ylabel('Kalkış İtkisi  [kN]','FontSize',13);
title('Ters Brayton Çevrimli Motor – Kalkış İtkisi', ...
      'FontSize',14,'FontWeight','bold');
subtitle({'Soğutma Bölümü Toplam Sıcaklık Düşüşüne Bağlı Değişim', ...
          sprintf('\\dot{m} = %d kg/s,  M = 0,  Deniz Seviyesi', m_dot)}, ...
         'FontSize',11);

%% ============================================================
%  KONSOL ÖZET TABLOSU  (Tasarım noktası M=0)
%% ============================================================
fprintf('\n======================================================\n');
fprintf('  KALKIŞ (M=0) TASARIM NOKTASI ÖZET\n');
fprintf('======================================================\n');
fprintf('%-25s %12s %12s\n','Motor','Isp [N/(kg/s)]','TSFC [g/(kN·s)]');
fprintf('------------------------------------------------------\n');

% TJ kalkış
if ~isempty(Isp_TJ_plot)
    fprintf('%-25s %12.1f %12.2f\n','Turbojet', Isp_TJ_plot(1), TSFC_TJ_plot(1));
end
% RJ kalkış (pratik olarak sıfır veya çok düşük)
fprintf('%-25s %12s %12s\n','Ramjet','(kalkışta sıfır)','—');
% TB kalkış
if ~isempty(Isp_TB_plot)
    fprintf('%-25s %12.1f %12.2f\n','Ters Brayton', Isp_TB_plot(1), TSFC_TB_plot(1));
end
fprintf('======================================================\n\n');