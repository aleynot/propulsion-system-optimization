%% =========================================================
%  İTKİ SİSTEMLERİ PERFORMANS OPTİMİZASYONU - PROJE 2
%  3 Şaftlı Turbofan Motor PSO Optimizasyonu
%
%  ÇEVRIM MODELİ (GasTurb notasyonu):
%  --------------------------------------------------------
%  İstasyon 2  : Giriş çıkışı (inlet exit)
%  İstasyon 13 : Fan bypass çıkışı
%  İstasyon 21 : Fan core çıkışı (LPC girişi)
%  İstasyon 25 : LPC çıkışı (IPC girişi) -- 3 şaftlı: LPC=IPC
%  İstasyon 3  : HPC çıkışı (yanma odası girişi)
%  İstasyon 4  : Türbin girişi (yanma odası çıkışı)
%  İstasyon 45 : IPT (orta basınç türbini) çıkışı
%  İstasyon 5  : LPT (düşük basınç türbini) çıkışı
%  İstasyon 8  : Core lüle çıkışı
%  İstasyon 18 : Bypass lüle çıkışı
%
%  3 ŞAFT:
%   - LP şaftı  : Fan + LPT
%   - IP şaftı  : LPC + IPT
%   - HP şaftı  : HPC + HPT
%
%  GAZ ÖZELLİKLERİ:
%   Soğuk bölge (fan, komp.): gamma_c=1.4, cp_c=1005 J/(kg.K)
%   Sıcak bölge (türbin, yanma): gamma_h=1.333, cp_h=1150 J/(kg.K)
%
%  TASARIM PARAMETRELERİ (6 adet):
%   x(1) = eta_fan_p   Fan politropik verimi    [0.80, 0.95]  (kısıt <1 birim değil - bkz. not)
%   x(2) = PR_fan      Fan basınç oranı         [1.2, 2.0]    (kısıt <3)
%   x(3) = PR_LPC      LPC basınç oranı         [1.5, 10]
%   x(4) = PR_HPC      HPC basınç oranı         [2.0, 10]
%   x(5) = TIT_TO      Türbin giriş sıcaklığı kalkış [1400, 2000] K
%   x(6) = BPR         Bypass oranı             [3, 12]
%
%  NOT: Tabloda "Fan politropik verimi < 3" yazıyor.
%       Bu büyük ihtimalle "Fan basınç oranı < 3" olmalı (typo),
%       ya da "Fan politropik verimi" yanlış hücreye girmiş.
%       Mevcut tabloya göre fan politropik verimi sabit=0.89 (Tablo 6),
%       fan basınç oranı < 3 kısıtı uygulanacak.
%% =========================================================
clear; clc; close all;

%% --- Sabitler ---
gamma_c = 1.4;   cp_c = 1005;   % Soğuk bölge
gamma_h = 1.333; cp_h = 1150;   % Sıcak bölge
R       = 287;
LHV     = 42.8e6;

% Verimler (Tablo 6)
eta_inlet  = 0.98;   % Hava alığı basınç oranı
eta_comb   = 0.985;  % Yanma verimi
PR_comb    = 0.97;   % Yanma odası basınç oranı
eta_fan_p  = 0.89;   % Fan politropik verimi (SABİT - Tablo 6)
eta_LPC_p  = 0.90;   % LPC politropik verimi
eta_HPC_p  = 0.91;   % HPC politropik verimi
eta_LPT_p  = 0.93;   % LPT politropik verimi
eta_IPT_p  = 0.91;   % IPT (orta basınç) politropik verimi
eta_HPT_p  = 0.90;   % HPT politropik verimi
eta_mHP    = 0.99;   % HP şaft mekanik verimi
eta_mIP    = 0.99;   % IP şaft mekanik verimi
eta_mLP    = 0.995;  % LP şaft mekanik verimi
PR_byp_noz = 0.98;   % Bypass lüle basınç oranı
PR_cor_noz = 0.98;   % Core lüle basınç oranı

%% --- Kalkış koşulları (SL, M=0) ---
T0_TO  = 288.15;
P0_TO  = 101325;
M_TO   = 0;
T02_TO = T0_TO * (1 + (gamma_c-1)/2*M_TO^2);
P02_TO = P0_TO * (T02_TO/T0_TO)^(gamma_c/(gamma_c-1));
T02_TO = T02_TO * eta_inlet + T02_TO*(1-eta_inlet);  % Basitleştirilmiş
P02_TO = P02_TO * eta_inlet;
mdot_TO = 1000;   % kg/s (Tablo 4)

%% --- Seyir koşulları (11 km, M=0.85) ---
T0_cr  = 216.65;   % ISA 11 km [K]
P0_cr  = 22632;    % ISA 11 km [Pa]
M_cr   = 0.85;
a_cr   = sqrt(gamma_c*R*T0_cr);
V_cr   = M_cr*a_cr;
T02_cr = T0_cr*(1+(gamma_c-1)/2*M_cr^2);
P02_cr = P0_cr*(T02_cr/T0_cr)^(gamma_c/(gamma_c-1));
P02_cr = P02_cr * eta_inlet;
T02_cr_total = T02_cr;

fprintf('=== GİRİŞ KOŞULLARI ===\n');
fprintf('Kalkış: T02=%.2f K, P02=%.1f Pa\n',T02_TO,P02_TO);
fprintf('Seyir : T02=%.2f K, P02=%.1f Pa\n',T02_cr,P02_cr);
fprintf('V_cr  = %.1f m/s\n\n',V_cr);

%% Motor ağırlık modeli (Torenbeek/Jenkinson yaklaşımı)
% W_engine [kg] ~ k * (mdot)^0.9 * (OPR)^0.4 * (TIT/1000)^(-0.5)
% Basit model: W ~ 0.084 * mdot_TO^0.75 * OPR^0.5
% Lüle alanı modeli: A ~ mdot / (rho * V)_exit

%% =========================================================
%  PSO AYARLARI
%% =========================================================
N_par  = 120;
N_iter = 500;
w_max  = 0.9; w_min = 0.4;
c1 = 1.6; c2 = 1.6;
n_var = 6;

% [PR_fan, PR_LPC, PR_HPC, TIT_TO, BPR, TIT_cr_frac]
%  TIT_cr = TIT_TO * TIT_cr_frac (<1920 K kısıtı)
lb = [1.2,  1.5,  2.0,  1400, 3,  0.7];
ub = [2.99, 9.99, 9.99, 2000, 12, 1.0];

rng(42);
X = lb + rand(N_par,n_var).*(ub-lb);
V = 0.05*(rand(N_par,n_var)-0.5).*(ub-lb);
pBest = X; pBest_f = inf(N_par,1);

for i=1:N_par
    pBest_f(i) = turbofan(X(i,:), T02_TO, P02_TO, T02_cr, P02_cr, ...
        mdot_TO, gamma_c, gamma_h, cp_c, cp_h, LHV, V_cr, P0_TO, P0_cr, ...
        eta_fan_p, eta_LPC_p, eta_HPC_p, eta_HPT_p, eta_IPT_p, eta_LPT_p, ...
        eta_mHP, eta_mIP, eta_mLP, eta_comb, PR_comb, PR_byp_noz, PR_cor_noz);
end
[gBest_f,idx] = min(pBest_f); gBest = pBest(idx,:);
hist_f = nan(N_iter,1); hist_x = nan(N_iter,n_var);

fprintf('%-8s %-14s %-12s %-12s %-10s %-8s\n','İter','SFC','İtki_TO','İtki_CR','BPR','OPR');
fprintf('%s\n',repmat('-',68,1));

for iter=1:N_iter
    w = w_max-(w_max-w_min)*iter/N_iter;
    r1=rand(N_par,n_var); r2=rand(N_par,n_var);
    V = w*V+c1*r1.*(pBest-X)+c2*r2.*(repmat(gBest,N_par,1)-X);
    Vmax=0.25*(ub-lb); V=max(min(V,Vmax),-Vmax);
    X=X+V;
    for d=1:n_var
        lo=X(:,d)<lb(d); hi=X(:,d)>ub(d);
        X(lo,d)=lb(d)+rand(sum(lo),1).*(ub(d)-lb(d))*0.05;
        X(hi,d)=ub(d)-rand(sum(hi),1).*(ub(d)-lb(d))*0.05;
        V(lo,d)=abs(V(lo,d)); V(hi,d)=-abs(V(hi,d));
    end
    for i=1:N_par
        fi = turbofan(X(i,:), T02_TO, P02_TO, T02_cr, P02_cr, ...
            mdot_TO, gamma_c, gamma_h, cp_c, cp_h, LHV, V_cr, P0_TO, P0_cr, ...
            eta_fan_p, eta_LPC_p, eta_HPC_p, eta_HPT_p, eta_IPT_p, eta_LPT_p, ...
            eta_mHP, eta_mIP, eta_mLP, eta_comb, PR_comb, PR_byp_noz, PR_cor_noz);
        if fi<pBest_f(i); pBest_f(i)=fi; pBest(i,:)=X(i,:); end
        if fi<gBest_f;    gBest_f=fi;    gBest=X(i,:);       end
    end
    hist_f(iter)=gBest_f; hist_x(iter,:)=gBest;
    if mod(iter,50)==0||iter==1
        [~,tTO,tCR,sfc_val] = turbofan(gBest, T02_TO, P02_TO, T02_cr, P02_cr, ...
            mdot_TO, gamma_c, gamma_h, cp_c, cp_h, LHV, V_cr, P0_TO, P0_cr, ...
            eta_fan_p, eta_LPC_p, eta_HPC_p, eta_HPT_p, eta_IPT_p, eta_LPT_p, ...
            eta_mHP, eta_mIP, eta_mLP, eta_comb, PR_comb, PR_byp_noz, PR_cor_noz);
        OPR = gBest(1)*gBest(2)*gBest(3);
        fprintf('%-8d %-14.4e %-12.2f %-12.2f %-10.2f %-8.1f\n',...
                iter,sfc_val,tTO/1e3,tCR/1e3,gBest(5),OPR);
    end
end

[~,thr_TO,thr_CR,sfc_opt,~,W_eng,A_total,TW,OPR_opt] = turbofan(gBest, ...
    T02_TO, P02_TO, T02_cr, P02_cr, ...
    mdot_TO, gamma_c, gamma_h, cp_c, cp_h, LHV, V_cr, P0_TO, P0_cr, ...
    eta_fan_p, eta_LPC_p, eta_HPC_p, eta_HPT_p, eta_IPT_p, eta_LPT_p, ...
    eta_mHP, eta_mIP, eta_mLP, eta_comb, PR_comb, PR_byp_noz, PR_cor_noz);

TIT_TO = gBest(4);
TIT_cr = gBest(4)*gBest(6);

fprintf('\n=== OPTİMUM TASARIM PARAMETRELERİ ===\n');
nms={'PR_fan  [-]','PR_LPC  [-]','PR_HPC  [-]','TIT_TO  [K]','BPR     [-]','TIT_cr/TIT_TO'};
for k=1:n_var; fprintf('  %-22s = %10.4f\n',nms{k},gBest(k)); end
fprintf('  %-22s = %10.4f\n','TIT_cr  [K]',TIT_cr);
fprintf('  %-22s = %10.4f\n','OPR     [-]',OPR_opt);

fprintf('\n=== PERFORMANS ===\n');
fprintf('  SFC (seyir)           = %.4e kg/(N.s)\n',sfc_opt);
fprintf('  İtki (kalkış)         = %.2f kN\n',thr_TO/1e3);
fprintf('  İtki (seyir)          = %.2f kN\n',thr_CR/1e3);
fprintf('  Motor ağırlığı        = %.2f ton\n',W_eng/1000);
fprintf('  İtki/ağırlık oranı    = %.2f\n',TW);
fprintf('  Lüle alanları A8+A18  = %.3f m2\n',A_total);

fprintf('\n=== KISIT KONTROL ===\n');
ok=@(c)ternary(c,'✓ SAĞLANDI','✗ SAĞLANMADI');
fprintf('  İtki_TO > 350 kN      : %s (%.1f kN)\n',ok(thr_TO>350e3),thr_TO/1e3);
fprintf('  İtki_CR > 90 kN       : %s (%.1f kN)\n', ok(thr_CR>90e3),thr_CR/1e3);
fprintf('  TW > 5                : %s (%.2f)\n',     ok(TW>5),TW);
fprintf('  W < 6 ton             : %s (%.2f ton)\n', ok(W_eng<6000),W_eng/1000);
fprintf('  A8+A18 < 6 m2         : %s (%.3f m2)\n',  ok(A_total<6),A_total);
fprintf('  OPR < 70 (kalkış)     : %s (%.1f)\n',     ok(OPR_opt<70),OPR_opt);
fprintf('  TIT_TO < 2000 K       : %s (%.1f K)\n',   ok(TIT_TO<2000),TIT_TO);
fprintf('  TIT_cr < 1920 K       : %s (%.1f K)\n',   ok(TIT_cr<1920),TIT_cr);
fprintf('  PR_fan < 3            : %s (%.2f)\n',      ok(gBest(1)<3),gBest(1));

%% GRAFİKLER
figure('Name','PSO Yakınsama','Position',[50 50 1400 900]);
subplot(3,3,1); plot(1:N_iter,hist_f,'b-','LineWidth',1.8);
xlabel('İterasyon'); ylabel('SFC [kg/(N·s)]'); title('Global En İyi SFC'); grid on;
lbl={'PR_{fan}','PR_{LPC}','PR_{HPC}','TIT_{TO} [K]','BPR','TIT_{cr}/TIT_{TO}'};
clr=lines(6);
for k=1:6
    subplot(3,3,k+1);
    plot(1:N_iter,hist_x(:,k),'Color',clr(k,:),'LineWidth',1.5);
    xlabel('İterasyon'); ylabel(lbl{k}); title(lbl{k}); grid on;
end
sgtitle('PSO Yakınsama Grafikleri','FontSize',12,'FontWeight','bold');

figure('Name','Motor Performans','Position',[100 100 1000 420]);
subplot(1,2,1);
b=bar([thr_TO/1e3,thr_CR/1e3],'FaceColor','flat');
b.CData=[0.2 0.4 0.8;0.1 0.6 0.3];
set(gca,'XTickLabel',{'Kalkış','Seyir'});
ylabel('İtki [kN]'); title('Motor İtkisi'); hold on;
yline(350,'r--','LineWidth',1.5,'Label','350 kN');
yline(90,'g--','LineWidth',1.5,'Label','90 kN'); grid on;
subplot(1,2,2);
v=~isnan(hist_f)&hist_f<1e9;
if any(v); semilogy(find(v),hist_f(v),'m-','LineWidth',2); end
xlabel('İterasyon'); ylabel('SFC (log)'); title('SFC Yakınsaması'); grid on;
sgtitle('Motor Performans Özeti','FontSize',12,'FontWeight','bold');
fprintf('\nBitti.\n');

%% =========================================================
%  TURBOFAN FİZİĞİ FONKSİYONU
%% =========================================================
function [f, thrust_TO, thrust_CR, sfc, pen_out, W_eng, A_total, TW, OPR] = ...
    turbofan(x, T02_TO, P02_TO, T02_cr, P02_cr, mdot_TO, ...
             gamma_c, gamma_h, cp_c, cp_h, LHV, V_cr, P0_TO, P0_cr, ...
             eta_fan_p, eta_LPC_p, eta_HPC_p, eta_HPT_p, eta_IPT_p, eta_LPT_p, ...
             eta_mHP, eta_mIP, eta_mLP, eta_comb, PR_comb, PR_byp_noz, PR_cor_noz)

    f=1e10; thrust_TO=0; thrust_CR=0; sfc=inf; pen_out=0;
    W_eng=inf; A_total=inf; TW=0; OPR=0;

    PR_fan   = x(1);
    PR_LPC   = x(2);
    PR_HPC   = x(3);
    TIT_TO   = x(4);
    BPR      = x(5);
    TIT_frac = x(6);   % TIT_cr = TIT_TO * TIT_frac
    TIT_cr   = TIT_TO * TIT_frac;

    pen=0; BIG=1e6;
    n_c = gamma_c/(gamma_c-1);  % = 3.5
    n_h = gamma_h/(gamma_h-1);  % ≈ 4.0

    % Politropik → isentropik verim dönüşümü için üsler
    % Kompresör: eta_is = (PR^((g-1)/g) - 1) / (PR^((g-1)/(g*eta_p)) - 1)
    % Türbin:    eta_is = (1 - PR^(eta_p*(g-1)/g)) / (1 - PR^((g-1)/g))
    % Sıcaklık oranı kompresör: T_out/T_in = PR^((g-1)/(g*eta_p))
    % Sıcaklık oranı türbin:    T_out/T_in = PR^(eta_p*(g-1)/g)

    %% ====== KALKI$ MODU ======
    % --- Fan ---
    tau_fan = PR_fan^((gamma_c-1)/(gamma_c*eta_fan_p));
    T13_TO  = T02_TO * tau_fan;   % Bypass akışı
    T21_TO  = T13_TO;             % Core fan çıkışı (aynı PR)
    P13_TO  = P02_TO * PR_fan;
    P21_TO  = P13_TO;

    % --- LPC ---
    tau_LPC = PR_LPC^((gamma_c-1)/(gamma_c*eta_LPC_p));
    T25_TO  = T21_TO * tau_LPC;
    P25_TO  = P21_TO * PR_LPC;

    % --- HPC ---
    tau_HPC = PR_HPC^((gamma_c-1)/(gamma_c*eta_HPC_p));
    T3_TO   = T25_TO * tau_HPC;
    P3_TO   = P25_TO * PR_HPC;

    OPR = PR_fan * PR_LPC * PR_HPC;
    if OPR > 70; pen=pen+BIG*(OPR-70)/10; end

    % --- Yanma odası ---
    P4_TO = P3_TO * PR_comb;
    if TIT_TO < T3_TO
        pen=pen+BIG*(T3_TO-TIT_TO+1); TIT_TO=T3_TO+1;
    end
    mdot_core_TO = mdot_TO/(1+BPR);
    mf_TO = mdot_core_TO * cp_h * (TIT_TO - T3_TO) / (eta_comb * LHV);
    if mf_TO < 0; pen=pen+BIG; mf_TO=0; end
    mdot4_TO = mdot_core_TO + mf_TO;

    % --- HPT (HPC'yi besler) ---
    W_HPC = mdot_core_TO * cp_c * (T3_TO - T25_TO);  % HPC güç gereksinimi
    % HPT entalpisi: W_HPT = eta_mHP * mdot4 * cp_h * (T4 - T45)
    dT_HPT = W_HPC / (eta_mHP * mdot4_TO * cp_h);
    T45_TO  = TIT_TO - dT_HPT;
    if T45_TO < 800; pen=pen+BIG*(800-T45_TO); T45_TO=800; end
    % HPT basınç oranı (türbin genleşme):
    PR_HPT = (T45_TO/TIT_TO)^(n_h/eta_HPT_p);
    P45_TO  = P4_TO * PR_HPT;

    % --- IPT (LPC'yi besler) ---
    W_LPC = mdot_core_TO * cp_c * (T25_TO - T21_TO);
    dT_IPT = W_LPC / (eta_mIP * mdot4_TO * cp_h);
    T49_TO  = T45_TO - dT_IPT;
    if T49_TO < 700; pen=pen+BIG*(700-T49_TO); T49_TO=700; end
    PR_IPT = (T49_TO/T45_TO)^(n_h/eta_IPT_p);
    P49_TO  = P45_TO * PR_IPT;

    % --- LPT (Fan'ı besler) ---
    W_fan = (mdot_TO) * cp_c * (T13_TO - T02_TO);  % Tüm akış fan'dan geçer
    dT_LPT = W_fan / (eta_mLP * mdot4_TO * cp_h);
    T5_TO   = T49_TO - dT_LPT;
    if T5_TO < 500; pen=pen+BIG*(500-T5_TO); T5_TO=500; end
    PR_LPT = (T5_TO/T49_TO)^(n_h/eta_LPT_p);
    P5_TO   = P49_TO * PR_LPT;

    % --- Core lüle (İstasyon 8) ---
    P8_TO  = P5_TO * PR_cor_noz;
    if P8_TO < P0_TO; pen=pen+BIG*(P0_TO-P8_TO)/P0_TO*100; end
    P8e_TO = min(P8_TO, P0_TO);
    T8e_TO = T5_TO * (P8e_TO/P8_TO)^((gamma_h-1)/gamma_h);
    V8_TO  = sqrt(max(0, 2*cp_h*(T5_TO-T8e_TO)));

    % --- Bypass lüle (İstasyon 18) ---
    mdot_bypass_TO = mdot_TO * BPR/(1+BPR);
    P18_TO  = P13_TO * PR_byp_noz;
    if P18_TO < P0_TO; pen=pen+BIG*(P0_TO-P18_TO)/P0_TO*100; end
    P18e_TO = min(P18_TO, P0_TO);
    T18e_TO = T13_TO * (P18e_TO/P18_TO)^((gamma_c-1)/gamma_c);
    V18_TO  = sqrt(max(0, 2*cp_c*(T13_TO-T18e_TO)));

    % --- Toplam kalkış itkisi ---
    thrust_TO = (mdot_core_TO+mf_TO)*V8_TO + mdot_bypass_TO*V18_TO;  % V0=0
    if thrust_TO < 350e3; pen=pen+BIG*(350e3-thrust_TO)/1e3; end

    %% ====== SEYİR MODU ======
    % Seyirde hava debisi: mdot_cr hesabı
    % Sabit giriş alanı yaklaşımı: A_inlet = mdot_TO/rho_TO/V_eff
    % Kalkışta V=0 -> tanımsız. Pratik kabul:
    % mdot_cr = mdot_TO * (rho0_cr * V_cr) / (rho0_TO * a_TO * 0.3)
    % [kalkışta ortalama giriş hızı ~0.3*a_TO varsayımı]
    a_TO = sqrt(gamma_c * 287 * T02_TO);
    rho0_TO = P02_TO / (287 * T02_TO);
    rho0_cr = P02_cr / (287 * T02_cr);
    V_inlet_TO = 0.3 * a_TO;  % ~100 m/s kalkış giriş hızı
    A_inlet = mdot_TO / (rho0_TO * V_inlet_TO);
    mdot_cr = rho0_cr * V_cr * A_inlet;

    T21_cr = T02_cr * PR_fan^((gamma_c-1)/(gamma_c*eta_fan_p));
    P21_cr = P02_cr * PR_fan;
    T25_cr = T21_cr * PR_LPC^((gamma_c-1)/(gamma_c*eta_LPC_p));
    P25_cr = P21_cr * PR_LPC;
    T3_cr  = T25_cr * PR_HPC^((gamma_c-1)/(gamma_c*eta_HPC_p));
    P3_cr  = P25_cr * PR_HPC;
    P4_cr  = P3_cr * PR_comb;

    if TIT_cr < T3_cr; pen=pen+BIG*(T3_cr-TIT_cr+1); TIT_cr=T3_cr+1; end
    if TIT_cr > 1920;  pen=pen+BIG*(TIT_cr-1920);     TIT_cr=1920;    end
    mdot_core_cr = mdot_cr/(1+BPR);
    mf_cr = mdot_core_cr * cp_h * (TIT_cr - T3_cr) / (eta_comb * LHV);
    if mf_cr<0; pen=pen+BIG; mf_cr=0; end
    mdot4_cr = mdot_core_cr + mf_cr;

    W_HPC_cr = mdot_core_cr*cp_c*(T3_cr-T25_cr);
    dT_HPT_cr = W_HPC_cr/(eta_mHP*mdot4_cr*cp_h);
    T45_cr = TIT_cr - dT_HPT_cr;
    if T45_cr<600; pen=pen+BIG*(600-T45_cr); T45_cr=600; end
    PR_HPT_cr = (T45_cr/TIT_cr)^(n_h/eta_HPT_p);
    P45_cr = P4_cr * PR_HPT_cr;

    W_LPC_cr = mdot_core_cr*cp_c*(T25_cr-T21_cr);
    dT_IPT_cr = W_LPC_cr/(eta_mIP*mdot4_cr*cp_h);
    T49_cr = T45_cr - dT_IPT_cr;
    if T49_cr<500; pen=pen+BIG*(500-T49_cr); T49_cr=500; end
    PR_IPT_cr = (T49_cr/T45_cr)^(n_h/eta_IPT_p);
    P49_cr = P45_cr * PR_IPT_cr;

    W_fan_cr = mdot_cr*cp_c*(T21_cr-T02_cr);
    dT_LPT_cr = W_fan_cr/(eta_mLP*mdot4_cr*cp_h);
    T5_cr = T49_cr - dT_LPT_cr;
    if T5_cr<400; pen=pen+BIG*(400-T5_cr); T5_cr=400; end
    PR_LPT_cr = (T5_cr/T49_cr)^(n_h/eta_LPT_p);
    P5_cr = P49_cr * PR_LPT_cr;

    P8_cr  = P5_cr*PR_cor_noz;
    if P8_cr<P0_cr; pen=pen+BIG*(P0_cr-P8_cr)/P0_cr*100; end
    P8e_cr = min(P8_cr,P0_cr);
    T8e_cr = T5_cr*(P8e_cr/P8_cr)^((gamma_h-1)/gamma_h);
    V8_cr  = sqrt(max(0,2*cp_h*(T5_cr-T8e_cr)));

    mdot_bypass_cr = mdot_cr*BPR/(1+BPR);
    P18_cr  = P21_cr*PR_byp_noz;  % P21 seyir
    if P18_cr<P0_cr; pen=pen+BIG*(P0_cr-P18_cr)/P0_cr*100; end
    P18e_cr = min(P18_cr,P0_cr);
    T18e_cr = T21_cr*(P18e_cr/P18_cr)^((gamma_c-1)/gamma_c);
    V18_cr  = sqrt(max(0,2*cp_c*(T21_cr-T18e_cr)));

    thrust_CR = (mdot_core_cr+mf_cr)*V8_cr + mdot_bypass_cr*V18_cr - mdot_cr*V_cr;
    if thrust_CR<90e3; pen=pen+BIG*(90e3-thrust_CR)/1e3; end

    mf_tot_cr = mf_cr;
    if thrust_CR>0; sfc=mf_tot_cr/thrust_CR; else; sfc=inf; pen=pen+BIG; end

    %% Ağırlık modeli (basit parametrik)
    % Brayton verisi: W[kg] = 0.084 * (Thrust_TO[N])^0.75 / g^0.75 * OPR_factor
    % Daha yaygın: W ~ thrust_TO / (g * TW_ref) için TW_ref = 5-7
    % Kullanacağımız: Jenkinson 2-spool turbofan modeli
    % W_eng [kg] = 14.7 * (mdot_core_TO)^0.85 * (OPR)^0.3 * (1+BPR)^0.5
    W_eng = 14.7 * (mdot_core_TO)^0.85 * OPR^0.3 * (1+BPR)^0.5;
    if W_eng > 6000; pen=pen+BIG*(W_eng-6000)/1000; end

    % İtki/ağırlık oranı
    TW = thrust_TO / (W_eng * 9.81);
    if TW < 5; pen=pen+BIG*(5-TW); end

    % Lüle alanları (sürekli akış: A = mdot/(rho*V))
    % Core lüle: A8 = (mdot_core+mf)/rho8/V8
    rho8 = P8e_TO/(287*T8e_TO);
    A8   = (mdot_core_TO+mf_TO)/max(rho8*V8_TO,1e-6);
    % Bypass lüle: A18
    rho18 = P18e_TO/(287*T18e_TO);
    A18   = mdot_bypass_TO/max(rho18*V18_TO,1e-6);
    A_total = A8 + A18;
    if A_total > 6; pen=pen+BIG*(A_total-6); end

    pen_out = pen;
    f = sfc + pen;
end
function s=ternary(c,a,b); if c; s=a; else; s=b; end; end