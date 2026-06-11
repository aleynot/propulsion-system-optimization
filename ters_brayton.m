%% =========================================================
%  İTKİ SİSTEMLERİ PERFORMANS OPTİMİZASYONU - PROJE 1
%  Ters Brayton Çevrimli Motor PSO Optimizasyonu
%
%  MODELİN FİZİKSEL TEMELİ:
%  --------------------------------------------------------
%  Kalkış (M=0, SL, mdot=100 kg/s):
%    Preburner: sabit hacim yanma modeli -> P03 = P02*(T03/T02)
%    [İdeal motor kabulü: yanma odası basınç artışını modeller]
%    Türbin: P04 = PRt * P03 (genleşme, PRt < 1)
%    Soğutma: T05 = T04 - dT, P05 = P04
%    Afterburner: T07 > T05, P07 = P05
%    Lüle: optimum genleşme P9 = P_amb
%
%  Seyir (M=5.5, 25km, mdot=100 kg/s sabit):
%    RAM basıncı: P02_cr = 23.4 atm >> P0_cr
%    Preburner: sabit basınç (P03=P02), T03_cr > T02_cr=1562 K
%    [Not: T03_cr kısıtı <1600 K -> sadece 38 K ısıtma imkanı]
%    Türbin, soğutma, afterburner, lüle aynı şema
%
%  TASARIM PARAMETRELERİ (8 adet, mdot sabit=100 kg/s):
%    x(1) T03_TO   preburner çıkış sıc - kalkış [K]   [T02+50, 1800]
%    x(2) T07_TO   afterburner çıkış sıc - kalkış [K] [800, 2000]
%    x(3) dTcl_TO  soğutma ΔT - kalkış [K]            [10, 500]
%    x(4) PRt_TO   türbin P04/P03 - kalkış [-]         [0.10, 0.99]
%    x(5) T03_cr   preburner çıkış sıc - seyir [K]    [T02_cr, 1600]
%    x(6) T07_cr   afterburner çıkış sıc - seyir [K]  [900, 1700]
%    x(7) dTcl_cr  soğutma ΔT - seyir [K]             [10, 300]
%    x(8) PRt_cr   türbin P04/P03 - seyir [-]          [0.10, 0.99]
%
%  NOT: mdot=100 kg/s ile seyirde 300 kN kısıtı matematiksel
%  olarak aşılması güç bir sınırdır (V9_max~1712 m/s, V_cr=1641 m/s,
%  net momentum farkı ~71 m/s). PSO fiziksel optimumu bulur.
%% =========================================================
clear; clc; close all;

gamma=1.4; R=287; cp=1005; LHV=42.8e6;
T0_TO=288.15; P0_TO=101325;
T02_TO=T0_TO; P02_TO=P0_TO;
mdot=100;  % Her iki modda sabit hava debisi

T0_cr=221.55; P0_cr=2549.2; M_cr=5.5;
T02_cr=T0_cr*(1+(gamma-1)/2*M_cr^2);
P02_cr=P0_cr*(T02_cr/T0_cr)^(gamma/(gamma-1));
a_cr=sqrt(gamma*R*T0_cr); V_cr=M_cr*a_cr;

fprintf('=== GİRİŞ KOŞULLARI ===\n');
fprintf('Kalkış: T02=%.2f K, P02=%.1f Pa\n',T02_TO,P02_TO);
fprintf('Seyir : T02=%.2f K, P02=%.1f Pa\n',T02_cr,P02_cr);
fprintf('V_cr  = %.1f m/s,  mdot = %.0f kg/s\n\n',V_cr,mdot);

%% PSO AYARLARI
N_par=150; N_iter=600; w_max=0.9; w_min=0.4; c1=1.6; c2=1.6; n_var=8;

lb=[T02_TO+50, 800,  10, 0.10, T02_cr,  900,  10, 0.10];
ub=[1800,      2000, 500, 0.99, 1600,   1700, 300, 0.99];

rng(42);
X=lb+rand(N_par,n_var).*(ub-lb);
V=0.05*(rand(N_par,n_var)-0.5).*(ub-lb);
pBest=X; pBest_f=inf(N_par,1);
for i=1:N_par
    pBest_f(i)=engine(X(i,:),T02_TO,P02_TO,T02_cr,P02_cr,...
                      mdot,gamma,cp,LHV,V_cr,P0_TO,P0_cr);
end
[gBest_f,idx]=min(pBest_f); gBest=pBest(idx,:);
hist_f=nan(N_iter,1); hist_x=nan(N_iter,n_var);

fprintf('%-8s %-16s %-14s %-14s\n','İter','gBest_f','İtki_TO[kN]','İtki_CR[kN]');
fprintf('%s\n',repmat('-',56,1));

for iter=1:N_iter
    w=w_max-(w_max-w_min)*iter/N_iter;
    r1=rand(N_par,n_var); r2=rand(N_par,n_var);
    V=w*V+c1*r1.*(pBest-X)+c2*r2.*(repmat(gBest,N_par,1)-X);
    Vmax=0.25*(ub-lb); V=max(min(V,Vmax),-Vmax);
    X=X+V;
    for d=1:n_var
        lo=X(:,d)<lb(d); hi=X(:,d)>ub(d);
        X(lo,d)=lb(d)+rand(sum(lo),1).*(ub(d)-lb(d))*0.05;
        X(hi,d)=ub(d)-rand(sum(hi),1).*(ub(d)-lb(d))*0.05;
        V(lo,d)=abs(V(lo,d)); V(hi,d)=-abs(V(hi,d));
    end
    for i=1:N_par
        fi=engine(X(i,:),T02_TO,P02_TO,T02_cr,P02_cr,...
                  mdot,gamma,cp,LHV,V_cr,P0_TO,P0_cr);
        if fi<pBest_f(i); pBest_f(i)=fi; pBest(i,:)=X(i,:); end
        if fi<gBest_f;    gBest_f=fi;    gBest=X(i,:);       end
    end
    hist_f(iter)=gBest_f; hist_x(iter,:)=gBest;
    if mod(iter,50)==0||iter==1
        [~,tTO,tCR]=engine(gBest,T02_TO,P02_TO,T02_cr,P02_cr,...
                           mdot,gamma,cp,LHV,V_cr,P0_TO,P0_cr);
        fprintf('%-8d %-16.6g %-14.2f %-14.2f\n',iter,gBest_f,tTO/1e3,tCR/1e3);
    end
end

[~,thr_TO,thr_CR,sfc_opt,T05_TO]=engine(gBest,T02_TO,P02_TO,T02_cr,P02_cr,...
    mdot,gamma,cp,LHV,V_cr,P0_TO,P0_cr);

fprintf('\n=== OPTİMUM TASARIM PARAMETRELERİ ===\n');
nms={'T03_TO  [K]','T07_TO  [K]','dTcl_TO [K]','PRt_TO  [-]',...
     'T03_cr  [K]','T07_cr  [K]','dTcl_cr [K]','PRt_cr  [-]'};
for k=1:n_var; fprintf('  %-18s = %10.4f\n',nms{k},gBest(k)); end
fprintf('\n=== PERFORMANS ===\n');
fprintf('  SFC (seyir)         = %.4e kg/(N.s)\n',sfc_opt);
fprintf('  İtki (kalkış)       = %.2f kN\n',thr_TO/1e3);
fprintf('  İtki (seyir)        = %.2f kN\n',thr_CR/1e3);
fprintf('  T05 soğutma - TO    = %.2f K\n',T05_TO);
fprintf('\n=== KISIT KONTROL ===\n');
ok=@(c)ternary(c,'✓ SAĞLANDI','✗ SAĞLANMADI');
fprintf('  İtki_TO > 120 kN    : %s (%.1f kN)\n',ok(thr_TO>120e3),thr_TO/1e3);
fprintf('  İtki_CR > 300 kN    : %s (%.1f kN)\n',ok(thr_CR>300e3),thr_CR/1e3);
fprintf('  T05_TO  > 273 K     : %s (%.1f K)\n', ok(T05_TO>273),T05_TO);
fprintf('  T03_TO  < 1800 K    : %s\n',ok(gBest(1)<1800));
fprintf('  T07_TO  < 2000 K    : %s\n',ok(gBest(2)<2000));
fprintf('  dTcl_TO < 500 K     : %s\n',ok(gBest(3)<500));
fprintf('  PRt_TO  > 0.1       : %s\n',ok(gBest(4)>0.1));
fprintf('  T03_cr  < 1600 K    : %s\n',ok(gBest(5)<1600));
fprintf('  T07_cr  < 1700 K    : %s\n',ok(gBest(6)<1700));
fprintf('  dTcl_cr < 300 K     : %s\n',ok(gBest(7)<300));

%% GRAFİKLER
figure('Name','PSO Yakınsama','Position',[50 50 1400 900]);
subplot(3,3,1); plot(1:N_iter,hist_f,'b-','LineWidth',1.8);
xlabel('İterasyon'); ylabel('Amaç Fonk.'); title('Global En İyi SFC'); grid on;
lbl={'T_{03,TO} [K]','T_{07,TO} [K]','\DeltaT_{cl,TO} [K]','PR_{t,TO} [-]',...
     'T_{03,cr} [K]','T_{07,cr} [K]','\DeltaT_{cl,cr} [K]','PR_{t,cr} [-]'};
clr=lines(8);
for k=1:8
    subplot(3,3,k+1);
    plot(1:N_iter,hist_x(:,k),'Color',clr(k,:),'LineWidth',1.5);
    xlabel('İterasyon'); ylabel(lbl{k}); title(lbl{k}); grid on;
end
sgtitle('PSO Yakınsama Grafikleri','FontSize',12,'FontWeight','bold');

figure('Name','Performans','Position',[100 100 900 420]);
subplot(1,2,1);
b=bar([thr_TO/1e3,thr_CR/1e3],'FaceColor','flat');
b.CData=[0.2 0.4 0.8;0.1 0.6 0.3];
set(gca,'XTickLabel',{'Kalkış','Hipersonik Seyir'});
ylabel('İtki [kN]'); title('Motor İtkisi'); hold on;
yline(120,'r--','LineWidth',1.5,'Label','120 kN');
yline(300,'g--','LineWidth',1.5,'Label','300 kN'); grid on;
subplot(1,2,2);
v=~isnan(hist_f)&hist_f<1e9;
if any(v); semilogy(find(v),hist_f(v),'m-','LineWidth',2); end
xlabel('İterasyon'); ylabel('SFC (log)'); title('SFC Yakınsaması'); grid on;
sgtitle('Motor Performans Özeti','FontSize',12,'FontWeight','bold');
fprintf('\nBitti.\n');

%% =========================================================
function [f, thrust_TO, thrust_CR, sfc, T05_TO] = ...
    engine(x, T02_TO, P02_TO, T02_cr, P02_cr, ...
           mdot, gamma, cp, LHV, V_cr, P0_TO, P0_cr)

    f=1e10; thrust_TO=0; thrust_CR=0; sfc=inf; T05_TO=0;
    T03_TO=x(1); T07_TO=x(2); dTcl_TO=x(3); PRt_TO=x(4);
    T03_cr=x(5); T07_cr=x(6); dTcl_cr=x(7); PRt_cr=x(8);
    pen=0; BIG=1e6;

    %% KALKI$
    if T03_TO<=T02_TO; pen=pen+BIG*(T02_TO-T03_TO+1); T03_TO=T02_TO+1; end
    P03_TO=P02_TO*(T03_TO/T02_TO);   % sabit hacim yanma
    mf_pre_TO=mdot*cp*(T03_TO-T02_TO)/LHV;

    P04_TO=PRt_TO*P03_TO;
    T04_TO=T03_TO*(P04_TO/P03_TO)^((gamma-1)/gamma);

    T05_TO=T04_TO-dTcl_TO; P05_TO=P04_TO;
    if T05_TO<273; pen=pen+BIG*(274-T05_TO); T05_TO=274; end

    mdot3_TO=mdot+mf_pre_TO;
    if T07_TO<=T05_TO; pen=pen+BIG; T07_TO=T05_TO+1; end
    mf_ab_TO=mdot3_TO*cp*(T07_TO-T05_TO)/LHV;
    if mf_ab_TO<0; pen=pen+BIG; mf_ab_TO=0; end
    P07_TO=P05_TO;

    if P07_TO<P0_TO; pen=pen+BIG*(P0_TO-P07_TO)/P0_TO*200; end
    T09_TO=T07_TO*(min(P07_TO,P0_TO)/P07_TO)^((gamma-1)/gamma);
    V9_TO=sqrt(max(0,2*cp*(T07_TO-T09_TO)));
    thrust_TO=(mdot+mf_pre_TO+mf_ab_TO)*V9_TO;
    if thrust_TO<120e3; pen=pen+BIG*(120e3-thrust_TO)/1e3; end

    %% SEYİR
    % T03_cr >= T02_cr = 1562 K (ram sıkıştırma sonucu sıcak gaz)
    if T03_cr<=T02_cr; T03_cr=T02_cr; end  % minimum T02_cr
    P03_cr=P02_cr;
    mf_pre_cr=mdot*cp*(T03_cr-T02_cr)/LHV;
    if mf_pre_cr<0; mf_pre_cr=0; end

    P04_cr=PRt_cr*P03_cr;
    T04_cr=T03_cr*(P04_cr/P03_cr)^((gamma-1)/gamma);

    T05_cr=T04_cr-dTcl_cr; P05_cr=P04_cr;

    mdot3_cr=mdot+mf_pre_cr;
    if T07_cr<=T05_cr; pen=pen+BIG; T07_cr=T05_cr+1; end
    mf_ab_cr=mdot3_cr*cp*(T07_cr-T05_cr)/LHV;
    if mf_ab_cr<0; pen=pen+BIG; mf_ab_cr=0; end
    P07_cr=P05_cr;

    if P07_cr<P0_cr; pen=pen+BIG*(P0_cr-P07_cr)/P0_cr*200; end
    T09_cr=T07_cr*(min(P07_cr,P0_cr)/P07_cr)^((gamma-1)/gamma);
    V9_cr=sqrt(max(0,2*cp*(T07_cr-T09_cr)));
    thrust_CR=(mdot+mf_pre_cr+mf_ab_cr)*V9_cr - mdot*V_cr;
    if thrust_CR<300e3; pen=pen+BIG*(300e3-thrust_CR)/1e3; end

    mf_tot_cr=mf_pre_cr+mf_ab_cr;
    if thrust_CR>0; sfc=mf_tot_cr/thrust_CR; else; sfc=inf; pen=pen+BIG; end
    f=sfc+pen;
end
function s=ternary(c,a,b); if c; s=a; else; s=b; end; end