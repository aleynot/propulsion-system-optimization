%% THREE-SHAFT TURBOFAN ENGINE OPTIMIZATION USING PSO
% Propulsion System Performance Optimization - Homework 3
% Objective : Minimize Specific Fuel Consumption (SFC)
% Constraint: Net Thrust >= 100 kN
%
% Design Variables:
%   x(1) = LPC total pressure ratio (pi_LPC) : [2, 10]
%   x(2) = HPC total pressure ratio (pi_HPC) : [2, 10]
%   x(3) = Fan total pressure ratio (pi_fan) : [1.1, 4]
%   x(4) = Turbine inlet total temperature T4 : [1600, 2000] K
%   x(5) = Bypass ratio (BPR)                 : [1, 20]
%
% Engine Architecture (3-Shaft):
%   LP shaft : Fan  <-- driven by LPT
%   IP shaft : LPC  <-- driven by IPT (MPT)
%   HP shaft : HPC  <-- driven by HPT
%
% Station numbering:
%   2   : Engine face (ambient inlet)
%   021 : After intake diffuser
%   13  : After fan (bypass stream + core feed)
%   25  : After LPC
%   3   : After HPC
%   4   : Combustor exit / HPT inlet
%   45  : After HPT
%   5   : After IPT
%   6   : After LPT  (core nozzle inlet)
%   8   : Core nozzle exit
%   18  : Bypass nozzle exit

clc; clear; close all;

fprintf('===========================================\n');
fprintf('  THREE-SHAFT TURBOFAN PSO OPTIMIZATION   \n');
fprintf('===========================================\n\n');

%% =========================================================
%  1. CONSTANTS & FLIGHT CONDITIONS  (ISA @ 10 km)
% ==========================================================
h       = 10000;          % [m]
M_inf   = 0.85;

% ISA @ 10 km
T_s     = 223.15;         % [K]  static temperature
P_s     = 26499.87;       % [Pa] static pressure

% Gas constants – cold / hot
Cp_c    = 1005;           % [J/kg/K]
Cp_h    = 1150;           % [J/kg/K]
R_gas   = 287;            % [J/kg/K]  (same both zones per Table 4)
gam_c   = Cp_c/(Cp_c - R_gas);   % ≈ 1.400
gam_h   = Cp_h/(Cp_h - R_gas);   % ≈ 1.333

% Total conditions at engine face (station 2)
T02 = T_s * (1 + (gam_c-1)/2 * M_inf^2);           % [K]
P02 = P_s * (1 + (gam_c-1)/2 * M_inf^2)^(gam_c/(gam_c-1));  % [Pa]
a_s = sqrt(gam_c * R_gas * T_s);
V0  = M_inf * a_s;                                  % [m/s] flight speed

fprintf('--- Atmospheric / Inlet Conditions ---\n');
fprintf('  T_static  = %.2f K   P_static  = %.2f Pa\n', T_s, P_s);
fprintf('  T02       = %.2f K   P02       = %.2f Pa\n', T02, P02);
fprintf('  V0 (flight)= %.2f m/s\n\n', V0);

%% =========================================================
%  2. ENGINE PARAMETERS
% ==========================================================
mdot_tot = 400;           % [kg/s] total intake mass flow
LHV      = 42.8e6;        % [J/kg] fuel lower heating value

%% =========================================================
%  3. COMPONENT EFFICIENCIES & PRESSURE RATIOS  (Table 5)
% ==========================================================
pi_diff   = 0.97;   % Intake total-pressure recovery
eta_comb  = 0.985;  % Combustor efficiency
pi_comb   = 0.97;   % Combustor pressure ratio
e_fan     = 0.88;   % Fan polytropic efficiency
e_LPC     = 0.89;   % LPC polytropic efficiency
e_HPC     = 0.90;   % HPC polytropic efficiency
e_LPT     = 0.90;   % LPT polytropic efficiency
e_IPT     = 0.91;   % IPT polytropic efficiency
e_HPT     = 0.92;   % HPT polytropic efficiency
eta_mHP   = 0.99;   % HP shaft mechanical efficiency
eta_mIP   = 0.995;  % IP shaft mechanical efficiency
eta_mLP   = 0.995;  % LP shaft mechanical efficiency
pi_byp_n  = 0.98;   % Bypass nozzle pressure ratio
pi_cor_n  = 0.98;   % Core nozzle pressure ratio

%% =========================================================
%  4. PSO HYPER-PARAMETERS
% ==========================================================
Np    = 80;     % swarm size
Nit   = 400;    % max iterations
w0    = 0.90;   % initial inertia
wf    = 0.40;   % final inertia
c1    = 2.0;    % cognitive
c2    = 2.0;    % social
lam   = 5e8;    % penalty multiplier  (Pa-scale: thrust in N)

% Bounds  [pi_LPC  pi_HPC  pi_fan  T4(K)  BPR]
lb = [2,    2,    1.1,  1600,  1 ];
ub = [10,   10,   4.0,  2000,  20];
nv = 5;

%% =========================================================
%  5. PACK PARAMETERS INTO STRUCT FOR CLEAN PASSING
% ==========================================================
eng.T02=T02; eng.P02=P02; eng.V0=V0; eng.P_s=P_s;
eng.mdot_tot=mdot_tot; eng.LHV=LHV;
eng.Cp_c=Cp_c; eng.Cp_h=Cp_h; eng.R=R_gas;
eng.gam_c=gam_c; eng.gam_h=gam_h;
eng.pi_diff=pi_diff; eng.eta_comb=eta_comb; eng.pi_comb=pi_comb;
eng.e_fan=e_fan; eng.e_LPC=e_LPC; eng.e_HPC=e_HPC;
eng.e_LPT=e_LPT; eng.e_IPT=e_IPT; eng.e_HPT=e_HPT;
eng.eta_mHP=eta_mHP; eng.eta_mIP=eta_mIP; eng.eta_mLP=eta_mLP;
eng.pi_byp_n=pi_byp_n; eng.pi_cor_n=pi_cor_n;

%% =========================================================
%  6. INITIALIZE SWARM
% ==========================================================
rng(7);
X  = lb + rand(Np,nv).*(ub-lb);
Vp = zeros(Np,nv);
Vmax = 0.25*(ub-lb);

Fval = zeros(Np,1);
for i=1:Np
    Fval(i) = penObj(X(i,:), eng, lam);
end

pB    = X;
pBval = Fval;
[gBval, gi] = min(Fval);
gB = X(gi,:);

hist_SFC    = nan(Nit,1);
hist_F      = nan(Nit,1);
hist_gBval  = nan(Nit,1);

%% =========================================================
%  7. PSO MAIN LOOP
% ==========================================================
fprintf('--- PSO Running  (Np=%d, Nit=%d) ---\n\n', Np, Nit);
for it = 1:Nit
    w = w0 - (w0-wf)*it/Nit;

    r1 = rand(Np,nv);
    r2 = rand(Np,nv);
    Vp = w*Vp + c1*r1.*(pB-X) + c2*r2.*(gB-X);
    Vp = max(-Vmax, min(Vmax, Vp));
    X  = X + Vp;

    % Clamp to bounds (absorbing)
    X = max(lb, min(ub, X));

    for i=1:Np
        f = penObj(X(i,:), eng, lam);
        if f < pBval(i)
            pBval(i) = f;
            pB(i,:)  = X(i,:);
        end
    end
    [cur,ci] = min(pBval);
    if cur < gBval
        gBval = cur;
        gB    = pB(ci,:);
    end

    % Evaluate true performance at global best
    [sfc_it, F_it] = truePerfQuick(gB, eng);
    hist_SFC(it)   = sfc_it*1e6;   % g/(kN·s)
    hist_F(it)     = F_it/1e3;     % kN
    hist_gBval(it) = gBval;

    if mod(it,50)==0 || it==1
        fprintf('  Iter %4d | SFC=%7.4f g/(kN·s) | Thrust=%7.2f kN\n', ...
            it, sfc_it*1e6, F_it/1e3);
    end
end

%% =========================================================
%  8. FINAL DETAILED EVALUATION
% ==========================================================
[sfc_opt, F_opt, perf] = fullEngine(gB, eng);

fprintf('\n===========================================\n');
fprintf('        OPTIMIZATION RESULTS               \n');
fprintf('===========================================\n');
fprintf('\n--- Optimal Design Variables ---\n');
fprintf('  pi_LPC (LPC pressure ratio)   : %8.4f\n', gB(1));
fprintf('  pi_HPC (HPC pressure ratio)   : %8.4f\n', gB(2));
fprintf('  pi_fan (Fan pressure ratio)   : %8.4f\n', gB(3));
fprintf('  T4     (Turbine inlet temp)   : %8.2f K\n', gB(4));
fprintf('  BPR    (Bypass ratio)         : %8.4f\n', gB(5));

fprintf('\n--- Station Total Temperatures ---\n');
fn = fieldnames(perf.T);
for k=1:numel(fn)
    fprintf('  T%-6s = %9.2f K\n', fn{k}, perf.T.(fn{k}));
end

fprintf('\n--- Station Total Pressures ---\n');
fp = fieldnames(perf.P);
for k=1:numel(fp)
    fprintf('  P%-6s = %12.2f Pa\n', fp{k}, perf.P.(fp{k}));
end

fprintf('\n--- Engine Performance ---\n');
fprintf('  Overall Pressure Ratio (OPR)  : %8.2f\n',  perf.OPR);
fprintf('  mdot_core                     : %8.2f kg/s\n', perf.mdot_core);
fprintf('  mdot_bypass                   : %8.2f kg/s\n', perf.mdot_byp);
fprintf('  mdot_fuel                     : %8.4f kg/s\n', perf.mdot_f);
fprintf('  Fuel-Air Ratio (FAR)          : %8.5f\n',  perf.FAR);
fprintf('  V8  (core nozzle exit)        : %8.2f m/s\n', perf.V8);
fprintf('  V18 (bypass nozzle exit)      : %8.2f m/s\n', perf.V18);
fprintf('  Net Thrust                    : %8.2f kN\n',   F_opt/1e3);
fprintf('  Specific Thrust               : %8.4f N/(kg/s)\n', F_opt/mdot_tot);
fprintf('  SFC                           : %8.4f g/(kN·s)\n', sfc_opt*1e6);
fprintf('  Thermal Efficiency            : %8.2f %%\n',  perf.eta_th*100);
fprintf('  Propulsive Efficiency         : %8.2f %%\n',  perf.eta_pr*100);
fprintf('  Overall Efficiency            : %8.2f %%\n',  perf.eta_ov*100);

fprintf('\n--- Constraint ---\n');
if F_opt >= 100e3
    fprintf('  Thrust SATISFIED : %.2f kN >= 100 kN  [OK]\n', F_opt/1e3);
else
    fprintf('  Thrust VIOLATED  : %.2f kN < 100 kN   [FAIL]\n', F_opt/1e3);
end

%% =========================================================
%  9. PLOTS
% ==========================================================
figure('Name','PSO Convergence','Position',[60 60 1300 440]);

subplot(1,3,1);
plot(hist_SFC,'b','LineWidth',1.6);
xlabel('Iteration'); ylabel('SFC  [g/(kN·s)]');
title('SFC Convergence'); grid on;
xlim([1 Nit]);

subplot(1,3,2);
plot(hist_F,'r','LineWidth',1.6); hold on;
yline(100,'k--','100 kN','LabelVerticalAlignment','bottom','LineWidth',1.2);
xlabel('Iteration'); ylabel('Net Thrust  [kN]');
title('Thrust Convergence'); grid on;
xlim([1 Nit]);

subplot(1,3,3);
plot(hist_gBval,'g','LineWidth',1.6);
xlabel('Iteration'); ylabel('Penalized Objective  [kg/N·s]');
title('Global Best (Penalized)'); grid on;
xlim([1 Nit]);

sgtitle('PSO Optimization – Three-Shaft Turbofan Engine','FontSize',13,'FontWeight','bold');
saveas(gcf,'PSO_convergence.png');
fprintf('\nFigure saved: PSO_convergence.png\n');

% Station temperature plot
Tstations = [perf.T.s02, perf.T.s021, perf.T.s13, perf.T.s25, ...
             perf.T.s3,  perf.T.s4,   perf.T.s45, perf.T.s5,  perf.T.s6];
labels = {'02','021','13','25','3','4','45','5','6'};
idx_s  = 0:numel(Tstations)-1;

figure('Name','Station Temperatures','Position',[60 560 750 500]);
plot(idx_s, Tstations,'bo-','LineWidth',1.8,'MarkerFaceColor','b','MarkerSize',8);
for k=1:numel(labels)
    text(idx_s(k)+0.06, Tstations(k)+15, labels{k},'FontSize',9);
end
xlabel('Station Index'); ylabel('Total Temperature  [K]');
title('Engine Station Total Temperatures at Optimal Point');
grid on;
saveas(gcf,'Engine_station_temperatures.png');
fprintf('Figure saved: Engine_station_temperatures.png\n');

fprintf('\n--- Done ---\n');

%%==========================================================
%% LOCAL FUNCTIONS
%%==========================================================

%----------------------------------------------------------
% penObj : penalised objective (SFC + penalty for F<100kN)
%----------------------------------------------------------
function val = penObj(x, e, lam)
    [sfc, F] = truePerfQuick(x, e);
    if isnan(sfc) || sfc<=0
        val = 1e15; return;
    end
    viol = max(0, 100e3 - F);       % [N]  violation amount
    val  = sfc + lam * viol^2;
end

%----------------------------------------------------------
% truePerfQuick : fast SFC & thrust (no full struct output)
%----------------------------------------------------------
function [sfc, F] = truePerfQuick(x, e)
    [sfc, F, ~] = fullEngine(x, e);
end

%----------------------------------------------------------
% fullEngine : complete 1-D thermodynamic cycle
%----------------------------------------------------------
function [sfc, Fnet, p] = fullEngine(x, e)

sfc = NaN; Fnet = NaN; p = struct();

pi_LPC = x(1);  pi_HPC = x(2);  pi_fan = x(3);
T4     = x(4);  BPR    = x(5);

% Convenience aliases
Cp_c=e.Cp_c; Cp_h=e.Cp_h; gc=e.gam_c; gh=e.gam_h;

try
    %--- Station 2 (engine face) ---
    T02 = e.T02;  P02 = e.P02;

    %--- Station 021 : diffuser ---
    P021 = e.pi_diff * P02;
    T021 = T02;               % adiabatic

    %--- Mass flow split ---
    mdot_c  = e.mdot_tot / (1 + BPR);
    mdot_bp = BPR * mdot_c;

    %--- Station 13 : after fan (polytropic) ---
    T13 = T021 * pi_fan^( (gc-1)/(gc*e.e_fan) );
    P13 = P021 * pi_fan;

    if ~isreal(T13) || T13<=T021, return; end

    %--- Station 25 : after LPC ---
    T25 = T13 * pi_LPC^( (gc-1)/(gc*e.e_LPC) );
    P25 = P13 * pi_LPC;

    if ~isreal(T25) || T25<=T13, return; end

    %--- Station 3 : after HPC ---
    T3 = T25 * pi_HPC^( (gc-1)/(gc*e.e_HPC) );
    P3 = P25 * pi_HPC;

    if ~isreal(T3) || T3<=T25, return; end

    %--- Station 4 : combustor exit ---
    P4 = e.pi_comb * P3;
    % T4 = design variable

    % Fuel mass flow (energy balance on combustor)
    %  mdot_c * Cp_h * T4 = mdot_c * Cp_c * T3 + mdot_f * eta * LHV
    % => mdot_f = mdot_c*(Cp_h*T4 - Cp_c*T3) / (eta*LHV)
    mdot_f = mdot_c * (Cp_h*T4 - Cp_c*T3) / (e.eta_comb * e.LHV);

    if mdot_f <= 0, return; end

    FAR    = mdot_f / mdot_c;
    mdot_g = mdot_c + mdot_f;   % hot gas flow through turbines

    %--- Station 45 : after HPT  (drives HPC via HP shaft) ---
    W_HPC  = mdot_c * Cp_c * (T3 - T25);
    dT_HPT = W_HPC / (e.eta_mHP * mdot_g * Cp_h);
    T45    = T4 - dT_HPT;

    if T45 <= 0 || T45 >= T4, return; end

    % HPT pressure ratio (polytropic turbine expansion)
    % For turbine: T_exit/T_in = (P_exit/P_in)^[(gam-1)*e_poly/gam]
    % => pi_HPT = (T45/T4)^[gam/(gam-1)/e_poly]
    pi_HPT = (T45/T4)^( gh / ((gh-1)*e.e_HPT) );
    P45    = P4 * pi_HPT;

    if P45 <= 0 || pi_HPT >= 1, return; end

    %--- Station 5 : after IPT  (drives LPC via IP shaft) ---
    W_LPC  = mdot_c * Cp_c * (T25 - T13);
    dT_IPT = W_LPC / (e.eta_mIP * mdot_g * Cp_h);
    T5     = T45 - dT_IPT;

    if T5 <= 0 || T5 >= T45, return; end

    pi_IPT = (T5/T45)^( gh / ((gh-1)*e.e_IPT) );
    P5     = P45 * pi_IPT;

    if P5 <= 0 || pi_IPT >= 1, return; end

    %--- Station 6 : after LPT  (drives Fan via LP shaft) ---
    W_fan  = e.mdot_tot * Cp_c * (T13 - T021);
    dT_LPT = W_fan / (e.eta_mLP * mdot_g * Cp_h);
    T6     = T5 - dT_LPT;

    if T6 <= 50 || T6 >= T5, return; end    % physical lower bound 50 K

    pi_LPT = (T6/T5)^( gh / ((gh-1)*e.e_LPT) );
    P6     = P5 * pi_LPT;

    if P6 <= 0 || pi_LPT >= 1, return; end

    %--- Core nozzle (station 8) ---
    P6t  = e.pi_cor_n * P6;   % nozzle inlet total pressure
    T6t  = T6;                 % adiabatic nozzle

    [V8, P8e] = convergentNozzle(T6t, P6t, e.P_s, gh, e.R, Cp_h);
    if isnan(V8) || V8<=0, return; end

    %--- Bypass nozzle (station 18) ---
    P13t = e.pi_byp_n * P13;
    T13t = T13;

    [V18, P18e] = convergentNozzle(T13t, P13t, e.P_s, gc, e.R, Cp_c);
    if isnan(V18) || V18<=0, return; end

    %--- Nozzle exit areas (from continuity) ---
    rho8  = P8e  / (e.R * (T6t*(P8e/P6t)^((gh-1)/gh)));
    A8    = mdot_g / (rho8 * V8);

    rho18 = P18e / (e.R * (T13t*(P18e/P13t)^((gc-1)/gc)));
    A18   = mdot_bp / (rho18 * V18);

    %--- Net Thrust ---
    F_core = mdot_g   * V8  - mdot_c * e.V0 + (P8e  - e.P_s)*A8;
    F_byp  = mdot_bp  * V18 - mdot_bp* e.V0 + (P18e - e.P_s)*A18;
    Fnet   = F_core + F_byp;

    if Fnet <= 0 || isnan(Fnet), return; end

    sfc = mdot_f / Fnet;   % kg/(N·s)

    %--- Efficiencies ---
    Q_in  = mdot_f * e.LHV;
    %KE_out = 0.5*mdot_g*V8^2 + 0.5*mdot_bp*V18^2;

    KE_out = 0.5*mdot_g*(V8^2 - e.V0^2) + 0.5*mdot_bp*(V18^2 - e.V0^2);

    eta_pr = Fnet * e.V0 / KE_out;

    KE_in  = 0.5*e.mdot_tot*e.V0^2;
    eta_th = (KE_out - KE_in) / Q_in;
    %eta_pr = Fnet * e.V0 / max(KE_out - KE_in, 1e-6);
    eta_ov = Fnet * e.V0 / Q_in;

    %--- Pack struct ---
    p.T.s02=T02; p.T.s021=T021; p.T.s13=T13; p.T.s25=T25;
    p.T.s3=T3;   p.T.s4=T4;    p.T.s45=T45; p.T.s5=T5; p.T.s6=T6;

    p.P.s02=P02; p.P.s021=P021; p.P.s13=P13; p.P.s25=P25;
    p.P.s3=P3;   p.P.s4=P4;    p.P.s45=P45; p.P.s5=P5; p.P.s6=P6;

    p.OPR      = pi_fan * pi_LPC * pi_HPC;
    p.mdot_core= mdot_c;
    p.mdot_byp = mdot_bp;
    p.mdot_f   = mdot_f;
    p.FAR      = FAR;
    p.V8       = V8;
    p.V18      = V18;
    p.eta_th   = eta_th;
    p.eta_pr   = eta_pr;
    p.eta_ov   = eta_ov;

catch
    sfc = NaN; Fnet = NaN;
end
end

%----------------------------------------------------------
% convergentNozzle
%   Inputs : T_t, P_t  (total), P_amb, gamma, R, Cp
%   Output : V_exit, P_exit
%----------------------------------------------------------
function [V, Pe] = convergentNozzle(Tt, Pt, Pamb, gam, R, Cp)
    V = NaN; Pe = NaN;
    if Pt <= 0 || Tt <= 0, return; end

    PR_crit = (1 + (gam-1)/2)^(gam/(gam-1));   % ≈ 1.893 for gam=1.4

    if Pt/Pamb >= PR_crit
        % Choked
        Te = Tt * 2/(gam+1);
        V  = sqrt(gam * R * Te);
        Pe = Pt / PR_crit;
    else
        % Fully expanded to ambient
        Te = Tt * (Pamb/Pt)^((gam-1)/gam);
        dH = Cp * (Tt - Te);
        if dH <= 0, return; end
        V  = sqrt(2 * dH);
        Pe = Pamb;
    end
end