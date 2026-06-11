% q3_output_feedback_fast.m
clear; clc; close all;

c1 = 2; c2 = 3;

a    = 20;      % dirty-derivative speed
khat = 20;      % z1hat filter gain
umax = 10;      % control saturation
satU = @(u) max(-umax, min(umax, u));

alpha = @(z1) (-(1+c1)*z1 - z1.^3);
dalph = @(z1) (-(1+c1) - 3*z1.^2);
f2    = @(z1) (z1 - 3*z1.^3);

% state: [z1; y; zeta; z1hat]
z10 = 0.8; y0 = -0.5;
s0 = [z10; y0; y0; 0.0];     % zeta0=y0

tspan = [0 10];

opts = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',1e-2, ...
              'OutputFcn', @odeprog);

dyn = @(t,s) cl_fast(t,s,alpha,dalph,f2,a,khat,satU);

[t,s] = ode23t(dyn, tspan, s0, opts);   % ode23t genelde daha hızlı

z1 = s(:,1); y = s(:,2); z1h = s(:,4);

figure; grid on; hold on;
plot(t,z1,'LineWidth',1.5);
plot(t,y ,'LineWidth',1.5);
legend({'z1=x1','z2=y'},'Location','best');
xlabel('t'); title('Durumlar');

figure; grid on;
plot(t,z1-z1h,'LineWidth',1.5);
xlabel('t'); title('z1 - z1hat');

% -----------------------------
function ds = cl_fast(~,s,alpha,dalph,f2,a,khat,satU)
    z1   = s(1);
    y    = s(2);
    zeta = s(3);
    z1h  = s(4);

    % dirty derivative
    dzeta = -a*zeta + a*y;
    ydot_hat = a*(y - zeta);

    % control guess using current z1hat
    e2 = y - alpha(z1h);
    u  = -f2(z1h) - dalph(z1h)*2*z1h - z1h + (dalph(z1h)-3)*e2; % c1=2,c2=3 inline
    u  = satU(u);

    % algebraic inversion: solve 3 z^3 - z + (ydot_hat - u) = 0
    w = ydot_hat - u;
    z1_alg = cubic_inverse_pick(w, z1h);

    % filter estimate
    dz1h = khat*(z1_alg - z1h);

    % plant
    dz1 = z1 + z1^3 + y;
    dy  = z1 - 3*z1^3 + u;

    ds = [dz1; dy; dzeta; dz1h];
end

function z = cubic_inverse_pick(w, zref)
    % 3 z^3 - z + w = 0  -> coefficients [3 0 -1 w]
    r = roots([3 0 -1 w]);          % 3 roots (possibly complex)
    r = r(abs(imag(r))<1e-8);       % keep real roots
    if isempty(r)
        z = zref;                   % fallback
        return;
    end
    [~,ix] = min(abs(r - zref));    % pick closest to current estimate
    z = real(r(ix));
    z = max(-10, min(10, z));       % safety clamp
end

function status = odeprog(t,~,flag)
    persistent t0
    status = 0;
    if isempty(flag)
        if isempty(t0), t0 = t(1); end
        if t(end) - t0 > 0.5
            fprintf('t = %.3f\n', t(end));
            t0 = t(end);
        end
    elseif strcmp(flag,'init')
        fprintf('Integrating...\n');
        t0 = [];
    elseif strcmp(flag,'done')
        fprintf('Done.\n');
        t0 = [];
    end
end
