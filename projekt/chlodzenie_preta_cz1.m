%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CHŁODZENIE PRĘTA – CZĘŚĆ 1
%   Wersja zgodna ze szkieletem z instrukcji
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chlodzenie_preta
clc; clear; close all;

%% ================================
%   PARAMETRY MODELU FIZYCZNEGO
% ================================
params.h  = 160;      
params.A  = 0.0109;   
params.mb = 0.2;      
params.mw = 2.5;      
params.cb = 3.85;     
params.cw = 4.1813;   

%% ================================
%   WARUNKI POCZĄTKOWE
% ================================
Tb0 = 1200;
Tw0 = 25;
x0 = [Tb0; Tw0];

%% ================================
%   PARAMETRY NUMERYCZNE
% ================================
tspan = [0 5];
h = 0.0025;               % wybrany „lepszy" krok

%% ================================
%   1. SYMULACJA – EULER
% ================================
[tE, yE] = euler(@f, tspan, h, x0, params);

%% ================================
%   2. SYMULACJA – ZMODYFIKOWANY EULER
% ================================
[tM, yM] = euler_modified(@f, tspan, h, x0, params);

%% ================================
%   3. SYMULACJA – ODE45 (referencja)
% ================================
odefun = @(t,x) f(t,x,params);
[t45, y45] = ode45(odefun, tspan, x0);

%% ================================
%   4. WYKRESY: Tb(t), Tw(t)
% ================================
figure;
subplot(2,1,1); hold on;
plot(tE, yE(:,1),'b','LineWidth',1.3,'DisplayName','Euler');
plot(tM, yM(:,1),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
plot(t45, y45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Temperatura pręta T_b(t)');
xlabel('t [s]'); ylabel('T_b [C]'); legend; grid on;

subplot(2,1,2); hold on;
plot(tE, yE(:,2),'b','LineWidth',1.3,'DisplayName','Euler');
plot(tM, yM(:,2),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
plot(t45, y45(:,2),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Temperatura oleju T_w(t)');
xlabel('t [s]'); ylabel('T_w [C]'); legend; grid on;


%% =============================================================
%   5. ANALIZA WRAŻLIWOŚCI NA KROK h
%% =============================================================
H = [0.01, 0.0025];

figure; hold on;
for hk = H
    [tT, yT] = euler(@f, tspan, hk, x0, params);
    plot(tT, yT(:,1),'LineWidth',1.2,'DisplayName',['Euler h=',num2str(hk)]);
end
plot(t45, y45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Wpływ kroku h na rozwiązanie'); xlabel('t'); ylabel('T_b');
legend; grid on;


%% =============================================================
%   6. WERYFIKACJA Z DANYMI POMIAROWYMI
%% =============================================================
meas = [
    1200 25 2.5 3 107.7 105.1
    800  25 2.5 3 79.1  78.0
    1100 70 2.5 3 142.1 139.1
    1200 25 2.5 5 105.7 105.5
    800  25 2.5 5 78.2  78.1
    1100 70 2.5 2 150.1 138.2
    1100 70 5   2 116.6 105.1
    1100 70 10  2 99.1  88.1
    1100 70 2.5 4 141.2 139.8
    1100 70 2.5 5 140.9 140.1
];

fprintf("\n===== WERYFIKACJA MODEL–POMIAR =====\n");
fprintf("Nr | Tb_mod | Tb_exp | Tw_mod | Tw_exp | ErrTb | ErrTw\n");

Tb_mod = zeros(10,1);
Tw_mod = zeros(10,1);

for i = 1:10
    
    params_i = params;
    params_i.mw = meas(i,3);   % zmiana masy oleju

    [tloc, yloc] = euler_modified(@f, [0 meas(i,4)], h, meas(i,1:2)', params_i);

    Tb_mod(i) = yloc(end,1);
    Tw_mod(i) = yloc(end,2);

    fprintf("%2d | %7.2f | %7.2f | %7.2f | %7.2f | %+6.2f | %+6.2f\n",...
        i, Tb_mod(i), meas(i,5), Tw_mod(i), meas(i,6), ...
        Tb_mod(i)-meas(i,5), Tw_mod(i)-meas(i,6));
end


%% =============================================================
%   7. WRAŻLIWOŚĆ NA BŁĘDY DANYCH WEJŚCIOWYCH
%% =============================================================

Tb0_base = 1200;
Tw0_base = 25;
mw_base  = 2.5;
t_obs = 3;

Tb_ref = simulate_once(Tb0_base,Tw0_base,mw_base,t_obs,params,h);

cases = {
    'Tb0+10', Tb0_base+10, Tw0_base,   mw_base;
    'Tb0-10', Tb0_base-10, Tw0_base,   mw_base;
    'Tw0+1',  Tb0_base,    Tw0_base+1, mw_base;
    'Tw0-1',  Tb0_base,    Tw0_base-1, mw_base;
    'mw+5%',  Tb0_base,    Tw0_base,   mw_base*1.05;
    'mw-5%',  Tb0_base,    Tw0_base,   mw_base*0.95;
};

nC = size(cases,1);
dTb = zeros(nC,1);

for i = 1:nC
    dTb(i) = simulate_once(cases{i,2}, cases{i,3}, cases{i,4}, t_obs, params, h) - Tb_ref;
end

figure;
bar(dTb); xticklabels(cases(:,1));
title('Wrażliwość na błędy danych wejściowych');
ylabel('\Delta T_b [C]');
grid on;

end % ================= koniec funkcji main =====================



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                FUNKCJE – ZGODNE ZE SZKIELETEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ========================================================================
%   RÓWNANIA STANU (ODE)
%% ========================================================================
function dx = f(~, x, params)

Tb = x(1);
Tw = x(2);

h  = params.h;
A  = params.A;
mb = params.mb;
mw = params.mw;
cb = params.cb;
cw = params.cw;

dTb = -(h*A)/(mb*cb) * (Tb - Tw);
dTw =  +(h*A)/(mw*cw) * (Tb - Tw);

dx = [dTb; dTw];

end


%% ========================================================================
%   MODEL WSPÓŁCZYNNIKA h(Tb,Tw)
%   (W CZĘŚCI 1 STAŁY – docelowo w Części 2 będzie nieliniowy)
%% ========================================================================
function h = fh(Tb, Tw)
h = 160; % stały dla Części 1
end


%% ========================================================================
%   METODA EULERA
%% ========================================================================
function [t, y] = euler(func, tspan, h, x0, params)

t = tspan(1):h:tspan(2);
N = length(t);

y = zeros(N,2);
y(1,:) = x0;

for i = 1:N-1
    dy = func(t(i), y(i,:).', params);
    y(i+1,:) = y(i,:) + h * dy.';
end

end


%% ========================================================================
%   METODA ZMODYFIKOWANA EULERA (HEUN)
%% ========================================================================
function [t, y] = euler_modified(func, tspan, h, x0, params)

t = tspan(1):h:tspan(2);
N = length(t);

y = zeros(N,2);
y(1,:) = x0;

for i = 1:N-1
    
    k1 = func(t(i), y(i,:).', params);

    y_mid = y(i,:).' + (h/2)*k1;

    k2 = func(t(i)+h/2, y_mid, params);

    y(i+1,:) = y(i,:) + h * k2.';
end

end



%% ========================================================================
%   SYMULACJA DO t_obs – używana w wrażliwości
%% ========================================================================
function Tb_end = simulate_once(Tb0,Tw0,mw,t_obs,params,h)

params2 = params;
params2.mw = mw;

[t, y] = euler_modified(@f, [0 t_obs], h, [Tb0;Tw0], params2);

Tb_end = y(end,1);

end
