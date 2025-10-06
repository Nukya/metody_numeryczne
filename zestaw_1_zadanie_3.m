
clc

syms C i R

P=(R^2 * i^2)/(C+R);

% pochodne cząstkowe
Pc=diff(P, C);
Pr=diff(P, R);
Pi=diff(P, i);

disp(Pr)
% surowe dane
rawC = 1.1; %mF
rawR = 1; %kOhm
rawI = 1.2; %mA

% dane zamienione na podstawowe jednostki SI
aC = 1.1 * 10^(-3);
aR = 1 *  10^3;
ai = 1.2 * 10^(-3);

% błąd argumentu
deltaC = 0.05 * rawC;
deltaR = 0.02 * rawR;
deltaI = 0.01 * rawI;

PcPodstawione = subs(Pc,{C, i, R},{aC, ai, aR});
PrPodstawione = subs(Pr,{C, i, R},{aC, ai, aR});
PiPodstawione = subs(Pi,{C, i, R},{aC, ai, aR});


% podstawienie do wzoru na błąd bezwzględny funkcji
deltaP = abs(PrPodstawione) * deltaR + abs(PiPodstawione) * deltaI + abs(PcPodstawione) * deltaC

% podstawienie do wzoru na błąd względny funkcji
% wzglednaDeltaP