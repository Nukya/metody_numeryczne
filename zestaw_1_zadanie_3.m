clc
clear

syms C i R

P=(R^2 * i^2)/(C+R);

% pochodne cząstkowe
Pc=diff(P, C);
Pr=diff(P, R);
Pi=diff(P, i);

% surowe dane
rawC = 1.1;     %mF
rawR = 1;       %kOhm
rawI = 1.2;     %mA

% dane zamienione na podstawowe jednostki SI
aC = 1.1 * 10^(-3);     %F
aR = 1 *  10^3;         %Ohm
ai = 1.2 * 10^(-3);     %A

% błąd argumentu
deltaC = 0.05 * aC;
deltaR = 0.02 * aR;
deltaI = 0.01 * ai;

PcPodstawione = subs(Pc,{C, i, R},{aC, ai, aR});
PrPodstawione = subs(Pr,{C, i, R},{aC, ai, aR});
PiPodstawione = subs(Pi,{C, i, R},{aC, ai, aR});

PPodstawione = subs(P,{C, i, R},{aC, ai, aR});

% podstawienie do wzoru na błąd bezwzględny funkcji
bladBezwzglednyP = abs(PrPodstawione) * deltaR + abs(PiPodstawione) * deltaI + abs(PcPodstawione) * deltaC;

% podstawienie do wzoru na błąd względny funkcji

bladWzglednyP = (bladBezwzglednyP)/(abs(PPodstawione));


fprintf('Blad bezwzgledny funkcji wynosi: %d\n', double(bladBezwzglednyP));
fprintf('Blad wzgledny funkcji wynosi: %d\n', double(bladWzglednyP));