clc;
clear;
close all;

% Funkcja P(C, i, R)
P = @(C, i, R) (R.^2 .* i.^2) ./ (C + R);

% Pochodne cząstkowe
Pc = @(C, i, R) -R.^2 .* i.^2 ./ (C + R).^2;
Pr = @(C, i, R) (2.*R.*i.^2.*(C + R) - R.^2.*i.^2) / (C + R).^2;
Pi = @(C, i, R) 2 .* R.^2 .* i ./ (C + R);

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

PcPodstawione = Pc(aC, ai, aR);
PrPodstawione = Pr(aC, ai, aR);
PiPodstawione = Pi(aC, ai, aR);

PPodstawione = P(aC, ai, aR);

% podstawienie do wzoru na błąd bezwzględny funkcji
bladBezwzglednyP = abs(PrPodstawione) * deltaR + abs(PiPodstawione) * deltaI + abs(PcPodstawione) * deltaC;

% podstawienie do wzoru na błąd względny funkcji

bladWzglednyP = (bladBezwzglednyP)/(abs(PPodstawione));


fprintf('Blad bezwzgledny funkcji wynosi: %d\n', double(bladBezwzglednyP));
fprintf('Blad wzgledny funkcji wynosi: %d\n', double(bladWzglednyP));
