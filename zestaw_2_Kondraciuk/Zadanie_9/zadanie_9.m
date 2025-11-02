clear;
close all;
clc;

%  Układ nieliniowy 
F = @(x) [x(1)*x(2) - x(1) - 2;
          x(1)^2 - x(2)^2 - 1];

%  Jakobian 
J = @(x) [x(2) - 1,     x(1);
          2*x(1),       -2*x(2)];

%  Punkt startowy 
%    x1    x2
x = [-0.6; -0.65];

%  Iteracyjna metoda Newtona 
i = 1;
ff = 30;

for i = 1:1000
    ff = F(x);
    if norm(ff) < 1e-8
        break
    end
    x = x - J(x)\ff;
end


%  Siatka dla powierzchni 
x1 = -2:0.1:2;
x2 = -2:0.1:2;

ff1 = zeros(length(x1), length(x2));
ff2 = zeros(length(x1), length(x2));

for ii = 1:length(x1)
    for jj = 1:length(x2)
        ff = F([x1(ii); x2(jj)]);
        ff1(ii,jj) = ff(1);
        ff2(ii,jj) = ff(2);
    end
end

[X1, X2] = meshgrid(x1, x2);

%  Rysowanie powierzchni 
attr = {'Interpreter','latex'};

figure;
surf(X1, X2, ff1');
hold on;
surf(X1, X2, ff2');

scatter3(x(1), x(2), 0, 60, 'r', 'filled');

annotation('textarrow', [0.52 0.36], [0.82 0.65], ...
    'String', sprintf('(%.3f, %.3f)', x(1), x(2)), ...
    'FontSize', 10, 'FontWeight', 'bold', attr{:});

xlabel('$x_1$', attr{:});
ylabel('$x_2$', attr{:});
zlabel('$f(x_1,x_2)$', attr{:});
view([-312.5 21.3]);
grid on;
hold off;

% saveas(gcf,'zadanie_9.png')