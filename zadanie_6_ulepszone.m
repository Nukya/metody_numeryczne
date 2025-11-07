close all
clear
clc

%  Dane 
dane = [
0.2 0
0.3 0
0.4 1
0.5 1
0.6 0
0.7 0
0.8 0.5
0.9 1
1.0 1.5
1.1 0.5];

x = dane(:,1);
y = dane(:,2);

n = length(x) - 1;
h = 0.1;
l = ((n+1)*h)/2;   % pół długości przedziału
c = pi / l;        % współczynnik skali dla x

%  Współczynniki 

MTY = [ sum(y); 
        sum(y .* cos(1*c*x)); 
        sum(y .* sin(1*c*x)); 
        sum(y .* cos(2*c*x)); 
        sum(y .* sin(2*c*x)); 
        sum(y .* cos(3*c*x)); 
        sum(y .* sin(3*c*x))]



MTM = diag([n+1, (n+1)/2, (n+1)/2, (n+1)/2, (n+1)/2, (n+1)/2, (n+1)/2]);

A = MTM \ MTY

a0 = A(1);
a1 = A(2);
b1 = A(3);
a2 = A(4);
b2 = A(5);
a3 = A(6);
b3 = A(7);

%  Aproksymacja 
xd = 0:0.01:max(x);
yd = a0 ...
   + a1*cos(1*c*xd) + b1*sin(1*c*xd) ...
   + a2*cos(2*c*xd) + b2*sin(2*c*xd) ...
   + a3*cos(3*c*xd) + b3*sin(3*c*xd);


%  Wykres 
plot(x, y, 'or', 'MarkerFaceColor', 'r'); hold on;
plot(xd, yd, '--b', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('Aproksymacja trygonometryczna T_m(x)');
legend('Dane pomiarowe', 'T_m(x)', 'Location', 'best');
grid on;
