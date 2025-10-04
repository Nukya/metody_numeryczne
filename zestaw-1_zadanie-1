A = 41/4;
t = linspace(0, 2*pi, 100);
f1 = A * sin(A * t);
f2 = A * exp(-t) .* cos(A * t);

plot(t, f1, '-o', 'DisplayName', 'f_1(t) = Asin(At)'); 
hold on; 
plot(t, f2, '-o', 'DisplayName', 'f_2(t) = Aexp(-t)cos(At)'); 
hold off; 

xlim([0 3.5])

xlabel('t [s]');
ylabel('sin(t)');
title('Kinga Kondraciuk');
legend show; 
grid on;
