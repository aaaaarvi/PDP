clear;

% Data vectors
N_vect = [50,      100,     200,     400,     800,     1600,    3200,    6400];
e_vect = [1.74e-2, 2.36e-3, 2.01e-4, 3.35e-5, 3.71e-6, 2.62e-7, 1.68e-8, 1.06e-9];

% Plot the error
figure;
loglog(1./N_vect, e_vect, '*b');
xlim([1e-4, 4e-2]);
ylim([1e-10, 1e-1]);
grid on;
xlabel('step size dx');
ylabel('error');

figure;
plot(log10(1./N_vect), log10(e_vect));
grid on;