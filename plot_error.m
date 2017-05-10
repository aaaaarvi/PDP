clear;

% Data vectors
N_vect = [50,      100,     200,     400,     800,     1600];
e_vect = [9.54e-3, 2.36e-3, 5.87e-4, 1.46e-4, 3.65e-5, 9.13e-6];

% Fit a linear function to the log values
p = polyfit(log(1./N_vect), log(e_vect), 1);

% Plot the error
figure;
loglog(1./N_vect, e_vect, '*b');
xlim([3e-4, 4e-2]);
ylim([3e-6, 2e-2]);
grid on;
xlabel('step size dx');
ylabel('error');
legend(['convergence rate = ' num2str(p(1))], 'Location', 'northwest');

