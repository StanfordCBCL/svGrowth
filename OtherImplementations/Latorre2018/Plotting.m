clc; clear all; clf;

data = load('GnR_out.txt');
time = -140:1:560;
a = data(1:end - 1,1);
h = data(1:end - 1,2);
rho_m = data(1:end - 1,3);
rho_c = data(1:end - 1,4);
ups_c = data(1:end - 1,5);
perturb = data(1:end - 1,6);

% rho_m = data(1:end - 1,3);
% rho_c = data(1:end - 1,4);
a_e = data(end,1);
h_e = data(end,2);
rho_m_e = data(end,3);
rho_c_e = data(end,4);
ups_c_e = data(end,5);
perturb_e = data(end,6);

figure(1)
subplot(2,3,1)
hold on
plot(time,a/a(1))
plot(time(end), a_e/a(1), 's')
xlabel('Time (days)'); ylabel('radius (-)')
axis([ -50 560 0.95 1.2 ])

subplot(2,3,2)
hold on
plot(time,h/h(1))
plot(time(end), h_e/h(1), 's')
xlabel('Time (days)'); ylabel('thickness (-)')
axis([ -50 560 0.8 2.0])

subplot(2,3,3)
hold on
plot(time,ups_c/ups_c(1))
plot(time(end), ups_c_e/ups_c(1), 's')
xlabel('Time (days)'); ylabel('Ups_c (-)')
axis([ -50 560 0.95 1.3])

subplot(2,3,4)
hold on
plot(time,rho_m/rho_m(1))
plot(time(end), rho_m_e/rho_m(1), 's')
xlabel('Time (days)'); ylabel('rhoR_m (-)')
axis([-50 560 0.75 3.5])

subplot(2,3,5)
hold on
plot(time,rho_c/rho_c(1))
plot(time(end), rho_c_e/rho_c(1), 's')
xlabel('Time (days)'); ylabel('rhoR_c (-)')
axis([ -50 560 0.75 3.5])

subplot(2,3,6)
hold on
plot(time,perturb/perturb(1))
plot(time(end), perturb_e/perturb(1), 's')
xlabel('Time (days)'); ylabel('Perturb (-)')
axis([ -50 560 0.95 1.6])
