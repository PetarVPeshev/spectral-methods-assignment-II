close all;
clear;
clc;

addpath('../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
dipole.L  = 0.25e-3;
dipole.W  = 0.25e-3;
dielectric.h  = 0.8e-3;
dielectric.er = 12;
medium.er = 1;
wave.f  = 28e9;
% FF parameters
R_FF = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% SPHERICAL COORDINATE SYSTEM
phi  = (eps:2:360) * pi / 180;
theta = linspace(eps, 90, 90) * pi / 180;
sph_grid = meshgrid_comb(theta, phi);
Z = R_FF * cos(sph_grid(:, :, 1));

%% WAVE VECTOR
kx = wave.k0 * sin(sph_grid(:, :, 1)) .* cos(sph_grid(:, :, 2));
ky = wave.k0 * sin(sph_grid(:, :, 1)) .* sin(sph_grid(:, :, 2));
krho = sqrt(kx .^ 2 + ky .^ 2);

%% STRATIFIED MEDIA VOLTAGE AND CURRENTS
[v_te, i_te, v_tm, i_tm] = stratified_media(wave.k0, krho, Z, ...
    'GroundSlab', dielectric.h, dielectric.er);

%% PLOT TM VOLTAGE
figure('Position', [250 250 800 400]);
plot(sin(theta), 20 * log10( abs(v_tm(1, :)) ), ...
    'LineWidth', 2.0, 'DisplayName', 'v_{tm}');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('x / m');
ylabel('v_{tm} / dB');
