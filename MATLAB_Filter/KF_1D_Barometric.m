%% Kalman Filter Design 1 - USST
%{
Description: 

This Kalman Filter is an attempt to estimate and track the
altitude of a model rocket over time using simulated barometric pressure
measurements.

The filter uses a 1D constant acceleration motion model to describe the
vertical dynamics. Its state vector consists of three variables: altitude,
vertical velocity, and vertical acceleration.

The prediction step propagates the estimated state forward in time using
the equations of motion. The correction step incorporates altitude
measurements derived from barometric pressure using the International
Standard Atmosphere (ISA) barometric equation.

Measurement noise is introduced to simulate imperfections in the pressure
sensor. The KF then combines the predicted state with the measured altitude
to produce a smoothed altitude estimate.

The KF currently has constant quantization noise. Future addition is to
randomize this noise and observe the behavior of the filter.

For reference equations and theory, refer to: 
https://kalmanfilter.net/multiSummary.html
%}
close all;
clear;
clc;

%% 1 Dimensional KF

%{
Motion Equations:
    a_{t + ∆t} = a_t 
    v_{t + ∆t} = v_t + a_t * ∆t    
    h_{t + ∆t}  = h_t + v_t * ∆t + a_t * (∆t^2)/2
Pressure Equations:
    P = P_0 * (1 - (h * L / T_0)) ^ (g * M / (R * L))
    h = (T_0 / L) * ((P / P_0) ^ -(R * L / (g * M)) - 1)
%}


%%%%%%%%%%%%%
% Variables %
%%%%%%%%%%%%%

% Time Vector
%------------
Tf = 15;        % Stop Time
Ts = 0;         % Start Time
dt = 1e-2;      % Time Step
t = Ts:dt:Tf;
N = length(t);

% Natural
%--------
g = 9.80665;    % Gravity - m/s^2
T0 = 288.15;    % Sea Level Temp (Ideal 288.15 K) - K
P0 = 101325;    % Sea Level Pressure (Ideal 101.325 kPa) - Pa
L = 6.5e-3;     % Lapse Rate (Ideal 6.5 C/km) - K/m
M = 0.0289644;  % Molar Mass of Air (Ideal) - kg/mol
R = 8.31446;    % Gas Constant (Ideal) - J/(mol*K)

% Generation Rocket SIM
%----------------------
h_true = zeros(size(t));    % Height
v_true = zeros(size(t));    % Velocity
a_true = zeros(size(t));    % Acceleration

for k = 2:length(t)
    if t(k) < 3
        a_true(k) = 30; % thrust phase --> ~3 seconds of acceleration
    else
        a_true(k) = -g; % gravity only
    end

    v_true(k) = v_true(k-1) + a_true(k-1)*dt;
    h_true(k) = h_true(k-1) + v_true(k-1)*dt + a_true(k-1)*dt;
end

%{

% Plot Test Data
%---------------
[y, x] = max(h_true);   % Max Height

subplot(3, 1, 1);
title("Height vs Time");
plot(t, h_true);
xlabel("Time (s)");
ylabel("Height (m)");
grid on;

subplot(3, 1, 2);
title("Velocity vs Time");
plot(t, v_true);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
grid on;

subplot(3, 1, 3);
title("Acceleration vs Time");
plot(t, a_true);
xlabel("Time (s)");
ylabel("Acceleration (m/s^2)");
grid on;

figure();
title("Test Data");
plot(t, h_true, "-r");
hold on;
plot(t, v_true, "--b");
plot(t, a_true, "-.m");
plot(t(x), y, "ok");
xlabel("Time (s)");
ylabel("Amplitude");
legend("Height", "Velocity", "Acceleration", ...
    sprintf("Max Height = %.2f", y), 'Location', 'northwest');
grid on;
%}


%%%%%%%%%%%%%%%%%
% Kalman Filter %
%%%%%%%%%%%%%%%%%

% Functions
%----------
Motion = @(x)[
    x(1) + x(2)*dt + 0.5*x(3)*dt^2;
    x(2) + x(3)*dt;
    x(3)];
Ptoh = @(P) (T0 / L) * ((P ./ P0).^-(R * L / (g * M)) - 1);

% Filter Variables
%-----------------
x0 = [0; 0; 0];  % x(1) - Height, x(2) - Velocity, x(3) - Acceleration
p0 = 8;     % Covariance - Pressure to Amplitude error
q = 0.3;    % Quntization Noise
r0 = 0.09;  % Measurment Variance - Quantization Noise Squared

x_pred = zeros(3, N);
x_update = zeros(3, N);
p_pred = zeros(size(t));
p_update = zeros(size(t));
K = 0;

x_pred(:, 1) = x0;
p_pred(1) = p0;
P_meas = P0 * (1 - (h_true * L / T0)).^(g * M / (R * L));
h_meas = Ptoh(P_meas) + q;

% Filter Loop
%------------
for k = 2:length(t)
    % Measure
    z = h_meas(k);

    % Update
    K = p_pred(k-1) / (p_pred(k-1) + r0);
    x_update(:, k-1) = x_pred(:, k-1) + K*(z - x_pred(:, k-1));
    p_update(k-1) = (1 - K)*p_pred(k-1);
    
    % Predict
    x_pred(:, k) = Motion(x_update(:, k-1));
    p_pred(k) = p_pred(k-1);
end

plot(t, h_true, "-r");
hold on;
plot(t, x_update(1, :), "--b")
legend("True height", "Predicted height");