% Given parameters
Cl_maxAirfoil = 2.5; % Maximum lift coefficient
Cl0_airfoil   = 0.65; % Lift coefficient at zero angle of attack
e             = 0.6; % Oswald efficiency factor
Cd0           = 0.056; % Parasitic drag coefficient
AR            = 12; % Aspect ratio of the wing, assumed value
dCl_dalpha    = 2 * pi; % Slope of CL-alpha curve in rad^-1

% Define range of alpha in degrees and convert to radians
alpha_deg = linspace(-17, 17, 35);
alpha_rad = deg2rad(alpha_deg);

% Calculate CL for each alpha
CL = Cl0_airfoil + dCl_dalpha * alpha_rad;

% Ensure CL does not exceed Cl_maxAirfoil
CL(CL > Cl_maxAirfoil)  = Cl_maxAirfoil;
CL(CL < -Cl_maxAirfoil) = -Cl_maxAirfoil;

% Calculate CD for each CL
CD = Cd0 + ((CL - Cl0_airfoil).^2) / (pi * AR * e);

% Plot CL-alpha curve
figure('Position', [100, 100, 800, 600]);
subplot(3, 1, 1);
plot(alpha_deg, CL, 'LineWidth', 1.5);
xlabel('Angle of Attack (degrees)');
ylabel('CL');
grid on;

% Plot CD-alpha curve
subplot(3, 1, 2);
plot(alpha_deg, CD, 'LineWidth', 1.5);
xlabel('Angle of Attack (degrees)');
ylabel('CD');
grid on;

% Plot CD-CL curve
subplot(3, 1, 3);
plot(CD, CL, 'LineWidth', 1.5);
ylabel('CL');
xlabel('CD');
grid on;






