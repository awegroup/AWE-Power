% Given parameters
Cl_maxAirfoil = 2.5; % Maximum lift coefficient
Cl0_airfoil   = 0.65; % Lift coefficient at zero angle of attack
Cd0           = 0.056; % Parasitic drag coefficient
AR            = 6; % Aspect ratio of the wing, assumed value
% e             = 0.6; % Oswald efficiency factor
e = 1.78*(1-0.045*AR^0.68)-0.64;
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
% CD = Cd0 + ((CL - Cl0_airfoil).^2) / (pi * AR * e);

% Considering Tether drag for a 500m long tether
CD = Cd0 + ((CL - Cl0_airfoil).^2) / (pi * AR * e) + (1/4)*1.2*0.0087*600/12;


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

% Plot CL/CD vs CL
figure()
hold on
grid on
plot(CL,CL./CD, 'LineWidth', 1.5)
xlabel('CL');
ylabel('CL/CD');
hold off




