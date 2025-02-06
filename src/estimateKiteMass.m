function [m_k] = estimateKiteMass(Ft_max, S, AR)
  % estimateKiteMass: Estimates the mass of a kite for a fixed-wing airborne wind energy system
  %
  % This function estimates the mass of a kite based on a mass model developed at Ampyx Power.
  % The model uses the maximum tether force (Ft_max), wing surface area (S), and aspect ratio (AR)
  % to calculate the kite mass. The model is derived from empirical data (AP3) and projected
  % data for larger systems (AP4-AP5).
  %
  % Inputs:
  %   Ft_max - Maximum tether force [N]
  %   S      - Wing surface area [m^2]
  %   AR     - Aspect ratio [-]
  %
  % Outputs:
  %   m_k    - Estimated kite mass [kg]
  %
  % The function applies predefined coefficients and reference values based on empirical relationships
  % between the input parameters, capturing the influence of tether force, surface area, and aspect ratio.


  % Coefficients and reference values
  a1     = 0.002415;
  a2     = 0.0090239;
  b1     = 0.17025;
  b2     = 3.2493;
  k1     = 5;
  c1     = 0.46608;
  d1     = 0.65962;
  k2     = 1.1935;
  AR_ref = 12;

  % Calculate intermediate values
  a      = a1 * (Ft_max/1000 / S) + a2;
  b      = b1 * (Ft_max/1000 / S) + b2;

  % Estimate kite mass
  m_k = 10 * (a * S^2 + b * S - k1) * (c1 * (AR / AR_ref)^2 - d1 * (AR / AR_ref) + k2);
end
