function [m_k] = estimateKiteMass_awePower(Ft_max, S, AR)

  % Vincent Bonnin's simple mass model developed at Ampyx Power. Based on AP3 data and projected data for larger systems (AP4-AP5)      
  a1     = 0.002415;       a2     = 0.0090239;       b1     = 0.17025;       b2     = 3.2493;
  k1     = 5;              c1     = 0.46608;         d1     = 0.65962;       k2     = 1.1935;
  AR_ref = 12;
  a = a1*(Ft_max/S) + a2;
  b = b1*(Ft_max/S) + b2;
  m_k = 10*(a*S^2 +b*S-k1)*(c1*(AR/AR_ref)^2-d1*(AR/AR_ref)+k2); 

end