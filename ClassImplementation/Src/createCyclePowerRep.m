function [cyclePowerRep] = createCyclePowerRep(ws, system, vw)
  % createCyclePowerRep: Create a representation of the power cycle for a given wind speed
  %
  % This function creates a structure `cyclePowerRep` containing relevant data
  % about the power cycle for a specific wind speed `ws` within a system `system`.
  %
  % Inputs:
  %   ws     - Wind speed of interest.
  %   system - A structure containing system parameters such as t1, t2, to_eff,
  %            ti_eff, P_e_o_eff, P_e_i_eff, P_m_o_eff, and P_m_i_eff.
  %   vw     - Vector of wind speeds corresponding to the system data.
  %
  % Outputs:
  %   cyclePowerRep - A structure containing the following fields:
  %     idx        - Index where the wind speed matches `ws`.
  %     t1         - Rounded t1 values at the specific wind speed.
  %     to_eff     - Total effective output time.
  %     t2         - Rounded t2 values at the specific wind speed.
  %     ti_eff     - Total effective input time.
  %     to         - Total output time (sum of t1 and to_eff).
  %     ti         - Total input time (sum of t2 and ti_eff).
  %     tCycle     - Total cycle time.
  %     t_inst     - Cumulative time points for each phase in the cycle.
  %     P_e_inst   - Instantaneous electrical power at key points in the cycle.
  %     P_m_inst   - Instantaneous mechanical power at key points in the cycle.
  %
  % All time and power values are rounded and provided in seconds and kilowatts, respectively.


  cyclePowerRep.ws     = ws;
  cyclePowerRep.idx    = find(vw==ws);
  cyclePowerRep.t1     = round(system.t1(cyclePowerRep.idx),2);
  cyclePowerRep.to_eff = round(sum(system.to_eff(cyclePowerRep.idx,:),2),2);
  cyclePowerRep.t2     = round(system.t2(cyclePowerRep.idx),2);
  cyclePowerRep.ti_eff = round(sum(system.ti_eff(cyclePowerRep.idx,:),2),2);
  cyclePowerRep.to     = cyclePowerRep.t1 + cyclePowerRep.to_eff; %[s]
  cyclePowerRep.ti     = cyclePowerRep.t2 + cyclePowerRep.ti_eff; %[s]
  cyclePowerRep.tCycle = cyclePowerRep.to + cyclePowerRep.to; %[s]

  cyclePowerRep.t_inst   = cumsum([0 cyclePowerRep.t1 (cyclePowerRep.to_eff-cyclePowerRep.t1) cyclePowerRep.t1 ...
    cyclePowerRep.t2 (cyclePowerRep.ti_eff-cyclePowerRep.t2) cyclePowerRep.t2]);

  cyclePowerRep.P_e_inst = [0 system.P_e_o_eff(cyclePowerRep.idx,1) system.P_e_o_eff(cyclePowerRep.idx,end) 0 ...
    -system.P_e_i_eff(cyclePowerRep.idx,end) -system.P_e_i_eff(cyclePowerRep.idx,1) 0]./10^3;

  cyclePowerRep.P_m_inst = [0 system.P_m_o_eff(cyclePowerRep.idx,1) system.P_m_o_eff(cyclePowerRep.idx,end) 0 ...
    -system.P_m_i_eff(cyclePowerRep.idx,end) -system.P_m_i_eff(cyclePowerRep.idx,1) 0]./10^3;

end