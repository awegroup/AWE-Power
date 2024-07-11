function [cyclePowerRep] = createCyclePowerRep(ws, system, vw)

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