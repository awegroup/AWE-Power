function [ts] = createTimeseries(ws, system)
  
    ts.ws     = ws;
    ts.t1     = round(system.t1(ts.ws),2);
    ts.tROeff = round(sum(system.tROeff(ts.ws,:),2),2);
    ts.t2     = round(system.t2(ts.ws),2);
    ts.tRIeff = round(sum(system.tRIeff(ts.ws,:),2),2);
    ts.tRO    = ts.t1 + ts.tROeff; %[s]
    ts.tRI    = ts.t2 + ts.tRIeff; %[s]
    ts.tCycle = ts.tRO+ts.tRI; %[s]
    
    ts.t_inst = cumsum([0 ts.t1 system.tROeff(ts.ws,:) ts.t1 ts.t2 system.tRIeff(ts.ws,:) ts.t2 0]);
    
    ts.P_e_inst = [0 system.PROeff_elec(ts.ws,1) system.PROeff_elec(ts.ws,:) 0 ...
                -flip(system.PRIeff_elec(ts.ws,1)) -flip(system.PRIeff_elec(ts.ws,:)) 0 0]./10^3;
              
    ts.P_m_inst = [0 system.PROeff_mech(ts.ws,1) system.PROeff_mech(ts.ws,:) 0 ...
                -flip(system.PRIeff_mech(ts.ws,1)) -flip(system.PRIeff_mech(ts.ws,:)) 0 0]./10^3;
    
end