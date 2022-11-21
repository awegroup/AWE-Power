function [d] = createTimeseries(ws, system)
  
    d.ws     = ws;
    d.t1     = round(system.t1(d.ws),2);
    d.tROeff = round(system.tROeff(d.ws),2);
    d.t2     = round(system.t2(d.ws),2);
    d.tRIeff = round(system.tRIeff(d.ws),2);
    d.tRO    = d.t1 + d.tROeff; %[s]
    d.tRI    = d.t2 + d.tRIeff; %[s]
    d.tCycle = d.tRO+d.tRI; %[s]
    
    d.t1a_inst     = 0:0.01:d.t1;
    d.tROeff_inst  = d.t1+0.01:0.01:d.tROeff;
    d.t1b_inst     = d.tROeff+0.01:0.01:d.tROeff+d.t1;
    d.t2a_inst     = d.t1b_inst(end)+0.01:0.01:d.t1b_inst(end)+d.t2;
    d.tRIeff_inst  = d.t2a_inst(end)+0.01:0.01:d.t2a_inst(end)+d.tRIeff-d.t2;
    d.t2b_inst     = d.tRIeff_inst(end)+0.01:0.01:d.tRIeff_inst(end)+d.t2;

    d.t_inst       = [d.t1a_inst d.tROeff_inst d.t1b_inst d.t2a_inst d.tRIeff_inst d.t2b_inst];
    
    q  = fix(system.numOfPatt(ws));
    r  = mod(system.numOfPatt(ws),1);
    a1 = repmat(system.PROeff_elec_osci(ws,2:end),1,q-1);
    a2 = system.PROeff_elec_osci(ws,2:round(system.numPattParts*r));
    p  = cat(2,system.PROeff_elec_osci(ws,:),a1,a2);
    d.PRO  = interp1(1:length(p), p, linspace(1, length(p), length(d.tROeff_inst)), 'linear')./10^3;
   
    d.PRI = system.PRIeff_elec(d.ws)/10^3; %[kW]

    for i=1:length(d.t1a_inst)
        d.P1a_inst(i) = d.PRO(1)/d.t1a_inst(end)*d.t1a_inst(i);
    end
    for i=1:length(d.tROeff_inst)
        d.PROeff_inst(i) = d.PRO(i);
    end
    for i=1:length(d.t1b_inst)
        d.P1b_inst(i) = -d.PRO(end)/(d.t1b_inst(end)-d.t1b_inst(1))*(d.t1b_inst(i)-d.t1b_inst(1))+d.PRO(end);
    end
    for i=1:length(d.t2a_inst)
        d.P2a_inst(i) = (-d.PRI)/(d.t2a_inst(end)-d.t2a_inst(1))*(d.t2a_inst(i)-d.t2a_inst(1));
    end
    for i=1:length(d.tRIeff_inst)
        d.PRIeff_inst(i) = -d.PRI;
    end
    for i=1:length(d.t2b_inst)
        d.P2b_inst(i) = d.PRI/(d.t2b_inst(end)-d.t2b_inst(1))*(d.t2b_inst(i)-d.t2b_inst(1))+(-d.PRI);
    end
    d.P_inst      = [d.P1a_inst d.PROeff_inst d.P1b_inst  d.P2a_inst d.PRIeff_inst d.P2b_inst];


end