% WESC2023 abstract plots

%% Wind distribution at average pattern height
ws        = 0:1:30;
windProb  = wblpdf(ws,11,2.2);
windOccur = 8760*windProb;
meanWS    = mean(windOccur)/31;

%% AEP, CF calculation
P_Rated = system.ratedPower
AEP     = 8760*sum(windProb(1:25).*system.Pcycle_elec)/10^6; %[MWh]
CF      = AEP/P_Rated/8760*10^6

