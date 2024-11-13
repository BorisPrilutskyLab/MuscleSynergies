TBurstSpeed = TBurst(TBurst.IpsiSpeed==0.4,:);
Muscle = 'GLU';
Function = 'Ext';
TBurstSpeedMus = TBurstSpeed(ismember(TBurstSpeed.Muscle,Muscle)&ismember(TBurstSpeed.Function,Function),:);

TBurstSpeedMusInt = TBurstSpeedMus(ismember(TBurstSpeedMus.Condition,'INTACT'),:);
TBurstSpeedMusSP = TBurstSpeedMus(ismember(TBurstSpeedMus.Condition,'SPINAL'),:);

mean(TBurstSpeedMusInt.BurstOnNorm)
mean(TBurstSpeedMusSP.BurstOnNorm)
