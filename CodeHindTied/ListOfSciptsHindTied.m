% 1. creates table 'TCleanSpeeds' based on talbe with all data'TWithOutliers'
% and list of outliers 'TOutListSpeeds'
CleanTable3  %v3- input from files - ended with 3, like: 'INTACT_TIED_1.0_TablesVicB3'
% output to files'INTACT_TIED_TablesVicB3' and 'SPINAL_TIED_TablesVicB3'
% output tables for each speed combination-'TClean', combined tables - 'TCleanSpeeds'

% 2. normalize EMG across speeds/conditions, saves as table variable 'EMG_Norm2' in file TIED_TablesVicB.mat
NormalizeAcrossTCond3  %v3- reads/writes from files - ended with 3
%input tables,'TCleanSpeeds', output tables-'TCleanTCond2'

% 3. synergies calculation:
GetSynCondTreadCon5 % calculate each Condition/Speeds for TIED using 'EMG_Norm2', v5- 15 muscles
GetSynAllCondSpeeds3 % synergies from combinations all Speeds/Conditions combined for TIED
%v3 - possible combinations,all speeds for one condition/ all conditions for one speed.

% 4. Plot/Sort Synergies, create tables for statistics
plotSynCondAllSpeedAll4 % v3 - add new combinations

% 5. rename 'SideSpeed' to 'IpsiSpeed'
renameTHTW_AllSpeedsAllCondr

% 6. create coeficient of determination (R2) tables for statistics
R2Tables2

% 7. Burst calc continue, change selection 'F', 'E' burst for 2-joint muscles
BusrtSelection4TwoJoint3 % remove FL, Select F/E bursts from all bursts for 2Joint muscles
%output - TCleanTCond15R with burst recalculated,%v3- reads/writes from files - ended with 3

% 8. create burst tables for statistics 
BurstTables3 %output -'TBursts2.txt', 'TCleanTCond3'; both with norm2

% 9. remove outlier bursts
CleanBurstTable %output-'TBurstClean2.txt', without ext burst of flexors

% 10. plot clusters
ClusterAnalysis

%11. plot muscles patterns
PlotMusPatterns

%12. compute and plot angles between weights (W) of different synergies 
GetAnglesBetweenW2 

%13. calculate correlations between temporal part (C or H) of different synergies 
CorrBetweenC2 % v2 -added extra R2

%14. Plot different combinations
PlotWAllCombined

%15. Create table with Kinematics for stat analysis
GetKinematicTable
