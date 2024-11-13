clear variables
DataFolder = '..\DataHindTied\';

doSaveFigures = true;%false;
TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';

TBurstsFilename = strcat(DataFolder,TreadmillCondition,'_','TBurstClean.txt');
TBurst = readtable(TBurstsFilename);

% TBurstVars = {'Muscle','Condition','IpsiSpeed','ContraSpeed','Function',...
%     'CycleTime','StanceTime','DF','BurstOnNorm','BurstOffNorm',...
%     'BurstAveNormMax2','BurstAveNormPeakEnv2'};

Conditions = unique(TBurst.Condition);
IpsiSpeeds = unique(TBurst.IpsiSpeed);

iCond = 1;
Condition = Conditions{iCond};
iSpeed = 1;
IpsiSpeed = IpsiSpeeds(iSpeed);

TCondSpeed = TBurst(ismember(TBurst.Condition,Condition)&ismember(TBurst.IpsiSpeed,IpsiSpeed),:);
Muscles = unique(TCondSpeed.Muscle);
Functions = unique(TBurst.Function);

Bursts = cell(1,length(Muscles)*length(Functions));
COnOffs = cell(1,length(Muscles)*length(Functions));%tables 'on/off' for each burst
for iMus = 1:length(Muscles)
    Muscle = Muscles{iMus};
    for iFunc = 1:length(Functions)
        Function = Functions{iFunc};
        BurstName = strcat(Muscle,{'.'},Function);
        iBurst = (iFunc-1)*length(Muscles) + iMus;
        Bursts(iBurst) = BurstName;
        COnOffs{iBurst} = TCondSpeed(ismember(TCondSpeed.Muscle,Muscle)&...
            ismember(TCondSpeed.Function,Function),{'BurstOnNorm','BurstOffNorm'});
    end
end

EmptyBurstsID = cellfun(@isempty,COnOffs);
Bursts(EmptyBurstsID) = [];
COnOffs(EmptyBurstsID) = [];

% Outliers
OutBurst = {'IP.Ext','TA.Ext','BFA.Flex','MG.Flex','VL.Flex','SOL.Flex'};

COnOffs(ismember(Bursts,OutBurst)) = [];
Bursts(ismember(Bursts,OutBurst)) = [];



%%

NSamplesMax = max(cellfun(@height,COnOffs));
cut_data = zeros(NSamplesMax,2*length(Bursts));
M = zeros(length(Bursts),2);
R = zeros(length(Bursts),2);
for iBurst =1:length(Bursts)
    ind = 2*iBurst-1:2*iBurst;
    OnOffs = table2array(COnOffs{iBurst}); % 2D clouds for each burst
    cut_data(1:height(COnOffs{iBurst}),ind) = OnOffs;% matrix raws - samples, columns- on/off for all burst = 2*Nbursts
    M(iBurst,:) = mean(OnOffs); % centers of clouds for each burst
    R(iBurst,:) = std(OnOffs); % axis of ellipse of clouds for each burst
end
[mus_table,v_len] = remove_zerosAK( cut_data );
% mus_table - all clouds of on/off reshaped (nAllSamples, 2); v_len -
% Nsamples for each burst

%plot parameters
Swing_b =  mean((TCondSpeed.CycleTime - TCondSpeed.StanceTime)./TCondSpeed.CycleTime);
Ext_b_ind = find(ismember(Bursts,'SOL.Ext'));
Flex_b_ind = find(ismember(Bursts,'IP.Flex'));
Flex_b = M(Flex_b_ind)-R(Flex_b_ind);
Ext_b = M(Ext_b_ind)-R(Ext_b_ind);

% connections parameters
w_thresh = 1.0;
c_frac = 1;                                          % cFrac parameter (see Krouchev et al, J Neurophys. 2006)
wc_thresh = .00;                                     % weak connection threshold
sc_thresh = .5;                                      % strong connection threshold

% overlapping clouds with centers in 'M'; 'R'- is 2D SD
% area of overlapping parts of cloudes normalized to mean radius of overlaping clouds
ovr = meancalc2( M, R, c_frac, wc_thresh );

n_eng = length( Bursts );

% define significant distance between vertices
graphB = 1./ovr;
graphB(isinf(graphB)|isnan(graphB)) = 0;
SG = sparse( graphB );
G = graph(SG,Bursts);
%
%S0 - number of strongly connected nodes,C0 - which claster node belongs to
[S0,C0] =  conncomp( G );%, 'Directed', false
%[S0,C0] =  graphconncomp( SG, 'Directed', false );
%h = view( graph( G, Bursts));%, 'ShowArrows','off','ShowWeights','on' 
figure('Name',strcat(Condition,'.',num2str(IpsiSpeed),'.','_GraphWithoutTreshold'))
plot(G)
% find the shortest path between vertices and calculate average distance between vertices
thresh = inf( n_eng,1 );
for eng = 1:n_eng
    %dist = graphshortestpath( SG, eng, 'Directed', true, 'Method', 'Dijkstra' );
    %    TR = shortestpathtree(G,s)
    %    P = shortestpath(G,s,t)
    TR = shortestpathtree(G,Bursts(eng) );
    dist = TR.Edges.Weight;
    thresh( eng ) = w_thresh* mean(dist( ~(isinf(dist)|dist ==0)));
end

ovr0 =ovr; %test
ovr(graphB > thresh)= 0;%set overlap to zero for distances > mean distance
ovr((graphB > thresh)') = 0;%for symmetry
ovr(ovr>0) = 1;%set overlap to '1' for all others

G1 = graph( ovr ,Bursts);
[S,C] =  conncomp(G1);
%
figure('Name',strcat(Condition,'.',num2str(IpsiSpeed),'.','_GraphUsingThreshold'))
h = plot(G1,'Layout','layered');
if doSaveFigures
    SpeedStr = num2str(IpsiSpeed);
    ax = gca;
    FigName = strcat('../Figures/Fig_Graph_',Condition,'_',SpeedStr,'.eps');
    exportgraphics(ax,FigName,'BackgroundColor','none','ContentType','vector')
    FigName = strcat('../Figures/Fig_Graph_',Condition,'_',SpeedStr,'.pdf');
    exportgraphics(ax,FigName,'BackgroundColor','none','ContentType','vector')
end

n_clust = length(C); 
Clustcolors = lines(n_clust );
for iClust = 1:n_clust
    highlight(h,iClust==S,'NodeColor',Clustcolors(iClust,:))
end
%n_clust = S;
[data,vlength] = remove_zerosAK( cut_data );

r_ellipse = zeros(100,2,n_clust);
for clust = 1:n_clust 
    iCmus = find(S==clust);
    %iCmus = find(C==clust); % index of muscles for this cluster
    idata = [];% index of muscles for this cluster in 'cut_data'
    for im = 1:length(iCmus)
        idata=[idata,iCmus(im)*2-1:2*iCmus(im)];
    end
    [merge_data,vlengthC] = remove_zerosAK( cut_data(:,idata) );
        
    %        r_ellipse
    covariance = cov(merge_data);
    [eigenvec, eigenval ] = eig(covariance);

    % Get the index of the largest eigenvector
    [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
    largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
    
    % Get the largest eigenvalue
    largest_eigenval = max(max(eigenval));
    
    % Get the smallest eigenvector and eigenvalue
    if(largest_eigenvec_ind_c == 1)
        smallest_eigenval = max(eigenval(:,2));
        smallest_eigenvec = eigenvec(:,2);
    else
        smallest_eigenval = max(eigenval(:,1));
        smallest_eigenvec = eigenvec(1,:);
    end
    
    % Calculate the angle between the x-axis and the largest eigenvector
    angle1 = atan2(largest_eigenvec(2), largest_eigenvec(1));
    
    % This angle is between -pi and pi.
    % Let's shift it such that the angle is between 0 and 2pi
    if(angle1 < 0)
        angle1 = angle1 + 2*pi;
    end
    
    % Get the coordinates of the data mean
    avg = mean(merge_data);
    
    % Get the 95% confidence interval error ellipse
    %Using link
    %https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
    % Set Given probability Q= 1-confidence interval;Degree of freedom d=2;
    % Will have "chisquare" value then take sqrt from it that will be chisquare_val
    %chisquare_val = 2.4477; %95% corresponds to 2 SD
    %chisquare_val = 3.0348; %99%
    %chisquare_val = 2.1459; %90%
    chisquare_val = 1.5096; %68% corresponds to 1 SD
    theta_grid = linspace(0,2*pi);
    phi = angle1;
    X0(clust)=avg(1);
    Y0(clust)=avg(2);
    a=chisquare_val*sqrt(largest_eigenval);
    b=chisquare_val*sqrt(smallest_eigenval);
    % the ellipse in x and y coordinates
    ellipse_x_r  = a*cos( theta_grid );
    ellipse_y_r  = b*sin( theta_grid );
    
    %Define a rotation matrix
    Rot = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    
    %let's rotate the ellipse to some angle phi
    r_ellipse(:,:,clust) = [ellipse_x_r;ellipse_y_r]' * Rot;
    
end
%%
assign_colors16
n_mus = length(Bursts);
figure('Name',strcat(Condition,'.',num2str(IpsiSpeed),'.','_Clusters'))
hold on
count = 1;
for( i = 1:n_mus ) 
    for( j = 1:v_len( i )) 
        if sh( i*2-1)=='c'
            hh(i)=plot( mus_table( count,1 ), mus_table( count,2 ),'Color',[0,0.7,0.9],'Marker', sh( 2*i ));
        else
            hh(i)=plot( mus_table( count,1 ), mus_table( count,2 ), sh( i*2-1:2*i ));
        end
        count=count+1;
    end
end
hleg1 =legend (hh(1:n_mus),Bursts(1:n_mus));
%set(hleg1,'Location','SouthEast','FontSize',12);
drawnow;
axis([-0.3 1.2 -0.1 1.2]);% equal;

hold on
for clust = 1: n_clust
    he(clust) = plot(r_ellipse(:,1,clust) + X0(clust),r_ellipse(:,2,clust) + Y0(clust),'k-','linewidth',2);
    he(clust).Annotation.LegendInformation.IconDisplayStyle = 'off';
end
hl(1) = line([0,0],[-0.2,1.2], 'LineStyle','-.', 'Color','b');
hl(2) = line([Swing_b,Swing_b],[-0.2,1.2], 'LineStyle','-.', 'Color','b');
hl(3) = line([1.0,1.0],[-0.2,1.2], 'LineStyle','-.', 'Color','b');
hl(4) = line([Flex_b,Flex_b],[-0.2,1.2], 'LineStyle','--', 'Color',[.7 .7 .7]);
hl(5) = line([Ext_b,Ext_b],[-0.2,1.2], 'LineStyle','--', 'Color',[.7 .7 .7]);

hl(6) = line([-0.3,1.3],[0,0], 'LineStyle','-.', 'Color','b');
hl(7) = line([-0.3,1.3],[Swing_b,Swing_b], 'LineStyle','-.', 'Color','b');
hl(8) = line([-0.3,1.3],[1.0,1.0], 'LineStyle','-.', 'Color','b');
hl(9) = line([-0.3,1.3],[Ext_b,Ext_b], 'LineStyle','--', 'Color',[.7 .7 .7]);
hl(10) = line([-0.3,1.3],[Flex_b+1,Flex_b+1], 'LineStyle','--', 'Color',[.7 .7 .7]);
for iLine =1:length(hl)
    hl(iLine).Annotation.LegendInformation.IconDisplayStyle = 'off';
end
hold off;

if doSaveFigures
    SpeedStr = num2str(IpsiSpeed);
    ax = gca;
    FigName = strcat('../Figures/Fig_Clusters_',Condition,'_',SpeedStr,'.eps');
    exportgraphics(ax,FigName,'BackgroundColor','none','ContentType','vector')
    FigName = strcat('../Figures/Fig_Clusters_',Condition,'_',SpeedStr,'.pdf');
    exportgraphics(ax,FigName,'BackgroundColor','none','ContentType','vector')
end
