clearvars

for i_n = [1,2]

fname = sprintf('Data_%d.mat', i_n);
load(fname)

% Transparency
T = 0.02; 
% TimeStamp number
TSn = 20; 

for i = 1:3
    
    X_TimeStamp = zeros(length(simData.DataTable_Averages(i,:)), TSn);
    Y_TimeStamp = zeros(length(simData.DataTable_Averages(i,:)), TSn);
    TimeRecord = zeros(length(simData.DataTable_Averages(i,:)), TSn);
    
    for j = 1:length(simData.DataTable_Averages(i,:)) 
    
        Traj = simData.DataTable_Averages{i,j};    
        TimeStamp = floor( linspace(1, length( Traj(1,:) ), TSn) );
        TimeRecord(j,:) = TimeStamp;
        X_TimeStamp(j,:) = Traj(1,TimeStamp);
        Y_TimeStamp(j,:) = Traj(2,TimeStamp); 
    
        if i_n == 1
            if j < 101      
                figure(1)      
                hold on       
                p1 = plot( Traj(1,:), Traj(2,:),'color','k','LineWidth',1 );       
                p1.Color(4) = T;
            end
            
            Record = Theory.Record{i,1}; 
            plot(Record(1,:),Record(2,:),'-','color','r', 'LineWidth',0.75)  
        end

        if i_n == 2
            if j < 101      
                figure(2)      
                hold on       
                p2 = plot( Traj(1,:), Traj(2,:),'color','k','LineWidth',1 );     
                p2.Color(4) = T;
            end                        
            Record = Theory.Record{i,1}; 
            plot(Record(1,:),Record(2,:),'-','color','r', 'LineWidth',0.75)
        end 

    end
    
    Mean_TimeStamp = zeros(2,TSn);
    
    for i_mean = 1:TSn
        Mean_TimeStamp(1,i_mean) = mean(X_TimeStamp(:,i_mean));
        Mean_TimeStamp(2,i_mean) = mean(Y_TimeStamp(:,i_mean));
    end
    
    if i_n == 1
        figure(1)
        hold on
        scatter(Mean_TimeStamp(1,:),Mean_TimeStamp(2,:),15,'MarkerEdgeColor','k',...
              'MarkerFaceColor','k')
        %line(Mean_TimeStamp(1,:),Mean_TimeStamp(2,:),'color','k')
    end

    if i_n == 2
        figure(2)
        hold on
        scatter(Mean_TimeStamp(1,:),Mean_TimeStamp(2,:),15,'MarkerEdgeColor','k',...
              'MarkerFaceColor','k')
          
    end
    
    MeanTrajectory = [];
    save(fname, 'MeanTrajectory', '-append')

    MeanTrajectory.X_TimeStamp = X_TimeStamp;
    MeanTrajectory.Y_TimeStamp = Y_TimeStamp;
    MeanTrajectory.Mean_TimeStamp = Mean_TimeStamp;
    MeanTrajectory.TimeRecord = TimeRecord;

    save(fname, 'MeanTrajectory', '-append')
end
end

%% Aesthetics
% Run this part for each figure obtained by running the loop above.

set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')

ax = gca;
ax.Box = 'off'; % switches off the surrounding box

H=gca;
H.LineWidth=1.2; %change to the desired value 
set(gca,'fontsize',14);
xticks(-3.2: 0.1: 0.25)
yticks(-3.2: 0.1: 0.25)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

xlim([-3.2,0.1])
ylim([-3.2,0.1])

outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

f = @(x,y) abs(x) + abs(y);
grayColor = [.7 .7 .7];
fc = fcontour( f,'LineColor', grayColor );
fc.LevelList = [4 2 1];
fc.LineStyle = '--';
fc.LineWidth = 1;

hold on
