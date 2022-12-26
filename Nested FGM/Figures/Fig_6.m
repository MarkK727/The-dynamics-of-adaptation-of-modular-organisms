%% Figure 6

greycolor = [.7, .7, .7];

for i_n = 1:3
fname = sprintf('Data_%d.mat', i_n);
load(fname)

figure(1)

Dep_Variable = simData.Summary(1,:);
std = simData.Summary(2,:);

xlim([0,pi/2])
ylim([0,pi/2])

if i_n == 1
    Indep_Variable = simData.Initial_Angle + 0.02;
    plt = errorbar(Indep_Variable,Dep_Variable,std,'o-','color','#56B4E9','MarkerFaceColor','#56B4E9','LineWidth', 3/4,'CapSize',0);
    xlim([0,pi/2])
    ylim([0,pi/2])
    
    ESL = [FixedParam.SelectionBias(2), FixedParam.SelectionBias(1)]';
    hold on
    yline(atan(abs(ESL(2)/ESL(1))), '--','LineWidth',1.5, 'color',greycolor);
    
    xticks(0*pi:(1/12)*pi:(1/2)*pi)
    xticklabels({'0\pi','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});
    yticks(0*pi:(1/12)*pi:(1/2)*pi)
    yticklabels({'','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});
end

if i_n == 2
    Indep_Variable = simData.Initial_Angle + 0.01;
    plt = errorbar(Indep_Variable,Dep_Variable,std,'o-','color','#009E73','MarkerFaceColor','#009E73','LineWidth', 3/4,'CapSize',0);
    xlim([0,pi/2])
    ylim([0,pi/2])
    
    ESL = [FixedParam.SelectionBias(2), FixedParam.SelectionBias(1)]';
    hold on
    yline(atan(abs(ESL(2)/ESL(1))), '--','LineWidth',1.5, 'color',greycolor);
    
    xticks(0*pi:(1/12)*pi:(1/2)*pi)
    xticklabels({'0\pi','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});
    yticks(0*pi:(1/12)*pi:(1/2)*pi)
    yticklabels({'','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});    
end

if i_n == 3
    Indep_Variable = simData.Initial_Angle - 0.0;
    plt = errorbar(Indep_Variable,Dep_Variable,std,'o-','color','#D55E00','MarkerFaceColor','#D55E00','LineWidth', 3/4,'CapSize',0);
    xlim([0,pi/2])
    ylim([0,pi/2])
    
    ESL = [FixedParam.SelectionBias(2), FixedParam.SelectionBias(1)]';
    hold on
    yline(atan(abs(ESL(2)/ESL(1))), '--','LineWidth',1.5, 'color',greycolor);
    
    xticks(0*pi:(1/12)*pi:(1/2)*pi)
    xticklabels({'0\pi','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});
    yticks(0*pi:(1/12)*pi:(1/2)*pi)
    yticklabels({'','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});
end

figure(2) 
hold on

Dep_Variable = simData.Summary(3,:);
std = simData.Summary(4,:);

if i_n == 1
    Indep_Variable = simData.Initial_Angle + 0.02;
    plt2 = errorbar(Indep_Variable,Dep_Variable,std,'o-','color','#56B4E9','MarkerFaceColor','#56B4E9','LineWidth', 3/4,'CapSize',0);
end

if i_n == 2
    Indep_Variable = simData.Initial_Angle + 0.01;
    plt2 = errorbar(Indep_Variable,Dep_Variable,std,'o-','color','#009E73','MarkerFaceColor','#009E73','LineWidth', 3/4,'CapSize',0);
end

if i_n == 3
    Indep_Variable = simData.Initial_Angle - 0.0;
    plt2 = errorbar(Indep_Variable,Dep_Variable,std,'o-','color','#D55E00','MarkerFaceColor','#D55E00','LineWidth', 3/4,'CapSize',0);
end

ylim([0.1*10^-3,1*10^-3])
xlim([0,pi/2])

xticks(0*pi:(1/12)*pi:(1/2)*pi)
xticklabels({'0\pi','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});

end

ESL = [FixedParam.SelectionBias(2), FixedParam.SelectionBias(1)]';
hold on
yline(atan(abs(ESL(2)/ESL(1))), '--','LineWidth',1.5, 'color',greycolor);

%% Aesthetics

figure(1)

xticks(0*pi:(1/12)*pi:(1/2)*pi)
ax = gca;
ax.TickLabelInterpreter = 'latex';
xticklabels({'$0\pi$','$\pi/12$','$\pi/6$','$\pi/4$','$\pi/3$', '$5\pi/12$', '$\pi/2$'});

yticks(0*pi:(1/12)*pi:(1/2)*pi)
ax.TickLabelInterpreter = 'latex';
yticklabels({'','$\pi/12$','$\pi/6$','$\pi/4$','$\pi/3$', '$5\pi/12$', '$\pi/2$'});

H=gca;
H.LineWidth=0.8; %change to the desired value 

xlabel('Initial Angles')
ylabel('Final Angles')

set(gca, 'FontName', 'Helvetica', 'FontSize',12)

figure(2)
xlim([0,pi/2])
ylim([4.3*10^-3, 0.0113])

xticks(0*pi:(1/12)*pi:(1/2)*pi)
ax = gca;
ax.TickLabelInterpreter = 'latex';
xticklabels({'$0\pi$','$\pi/12$','$\pi/6$','$\pi/4$','$\pi/3$', '$5\pi/12$', '$\pi/2$'});

H=gca;
H.LineWidth=0.8; %change to the desired value 

xlabel('Initial Angles')
ylabel('dF/dt')

set(gca, 'FontName', 'Helvetica', 'FontSize',12)