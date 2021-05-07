% Goddard Rocket Problem
%--------------------------------------------------------

clear all;close all;format compact;

%% Direct Collocation
[problem,guess]=GoddardRocket_coll;          % Fetch the problem definition
options= problem.settings(100);                  % Get options and solver settings 
[solution,~]=solveMyProblem( problem,guess,options);
sol.coll.solution=solution;

%% Direct Collocation with IRM/DAIR Solution Representation
[problem,guess]=GoddardRocket_DAIRrep;          % Fetch the problem definition
options= problem.settings(100);                  % Get options and solver settings 
[solution,~]=solveMyProblem( problem,guess,options);
sol.DAIRrep.solution=solution;

%% DAIR
[problem,guess]=GoddardRocket_DAIR;          % Fetch the problem definition
options= problem.settings(100);                  % Get options and solver settings 
[solution,~]=solveMyProblem( problem,guess,options);
sol.DAIR.solution=solution;

%% Multi-phase Direct Collocation
options.mp= settings_GoddardRocketThreePhase;                  % Get options and solver settings 
[problem,guess,options.phaseoptions]=GoddardRocketThreePhase;          % Fetch the problem definition
[solution,~]=solveMyProblem( problem,guess,options);
sol.coll_multiphase.solution=solution;

%% Figures
figure
hold on
xlabel('Time [s]','Interpreter','latex', 'FontSize',13);
grid on
ylabel('Thrust [lbf]','Interpreter','latex', 'FontSize',13);
    


xx=linspace(sol.coll.solution.t0,sol.coll.solution.tf,1000);
p2=plot(xx,speval(sol.coll.solution,'U',1,xx),'color',[0.976, 0.780, 0.745],'linewidth',2);

xx=linspace(sol.DAIRrep.solution.t0,sol.DAIRrep.solution.tf,1000);
p3=plot(xx,speval(sol.DAIRrep.solution,'U',1,xx),'color',[0.9290, 0.6940, 0.1250]	,'linewidth',2);

xx=linspace(sol.DAIR.solution.t0,sol.DAIR.solution.tf,1000);
p4=plot(xx,speval(sol.DAIR.solution,'U',1,xx),'color',[0, 0.4470, 0.7410],'linewidth',2);

for i=1:length(sol.coll_multiphase.solution.phaseSol)
    solp=sol.coll_multiphase.solution.phaseSol{i};
    xx=linspace(solp.t0,solp.tf,1000);
    p1=plot(xx,speval(solp,'U',1,xx),'linewidth',2,'color',[0 0 0],'linewidth',2);
    ylim([problem.phases{3}.inputs.ul problem.phases{1}.inputs.uu])
end

legend([p1,p2,p3,p4],{'Multi-phase direct collocation','Single-phase direct collocation','Single-phase direct collocation with solution representation by IRM','Single-phase AIRM with no early termination'},'location','northoutside','Interpreter','latex', 'FontSize',11);
