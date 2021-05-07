% High-index DAE System

%--------------------------------------------------------

clear all;close all;format compact;

%% Direct Collocation
[problem,guess]=HighOrderDAE_coll;          % Fetch the problem definition
options= problem.settings(8,5);                  % Get options and solver settings 
[solution,~]=solveMyProblem( problem,guess,options);
sol.coll.solution=solution;

%% DAIR
[problem,guess]=HighOrderDAE_DAIR;          % Fetch the problem definition
options= problem.settings(8,5);                  % Get options and solver settings 
[solution,~]=solveMyProblem( problem,guess,options);
sol.DAIR.solution=solution;

%% Figure
xx_coll=linspace(sol.coll.solution.T(1,1),sol.coll.solution.T(end,1),1000);
xx_DAIR=linspace(sol.DAIR.solution.T(1,1),sol.DAIR.solution.T(end,1),1000);

figure(1)
hold on
plot(speval(sol.coll.solution,'X',1,xx_coll),speval(sol.coll.solution,'X',3,xx_coll),'color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
plot(speval(sol.DAIR.solution,'X',1,xx_DAIR),speval(sol.DAIR.solution,'X',3,xx_DAIR),'color',[0, 0.4470, 0.7410],'LineWidth',2);

xlabel('State $x_1$ [m]','Interpreter','latex', 'FontSize',13);
ylabel('State $x_3$ [m]','Interpreter','latex', 'FontSize',13);
xlim([-1.3 1.1])
ylim([-1 1])
legend('Direct collocation','DAIR','Interpreter','latex', 'FontSize',11);
grid on
