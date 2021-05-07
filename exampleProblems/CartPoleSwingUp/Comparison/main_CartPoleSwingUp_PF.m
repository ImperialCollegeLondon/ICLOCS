% CartPoleSwingUp
%--------------------------------------------------------

clear all;close all;format compact;


%% prerun

[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
solveMyProblem( problem,guess,options);


%% Direct Collcation
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.coll.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll.solution=solution;
sol.coll.sim.tv=tv;
sol.coll.sim.xv=xv;
sol.coll.sim.uv=uv;


%% DAIR with no early termination
tic
[problem,guess]=CartPoleSwingUp_resmin;          % Fetch the problem definition
options= problem.settings(8);
options.resminEarlyStop=0;
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode45', 0.01 );
genSolutionPlots(options, solution);
sol.resmin.solution=solution;
sol.resmin.sim.tv=tv;
sol.resmin.sim.xv=xv;
sol.resmin.sim.uv=uv;
toc

%% Direct Collocation with Mesh Refinement 
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=sqrt(sol.resmin.solution.residuals.r)';
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr.solution=solution;
sol.coll_mr.sim.tv=tv;
sol.coll_mr.sim.xv=xv;
sol.coll_mr.sim.uv=uv;
sol.coll_mr.MRHistory=MRHistory;
sol.coll_mr.timeAll=toc;

%% DAIR Solution point 1
tic
[problem,guess]=CartPoleSwingUp_resmin;          % Fetch the problem definition
problem.states.xErrorTol_integral=[0.03 0.03 0.03 0.03];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_01.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_01.solution=solution;
sol.resmin_01.sim.tv=tv;
sol.resmin_01.sim.xv=xv;
sol.resmin_01.sim.uv=uv;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[0.03 0.03 0.03 0.03];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_01.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_01.solution=solution;
sol.coll_mr_01.sim.tv=tv;
sol.coll_mr_01.sim.xv=xv;
sol.coll_mr_01.sim.uv=uv;
sol.coll_mr_01.MRHistory=MRHistory;
sol.coll_mr_01.timeAll=toc;



%% DAIR Solution point 2
tic
[problem,guess]=CartPoleSwingUp_resmin;          % Fetch the problem definition
problem.states.xErrorTol_integral=[0.05 0.05 0.05 0.05];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_02.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_02.solution=solution;
sol.resmin_02.sim.tv=tv;
sol.resmin_02.sim.xv=xv;
sol.resmin_02.sim.uv=uv;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[0.05 0.05 0.05 0.05];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_02.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_02.solution=solution;
sol.coll_mr_02.sim.tv=tv;
sol.coll_mr_02.sim.xv=xv;
sol.coll_mr_02.sim.uv=uv;
sol.coll_mr_02.MRHistory=MRHistory;
sol.coll_mr_02.timeAll=toc;


%% DAIR Solution point 3
tic
[problem,guess]=CartPoleSwingUp_resmin;          % Fetch the problem definition
problem.states.xErrorTol_integral=[0.1 0.1 0.1 0.1];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_03.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_03.solution=solution;
sol.resmin_03.sim.tv=tv;
sol.resmin_03.sim.xv=xv;
sol.resmin_03.sim.uv=uv;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[0.1 0.1 0.1 0.1];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_03.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_03.solution=solution;
sol.coll_mr_03.sim.tv=tv;
sol.coll_mr_03.sim.xv=xv;
sol.coll_mr_03.sim.uv=uv;
sol.coll_mr_03.MRHistory=MRHistory;
sol.coll_mr_03.timeAll=toc;


%% DAIR Solution point 4
tic
[problem,guess]=CartPoleSwingUp_resmin_costonly;          % Fetch the problem definition
problem.states.xErrorTol_integral=[0.45 0.45 0.45 0.45];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_04.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_04.solution=solution;
sol.resmin_04.sim.tv=tv;
sol.resmin_04.sim.xv=xv;
sol.resmin_04.sim.uv=uv;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[0.45 0.45 0.45 0.45];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_04.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_04.solution=solution;
sol.coll_mr_04.sim.tv=tv;
sol.coll_mr_04.sim.xv=xv;
sol.coll_mr_04.sim.uv=uv;
sol.coll_mr_04.MRHistory=MRHistory;
sol.coll_mr_04.timeAll=toc;


%% DAIR Solution point 5
tic
[problem,guess]=CartPoleSwingUp_resmin_costonly;          % Fetch the problem definition
problem.states.xErrorTol_integral=[0.65 0.65 0.65 0.65];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_05.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_05.solution=solution;
sol.resmin_05.sim.tv=tv;
sol.resmin_05.sim.xv=xv;
sol.resmin_05.sim.uv=uv;
sol.resmin_05.timeAll=toc;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[0.65 0.65 0.65 0.65];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_05.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_05.solution=solution;
sol.coll_mr_05.sim.tv=tv;
sol.coll_mr_05.sim.xv=xv;
sol.coll_mr_05.sim.uv=uv;
sol.coll_mr_05.MRHistory=MRHistory;
sol.coll_mr_05.timeAll=toc;


%% DAIR Solution point 6
tic
[problem,guess]=CartPoleSwingUp_resmin_costonly;          % Fetch the problem definition
problem.states.xErrorTol_integral=[0.8 0.8 0.8 0.8];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_06.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_06.solution=solution;
sol.resmin_06.sim.tv=tv;
sol.resmin_06.sim.xv=xv;
sol.resmin_06.sim.uv=uv;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[0.8 0.8 0.8 0.8];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_06.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_06.solution=solution;
sol.coll_mr_06.sim.tv=tv;
sol.coll_mr_06.sim.xv=xv;
sol.coll_mr_06.sim.uv=uv;
sol.coll_mr_06.MRHistory=MRHistory;
sol.coll_mr_06.timeAll=toc;


%% DAIR Solution point 7
tic
[problem,guess]=CartPoleSwingUp_resmin_costonly;          % Fetch the problem definition
problem.states.xErrorTol_integral=[1 1 1 1];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_07.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_07.solution=solution;
sol.resmin_07.sim.tv=tv;
sol.resmin_07.sim.xv=xv;
sol.resmin_07.sim.uv=uv;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[1 1 1 1];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_07.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_07.solution=solution;
sol.coll_mr_07.sim.tv=tv;
sol.coll_mr_07.sim.xv=xv;
sol.coll_mr_07.sim.uv=uv;
sol.coll_mr_07.MRHistory=MRHistory;
sol.coll_mr_07.timeAll=toc;


%% DAIR Solution point 8
tic
[problem,guess]=CartPoleSwingUp_resmin_costonly;          % Fetch the problem definition
problem.states.xErrorTol_integral=[1.2 1.2 1.2 1.2];
options= problem.settings(8);
[solution,~]=solveMyProblem( problem,guess,options);
sol.resmin_08.timeAll=toc;
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
genSolutionPlots(options, solution);
sol.resmin_08.solution=solution;
sol.resmin_08.sim.tv=tv;
sol.resmin_08.sim.xv=xv;
sol.resmin_08.sim.uv=uv;
sol.resmin_08.timeAll=toc;

%% Direct Collocation with Mesh Refinement
tic
[problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
options= problem.settings(8);
options.meshstrategy='mesh_refinement';
problem.states.xErrorTol_integral=[1.2 1.2 1.2 1.2];
[solution,MRHistory]=solveMyProblem( problem,guess,options);
sol.coll_mr_08.timeAll=toc;
genSolutionPlots(options, solution);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
sol.coll_mr_08.solution=solution;
sol.coll_mr_08.sim.tv=tv;
sol.coll_mr_08.sim.xv=xv;
sol.coll_mr_08.sim.uv=uv;
sol.coll_mr_08.MRHistory=MRHistory;
sol.coll_mr_08.timeAll=toc;

%% figures

% Figure 1
xx=linspace(solution.T(1,1),solution.tf,100000);
fg=figure
set(fg,'color','w');
set(fg,'Position',[10 10 600 600])
subplot(2,2,1)
hold on
p1=plot(xx,speval(sol.coll.solution,'X',3,xx),'linewidth',2,'color',[0.6350, 0.0780, 0.1840])
p2=plot(xx,speval(sol.resmin.solution,'X',3,xx),'linewidth',2,'color',[0, 0.4470, 0.7410])
plot(sol.coll.sim.tv,sol.coll.sim.xv(:,3),'k--','linewidth',2,'color',[0.6350, 0.0780, 0.1840])
plot(sol.resmin.sim.tv,sol.resmin.sim.xv(:,3),'k--','linewidth',2,'color',[0, 0.4470, 0.7410])
xlim([0 solution.tf])
ylabel('$y_1$ [m]','Interpreter','latex', 'FontSize',13);
legend([p1,p2],{'direct collocation with direct interpolation ($J=54.2605$)','DAIR with no early termination ($J=170.1225$)'},'location','northoutside','Interpreter','latex', 'FontSize',11);
grid on


subplot(2,2,2)
hold on
plot(xx,speval(sol.coll.solution,'X',4,xx),'linewidth',2,'color',[0.6350, 0.0780, 0.1840])
plot(xx,speval(sol.resmin.solution,'X',4,xx),'linewidth',2,'color',[0, 0.4470, 0.7410])
plot(sol.coll.sim.tv,sol.coll.sim.xv(:,4),'k--','linewidth',2,'color',[0.6350, 0.0780, 0.1840])
plot(sol.resmin.sim.tv,sol.resmin.sim.xv(:,4),'k--','linewidth',2,'color',[0, 0.4470, 0.7410])
xlim([0 solution.tf])
ylabel('$\theta_1$ [rad]','Interpreter','latex', 'FontSize',13);
grid on



subplot(2,2,3)
hold on
plot(xx,speval(sol.coll.solution,'U',1,xx),'linewidth',2,'color',[0.6350, 0.0780, 0.1840])
plot(xx,speval(sol.resmin.solution,'U',1,xx),'linewidth',2,'color',[0, 0.4470, 0.7410])
plot(sol.coll.sim.tv,sol.coll.sim.uv(:,1),'k--','linewidth',2,'color',[0.6350, 0.0780, 0.1840])
plot(sol.resmin.sim.tv,sol.resmin.sim.uv(:,1),'k--','linewidth',2,'color',[0, 0.4470, 0.7410])
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]','Interpreter','latex', 'FontSize',13);
grid on
ylabel('$u$ [N]','Interpreter','latex', 'FontSize',13);


subplot(2,2,4)
hold on
stairs(sol.coll.solution.T_error,[max(sol.coll.solution.Error,[],2);0],'linewidth',2,'color',[0.6350, 0.0780, 0.1840])
stairs(sol.resmin.solution.T_error,[max(sol.resmin.solution.Error,[],2);0],'linewidth',2,'color',[0, 0.4470, 0.7410])
xlim([0,2])
xlabel('Time [s]','Interpreter','latex', 'FontSize',13);
grid on
ylabel('Error $\eta$','Interpreter','latex', 'FontSize',13);

% Figure 2
fg2=figure
ini_cost=66.667;
ini_res=4.077;
subplot(1,3,3)
hold on
grid on
p1=scatter(ini_res/2,ini_cost,'x','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'linewidth',2);
ylim([-5 180])
subplot(1,3,[1 2])
hold on
grid on
p1=scatter(ini_res/2,ini_cost,'x','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'linewidth',2);
p4=scatter(sum(sol.coll.solution.residuals.r)/2,sol.coll.solution.cost.J ,norm([abs(sol.coll.sim.xv(end,1)),abs(sol.coll.sim.xv(end,2)),abs(sol.coll.sim.xv(end,3)-1),abs(sol.coll.sim.xv(end,4)-pi)*2])*20,'o','MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840]);

p3=scatter(sum(sol.resmin_01.solution.residuals.r)/2,sol.resmin_01.solution.cost.J ,norm([abs(sol.resmin_01.sim.xv(end,1)),abs(sol.resmin_01.sim.xv(end,2)),abs(sol.resmin_01.sim.xv(end,3)-1),abs(sol.resmin_01.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);
scatter(sum(sol.resmin_02.solution.residuals.r)/2,sol.resmin_02.solution.cost.J ,norm([abs(sol.resmin_02.sim.xv(end,1)),abs(sol.resmin_02.sim.xv(end,2)),abs(sol.resmin_02.sim.xv(end,3)-1),abs(sol.resmin_02.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);
scatter(sum(sol.resmin_03.solution.residuals.r)/2,sol.resmin_03.solution.cost.J ,norm([abs(sol.resmin_03.sim.xv(end,1)),abs(sol.resmin_03.sim.xv(end,2)),abs(sol.resmin_03.sim.xv(end,3)-1),abs(sol.resmin_03.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);
scatter(sum(sol.resmin_04.solution.residuals.r)/2,sol.resmin_04.solution.cost.J ,norm([abs(sol.resmin_04.sim.xv(end,1)),abs(sol.resmin_04.sim.xv(end,2)),abs(sol.resmin_04.sim.xv(end,3)-1),abs(sol.resmin_04.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);
scatter(sum(sol.resmin_05.solution.residuals.r)/2,sol.resmin_05.solution.cost.J ,norm([abs(sol.resmin_05.sim.xv(end,1)),abs(sol.resmin_05.sim.xv(end,2)),abs(sol.resmin_05.sim.xv(end,3)-1),abs(sol.resmin_05.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);
scatter(sum(sol.resmin_06.solution.residuals.r)/2,sol.resmin_06.solution.cost.J ,norm([abs(sol.resmin_06.sim.xv(end,1)),abs(sol.resmin_06.sim.xv(end,2)),abs(sol.resmin_06.sim.xv(end,3)-1),abs(sol.resmin_06.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);
scatter(sum(sol.resmin_07.solution.residuals.r)/2,sol.resmin_07.solution.cost.J ,norm([abs(sol.resmin_07.sim.xv(end,1)),abs(sol.resmin_07.sim.xv(end,2)),abs(sol.resmin_07.sim.xv(end,3)-1),abs(sol.resmin_07.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);
scatter(sum(sol.resmin_08.solution.residuals.r)/2,sol.resmin_08.solution.cost.J ,norm([abs(sol.resmin_08.sim.xv(end,1)),abs(sol.resmin_08.sim.xv(end,2)),abs(sol.resmin_08.sim.xv(end,3)-1),abs(sol.resmin_08.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880] ,'linewidth',2);

text(sum(sol.resmin.solution.residuals.r)/2+0.01,sol.resmin.solution.cost.J+5,[num2str(round(sol.resmin.timeAll./sol.coll_mr.timeAll*100-100)) '%']);
text(sum(sol.resmin_01.solution.residuals.r)/2+0.01,sol.resmin_01.solution.cost.J+5,[num2str(round(sol.resmin_01.timeAll./sol.coll_mr_01.timeAll*100-100)) '%']);
text(sum(sol.resmin_04.solution.residuals.r)/2+0.01,sol.resmin_04.solution.cost.J+5,[num2str(round(sol.resmin_04.timeAll./sol.coll_mr_04.timeAll*100-100)) '%']);
text(sum(sol.resmin_07.solution.residuals.r)/2+0.01,sol.resmin_07.solution.cost.J+5,[num2str(round(sol.resmin_07.timeAll./sol.coll_mr_07.timeAll*100-100)) '%']);
xlim([0 0.3])
ylim([-5 180])

p2=scatter(sum(sol.resmin.solution.residuals.r)/2,sol.resmin.solution.cost.J ,norm([abs(sol.resmin.sim.xv(end,1)),abs(sol.resmin.sim.xv(end,2)),abs(sol.resmin.sim.xv(end,3)-1),abs(sol.resmin.sim.xv(end,4)-pi)])*20,'o','MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
xlabel('MIRNS Error','Interpreter','latex', 'FontSize',13);
ylabel('Objective Value','Interpreter','latex', 'FontSize',13);
legend([p1,p4,p2,p3],{'Initial guess','Direct collocation with direct interpolation','DAIR with no early termination','DAIR solutions for different requested error levels'},'location','northoutside','Interpreter','latex', 'FontSize',11);

