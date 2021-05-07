% CartPoleSwingUp
%--------------------------------------------------------

clear all;close all;format compact;

N_vect=[64 32 16 8];
r_vect=[14 5 1 0.1];


for i=1:length(N_vect)
    [problem,guess]=CartPoleSwingUp_directcoll;          % Fetch the problem definition
    options= problem.settings(N_vect(i));                  % Get options and solver settings 
    [solution,~]=solveMyProblem( problem,guess,options);
    genSolutionPlots(options, solution);
    [ tv, xv, uv ] = simulateSolution( problem, solution, 'ode45', 0.01 );
    sol.coll{i}.solution=solution;
    sol.coll{i}.sim.tv=tv;
    sol.coll{i}.sim.xv=xv;
    sol.coll{i}.sim.uv=uv;
    sol.resError.coll(i)=sum(sol.coll{i}.solution.residuals.r);

    [problem,guess]=CartPoleSwingUp_resmin;          % Fetch the problem definition
    options= problem.settings(N_vect(i));                  % Get options and solver settings 
    options.resminEarlyStop=0;
    options.minresRelaxPct=r_vect(i);
    [solution,MRHistory]=solveMyProblem( problem,guess,options);
    [ tv, xv, uv ] = simulateSolution( problem, solution, 'ode45', 0.01 );
    genSolutionPlots(options, solution);
    sol.resmin{i}.solution=solution;
    sol.resmin{i}.sim.tv=tv;
    sol.resmin{i}.sim.xv=xv;
    sol.resmin{i}.sim.uv=uv;
    sol.resError.resmin(i)=sum(sol.resmin{i}.solution.residuals.r);
end

%%
set(groot,'defaultAxesTickLabelInterpreter','latex');  

figure
p1=loglog(2./N_vect,sol.resError.coll/2,'o' ,'LineWidth',2,'color',[0.6350, 0.0780, 0.1840]);
hold on
p2=loglog(2./N_vect,sol.resError.resmin/2,'x' ,'LineWidth',2,'color',[0, 0.4470, 0.7410]);

xlim([min(2./N_vect), max(2./N_vect)])
xticklabels({'$2^{-5}$','$2^{-4}$','$2^{-3}$','$2^{-2}$'});
xticks(2./N_vect)
xlabel('Mesh Interval Size [$s$]','Interpreter','latex', 'FontSize',13);
ylabel('MIRNS Error','Interpreter','latex', 'FontSize',13);
legend([p1,p2],{'Direct collocation with direct interpolation','DAIR with no early termination'},'location','eastoutside','Interpreter','latex', 'FontSize',11);
grid on

