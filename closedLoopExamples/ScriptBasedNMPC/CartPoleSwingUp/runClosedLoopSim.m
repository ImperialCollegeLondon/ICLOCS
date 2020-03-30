initSim;

t_end=20;
x0=problem.states.x0;
T=[];
X=[];
U=[];
odesolver='ode113';
problem.sim.functions=@CartPoleSwingUp_Dynamics_Sim_Exact;
% problem.sim.functions=@CartPoleSwingUp_Dynamics_Sim_Inexact;

for simtime=0:tstep:6
    X=[X;x0];
    T=[T;simtime];
    tpx0=[simtime;x0'];
    [ tfpu ] = ICLOCS_NLPSolver(tpx0);
    U=[U;tfpu(3:end)'];
    x0 = simulateDynamics( problem, simtime, tfpu(3:end), x0, tstep, odesolver );
end

plotSimSolution;