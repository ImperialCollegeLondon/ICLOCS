initSim;

t_end=20;
x0=problem.states.x0;
simtime=0;
T=[];
X=[];
U=[];
odesolver='ode113';
problem.sim.functions=@MinTimeClimbBryson_Dynamics_Sim_Phase1;

while x0<problem.data.altf-10
    X=[X;x0];
    T=[T;simtime];
    tpx0=[simtime;x0'];
    [ tfpu ] = ICLOCS_NLPSolver(tpx0);
    U=[U;tfpu(3:end)'];
    x0 = simulateDynamics( problem, simtime, tfpu(3:end), x0, tstep, odesolver );
    simtime=simtime+tstep;
end

plotSimSolution;