NoTurb_ASAP=load('NoTurb_ASAP.mat');
NoTurb_EBMPC=load('NoTurb_EBMPC.mat');
NoTurb_Fixed=load('NoTurb_Fixed.mat');
Turb_ASAP=load('Turb_ASAP.mat');
Turb_EBMPC=load('Turb_EBMPC.mat');
Turb_Fixed=load('Turb_Fixed.mat');


%%
figure
hold on
plot(NoTurb_Fixed.X(:,2)/100,NoTurb_Fixed.X(:,1),'g-' ,'linewidth',2)
plot(Turb_Fixed.X(:,2)/100,Turb_Fixed.X(:,1),'b-.' ,'linewidth',2)
xlabel('Airspeed [100 m/s]')
ylabel('Altitude [m]')
legend('Closed-loop Simulation (fixed-step 1s update) without turbulence','Closed-loop Simulation (fixed-step 1s update) with turbulence')
grid on


figure
hold on
plot(NoTurb_Fixed.T,NoTurb_Fixed.X(:,3)*180/pi,'g-' ,'linewidth',2)
plot(Turb_Fixed.T,Turb_Fixed.X(:,3)*180/pi,'b-.' ,'linewidth',2)
plot([NoTurb_Fixed.T(1); NoTurb_Fixed.T(end)],[NoTurb_Fixed.problem.states.xl(3), NoTurb_Fixed.problem.states.xl(3)]*180/pi,'r-' )
plot([NoTurb_Fixed.T(1); NoTurb_Fixed.T(end)],[NoTurb_Fixed.problem.states.xu(3), NoTurb_Fixed.problem.states.xu(3)]*180/pi,'r-' )
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
legend('Closed-loop Simulation (fixed-step 1s update) without turbulence','Closed-loop Simulation (fixed-step 1s update) with turbulence')
grid on


figure
hold on
plot(NoTurb_Fixed.T,NoTurb_Fixed.X(:,1),'g-' ,'linewidth',2)
plot(Turb_Fixed.T,Turb_Fixed.X(:,1),'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('Closed-loop Simulation (fixed-step 1s update) without turbulence','Closed-loop Simulation (fixed-step 1s update) with turbulence')
grid on

figure
hold on
plot(NoTurb_Fixed.T,NoTurb_Fixed.X(:,2)/100,'g-' ,'linewidth',2)
plot(Turb_Fixed.T,Turb_Fixed.X(:,2)/100,'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Velocity [100 m/s]')
legend('Closed-loop Simulation (fixed-step 1s update) without turbulence','Closed-loop Simulation (fixed-step 1s update) with turbulence')
grid on

figure
hold on
plot(NoTurb_Fixed.T,NoTurb_Fixed.X(:,4),'g-' ,'linewidth',2)
plot(Turb_Fixed.T,Turb_Fixed.X(:,4),'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
legend('Closed-loop Simulation (fixed-step 1s update) without turbulence','Closed-loop Simulation (fixed-step 1s update) with turbulence')
grid on

figure
hold on
plot(NoTurb_Fixed.T,NoTurb_Fixed.U(:,1)*180/pi,'g-' ,'linewidth',2)
plot(Turb_Fixed.T,Turb_Fixed.U(:,1)*180/pi,'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
legend('Closed-loop Simulation (fixed-step 1s update) without turbulence','Closed-loop Simulation (fixed-step 1s update) with turbulence')
grid on

%%
figure
hold on
plot(NoTurb_ASAP.X(:,2)/100,NoTurb_ASAP.X(:,1),'g-' ,'linewidth',2)
plot(Turb_ASAP.X(:,2)/100,Turb_ASAP.X(:,1),'b-.' ,'linewidth',2)
xlabel('Airspeed [100 m/s]')
ylabel('Altitude [m]')
legend('Closed-loop Simulation (ASAP update) without turbulence','Closed-loop Simulation (ASAP update) with turbulence')
grid on


figure
hold on
plot(NoTurb_ASAP.T,NoTurb_ASAP.X(:,3)*180/pi,'g-' ,'linewidth',2)
plot(Turb_ASAP.T,Turb_ASAP.X(:,3)*180/pi,'b-.' ,'linewidth',2)
plot([NoTurb_ASAP.T(1); NoTurb_ASAP.T(end)],[NoTurb_ASAP.problem.states.xl(3), NoTurb_ASAP.problem.states.xl(3)]*180/pi,'r-' )
plot([NoTurb_ASAP.T(1); NoTurb_ASAP.T(end)],[NoTurb_ASAP.problem.states.xu(3), NoTurb_ASAP.problem.states.xu(3)]*180/pi,'r-' )
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
legend('Closed-loop Simulation (ASAP update) without turbulence','Closed-loop Simulation (ASAP update) with turbulence')
grid on


figure
hold on
plot(NoTurb_ASAP.T,NoTurb_ASAP.X(:,1),'g-' ,'linewidth',2)
plot(Turb_ASAP.T,Turb_ASAP.X(:,1),'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('Closed-loop Simulation (ASAP update) without turbulence','Closed-loop Simulation (ASAP update) with turbulence')
grid on

figure
hold on
plot(NoTurb_ASAP.T,NoTurb_ASAP.X(:,2)/100,'g-' ,'linewidth',2)
plot(Turb_ASAP.T,Turb_ASAP.X(:,2)/100,'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Velocity [100 m/s]')
legend('Closed-loop Simulation (ASAP update) without turbulence','Closed-loop Simulation (ASAP update) with turbulence')
grid on

figure
hold on
plot(NoTurb_ASAP.T,NoTurb_ASAP.X(:,4),'g-' ,'linewidth',2)
plot(Turb_ASAP.T,Turb_ASAP.X(:,4),'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
legend('Closed-loop Simulation (ASAP update) without turbulence','Closed-loop Simulation (ASAP update) with turbulence')
grid on

figure
hold on
plot(NoTurb_ASAP.T,NoTurb_ASAP.U(:,1)*180/pi,'g-' ,'linewidth',2)
plot(Turb_ASAP.T,Turb_ASAP.U(:,1)*180/pi,'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
legend('Closed-loop Simulation (ASAP update) without turbulence','Closed-loop Simulation (ASAP update) with turbulence')
grid on

%%
figure
hold on
plot(NoTurb_EBMPC.X(:,2)/100,NoTurb_EBMPC.X(:,1),'g-' ,'linewidth',2)
plot(Turb_EBMPC.X(:,2)/100,Turb_EBMPC.X(:,1),'b-.' ,'linewidth',2)
xlabel('Airspeed [100 m/s]')
ylabel('Altitude [m]')
legend('Closed-loop Simulation (event based update) without turbulence','Closed-loop Simulation (event based update) with turbulence')
grid on


figure
hold on
plot(NoTurb_EBMPC.T,NoTurb_EBMPC.X(:,3)*180/pi,'g-' ,'linewidth',2)
plot(Turb_EBMPC.T,Turb_EBMPC.X(:,3)*180/pi,'b-.' ,'linewidth',2)
plot([NoTurb_EBMPC.T(1); NoTurb_EBMPC.T(end)],[NoTurb_EBMPC.problem.states.xl(3), NoTurb_EBMPC.problem.states.xl(3)]*180/pi,'r-' )
plot([NoTurb_EBMPC.T(1); NoTurb_EBMPC.T(end)],[NoTurb_EBMPC.problem.states.xu(3), NoTurb_EBMPC.problem.states.xu(3)]*180/pi,'r-' )
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
legend('Closed-loop Simulation (event based update) without turbulence','Closed-loop Simulation (event based update) with turbulence')
grid on


figure
hold on
plot(NoTurb_EBMPC.T,NoTurb_EBMPC.X(:,1),'g-' ,'linewidth',2)
plot(Turb_EBMPC.T,Turb_EBMPC.X(:,1),'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('Closed-loop Simulation (event based update) without turbulence','Closed-loop Simulation (event based update) with turbulence')
grid on

figure
hold on
plot(NoTurb_EBMPC.T,NoTurb_EBMPC.X(:,2)/100,'g-' ,'linewidth',2)
plot(Turb_EBMPC.T,Turb_EBMPC.X(:,2)/100,'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Velocity [100 m/s]')
legend('Closed-loop Simulation (event based update) without turbulence','Closed-loop Simulation (event based update) with turbulence')
grid on

figure
hold on
plot(NoTurb_EBMPC.T,NoTurb_EBMPC.X(:,4),'g-' ,'linewidth',2)
plot(Turb_EBMPC.T,Turb_EBMPC.X(:,4),'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
legend('Closed-loop Simulation (event based update) without turbulence','Closed-loop Simulation (event based update) with turbulence')
grid on

figure
hold on
plot(NoTurb_EBMPC.T,NoTurb_EBMPC.U(:,1)*180/pi,'g-' ,'linewidth',2)
plot(Turb_EBMPC.T,Turb_EBMPC.U(:,1)*180/pi,'b-.' ,'linewidth',2)
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
legend('Closed-loop Simulation (event based update) without turbulence','Closed-loop Simulation (event based update) with turbulence')
grid on

