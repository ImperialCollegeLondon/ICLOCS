

figure
plot(X(:,2)/100,X(:,1),'b-')
xlabel('Airspeed [100 m/s]')
ylabel('Altitude [m]')
grid on

figure
plot(T,X(:,3)*180/pi,'b-' )
hold on
plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)]*180/pi,'r-' )
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on

figure
plot(T,X(:,1),'b-' )
hold on
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

figure
plot(T,X(:,2)/100,'b-' )
hold on
xlabel('Time [s]')
ylabel('Velocity [100 m/s]')
grid on

figure
plot(T,X(:,4),'b-' )
hold on
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on

figure
plot(T,U(:,1)*180/pi,'b-')
hold on
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
grid on