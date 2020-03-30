data_exact=load('data_exact.mat');
data_inexact=load('data_inexact.mat');


figure
plot(data_exact.T,data_exact.X(:,3),'b-' ,'linewidth',2)
hold on
plot(data_inexact.T,data_inexact.X(:,3),'b--' ,'linewidth',2)
plot([data_exact.tswitch data_exact.tswitch],[min(data_inexact.X(:,3)) max(data_inexact.X(:,3))],'k-')
plot([data_inexact.tswitch data_inexact.tswitch],[min(data_inexact.X(:,3)) max(data_inexact.X(:,3))],'k--')
ylim([min(data_inexact.X(:,3)) max(data_inexact.X(:,3))])
xlabel('Time [s]')
ylabel('Position [m]')
legend('Closed-loop Simulation with exact dynamics','Closed-loop Simulation with modified dynamics','Controller switch (Sim with exact dynamics)','Controller switch time (Sim with modified dynamics)')
grid on

figure
plot(data_exact.T,data_exact.X(:,4),'b-'  ,'linewidth',2)
hold on
plot(data_inexact.T,data_inexact.X(:,4),'b--'  ,'linewidth',2)
plot([data_exact.tswitch data_exact.tswitch],[min(data_inexact.X(:,4)) max(data_inexact.X(:,4))],'k-')
plot([data_inexact.tswitch data_inexact.tswitch],[min(data_inexact.X(:,4)) max(data_inexact.X(:,4))],'k--')
ylim([min(data_inexact.X(:,4)) max(data_inexact.X(:,4))])
xlabel('Time [s]')
ylabel('Angle [rad]')
legend('Closed-loop Simulation with exact dynamics','Closed-loop Simulation with modified dynamics','Controller switch (Sim with exact dynamics)','Controller switch time (Sim with modified dynamics)')
grid on


figure
plot(data_exact.T,data_exact.U(:,1),'b-' ,'linewidth',2)
hold on
plot(data_inexact.T,data_inexact.U(:,1),'b--' ,'linewidth',2)
plot([data_exact.tswitch data_exact.tswitch],[min(data_inexact.U(:,1)) max(data_inexact.U(:,1))],'k-')
plot([data_inexact.tswitch data_inexact.tswitch],[min(data_inexact.U(:,1)) max(data_inexact.U(:,1))],'k--')
ylim([min(data_inexact.U(:,1)) max(data_inexact.U(:,1))])
xlabel('Time [s]')
ylabel('Control Input (Force) [N]')
legend('Closed-loop Simulation with exact dynamics','Closed-loop Simulation with modified dynamics','Controller switch (Sim with exact dynamics)','Controller switch time (Sim with modified dynamics)')
grid on