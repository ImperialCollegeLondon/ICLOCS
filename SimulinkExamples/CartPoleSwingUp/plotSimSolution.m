


figure
plot(T,X(:,3),'b-' ,'linewidth',2)
hold on
plot([tswitch tswitch],[min(X(:,3)) max(X(:,3))],'k--')

ylim([min(X(:,3)) max(X(:,3))])
xlabel('Time [s]')
ylabel('Position [m]')
grid on

figure
plot(T,X(:,4),'b-'  ,'linewidth',2)
hold on
plot([tswitch tswitch],[min(X(:,4)) max(X(:,4))],'k--')
ylim([min(X(:,4)) max(X(:,4))])
xlabel('Time [s]')
ylabel('Angle [rad]')
grid on


figure
plot(T,U(:,1),'b-' ,'linewidth',2)
hold on
plot([tswitch tswitch],[min(U(:,1)) max(U(:,1))],'k--')
ylim([min(U(:,1)) max(U(:,1))])
xlabel('Time [s]')
ylabel('Control Input (Force) [N]')
grid on