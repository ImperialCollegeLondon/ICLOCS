% main_CarParking - Main script to solve the Optimal Control Problem 
%
% Minimum Time Parallel Parking
% The problem was adapted from 
% B. Li, K. Wang, and Z. Shao, "Time-optimal maneuver planning in automatic parallel parking using a simultaneous dynamic optimization approach". IEEE Transactions on Intelligent Transportation Systems, 17(11), pp.3263-3274, 2016.
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%--------------------------------------------------------

clear all;close all;format compact;

[problem,guess]=CarParking;          % Fetch the problem definition
options= problem.settings(120);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01,[5;2]);

%%
xx=linspace(solution.T(1,1),solution.tf,200);
posx=speval(solution,'X',1,xx);
posy=speval(solution,'X',2,xx);
theta=speval(solution,'X',4,xx);


A_x=posx+(problem.data.auxdata.l_axes+problem.data.auxdata.l_front).*cos(theta)-problem.data.auxdata.b_width.*sin(theta);
B_x=posx+(problem.data.auxdata.l_axes+problem.data.auxdata.l_front).*cos(theta)+problem.data.auxdata.b_width.*sin(theta);
C_x=posx-problem.data.auxdata.l_rear.*cos(theta)+problem.data.auxdata.b_width.*sin(theta);
D_x=posx-problem.data.auxdata.l_rear.*cos(theta)-problem.data.auxdata.b_width.*sin(theta);

A_y=posy+(problem.data.auxdata.l_axes+problem.data.auxdata.l_front).*sin(theta)+problem.data.auxdata.b_width.*cos(theta);
B_y=posy+(problem.data.auxdata.l_axes+problem.data.auxdata.l_front).*sin(theta)-problem.data.auxdata.b_width.*cos(theta);
C_y=posy-problem.data.auxdata.l_rear.*sin(theta)-problem.data.auxdata.b_width.*cos(theta);
D_y=posy-problem.data.auxdata.l_rear.*sin(theta)+problem.data.auxdata.b_width.*cos(theta);


%%
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
xwall=linspace(-8,10,10000);
hold on
plot([-8,10],[problem.data.auxdata.CL problem.data.auxdata.CL],'r-')
pw=problem.data.penalty.values(end);
plot(xwall,(-tanh(pw*(xwall-5/pw))+tanh(pw*(xwall-problem.data.auxdata.SL+5/pw)))/2.*problem.data.auxdata.SW,'r-')
for i=1:1:size(A_x,1)
    plot([A_x(i) B_x(i) C_x(i) D_x(i) A_x(i)], [A_y(i) B_y(i) C_y(i) D_y(i) A_y(i)], 'b-')
end
plot(speval(solution,'X',1,xx),speval(solution,'X',2,xx),'g-','linewidth', 2)
plot(solution.X(:,1),solution.X(:,2),'go','linewidth', 2)
plot(xv(:,1),xv(:,2),'k-.')
xlim([-8 10])
xlabel('Position x [m]')
ylabel('Position y [m]')
grid on


figure
hold on
plot(xx,speval(solution,'X',1,xx),'b-' )
plot(tv,xv(:,1),'k-.')
ylabel('Position x [m]')
xlabel('Time [s]')
grid on


figure
hold on
plot(xx,speval(solution,'X',2,xx),'b-' )
plot(tv,xv(:,2),'k-.')
ylabel('Position y [m]')
xlabel('Time [s]')
grid on


%%
figure
hold on
plot(xx,speval(solution,'X',3,xx),'b-' )
plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)],'r-' )
xlim([0 solution.tf])
plot(tv,xv(:,3),'k-.')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on
% 
figure
hold on
plot(xx,speval(solution,'U',1,xx),'b-' )
plot(tv,uv(:,1),'k-.')
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
grid on

% 
figure
hold on
plot(xx,speval(solution,'X',4,xx)*180/pi,'b-' )
xlim([0 solution.tf])
plot(tv,xv(:,4)*180/pi,'k-.')
xlabel('Time [s]')
ylabel('Orientation Angle [deg]')
grid on

figure
hold on
plot(xx,speval(solution,'X',5,xx)*180/pi,'b-' )
plot(tv,uv(:,2)*180/pi,'k-.')
plot([solution.T(1,1); solution.tf],[problem.states.xl(5), problem.states.xl(5)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(5), problem.states.xu(5)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Steering Angle [deg]')
grid on

figure
if length(solution.T)==length(solution.dU)
    plot(solution.T,solution.dU(:,1)*180/pi,'b-')
else
    plot(solution.T(1:end-1,1),solution.dU(:,1)*180/pi,'b-')
end
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.url(1), problem.inputs.url(1)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uru(1), problem.inputs.uru(1)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Jerk [m/s^3]')
grid on

figure
plot(solution.T(:,1),speval(solution,'U',2,solution.T)*180/pi,'b-')
hold on
% plot(tv,uv(:,2)*180/pi,'k-.')
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(2), problem.inputs.ul(2)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(2), problem.inputs.uu(2)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Control Input (steering rate) [deg/s]')
grid on

% %% figure
% xx=linspace(solution.T(1,1),solution.tf,100000);
% 
% figure
% plot(speval(solution,'X',2,xx)/100,speval(solution,'X',1,xx),'b-')
% hold on
% plot(xv(:,2)/100,xv(:,1),'k-.')
% xlabel('Airspeed [100 m/s]')
% ylabel('Altitude [m]')
% grid on
% 
% figure
% % plot([solution.T(:,1); solution.tf],speval(solution.Xp,3,[solution.T(:,1); solution.tf])*180/pi,'bo' )
% hold on
% plot(xx,speval(solution,'X',3,xx)*180/pi,'b-' )
% plot(tv,xv(:,3)*180/pi,'k-.')
% plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)]*180/pi,'r-' )
% plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)]*180/pi,'r-' )
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Flight Path Angle [deg]')
% grid on
% 
% figure
% % plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,[solution.T(:,1); solution.tf]),'bo' )
% hold on
% plot(xx,speval(solution,'X',1,xx),'b-' )
% plot(tv,xv(:,1),'k-.')
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Altitude [m]')
% grid on
% 
% figure
% % plot([solution.T(:,1); solution.tf],speval(solution.Xp,2,[solution.T(:,1); solution.tf])/100,'bo' )
% hold on
% plot(xx,speval(solution,'X',2,xx)/100,'b-' )
% plot(tv,xv(:,2)/100,'k-.')
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Velocity [100 m/s]')
% grid on
% 
% figure
% % plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,[solution.T(:,1); solution.tf]),'bo' )
% hold on
% plot(xx,speval(solution,'X',4,xx),'b-' )
% plot(tv,xv(:,4),'k-.')
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Aircraft Mass [kg]')
% grid on
% 
% figure
% % plot(solution.T(:,1),speval(solution.Up,1,solution.T)*180/pi,'bo')
% hold on
% plot(xx,speval(solution,'U',1,xx),'b-' )
% plot(tv,uv(:,1),'k-.')
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Control Input (angle of attack) [deg]')
% grid on
% 
