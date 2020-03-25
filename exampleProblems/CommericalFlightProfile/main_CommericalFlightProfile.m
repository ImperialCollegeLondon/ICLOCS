% Main script to solve the Optimal Control Problem
%
% F50MinTimeFlight - Minimum Time Flight Profile for Commercial Aircraft
%
% The aerodynamic and propulsion data are obtained from
% "Performance model Fokker 50", Delft University of Technology, 2010
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2019
% iclocs@imperial.ac.uk


%--------------------------------------------------------

clear all;close all;format compact;

[problem,guess]=CommericalFlightProfile;          % Fetch the problem definition
options= problem.settings(60);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);

%%
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
center=[problem.data.auxdata.obs_epos_1 problem.data.auxdata.obs_npos_1];
obspos = [center-problem.data.auxdata.obs_r_1 2*problem.data.auxdata.obs_r_1 2*problem.data.auxdata.obs_r_1];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
hold on
center=[problem.data.auxdata.obs_epos_2 problem.data.auxdata.obs_npos_2];
obspos = [center-problem.data.auxdata.obs_r_2 2*problem.data.auxdata.obs_r_2 2*problem.data.auxdata.obs_r_2];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_3 problem.data.auxdata.obs_npos_3];
obspos = [center-problem.data.auxdata.obs_r_3 2*problem.data.auxdata.obs_r_3 2*problem.data.auxdata.obs_r_3];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_4 problem.data.auxdata.obs_npos_4];
obspos = [center-problem.data.auxdata.obs_r_4 2*problem.data.auxdata.obs_r_4 2*problem.data.auxdata.obs_r_4];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_5 problem.data.auxdata.obs_npos_5];
obspos = [center-problem.data.auxdata.obs_r_5 2*problem.data.auxdata.obs_r_5 2*problem.data.auxdata.obs_r_5];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_6 problem.data.auxdata.obs_npos_6];
obspos = [center-problem.data.auxdata.obs_r_6 2*problem.data.auxdata.obs_r_6 2*problem.data.auxdata.obs_r_6];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_7 problem.data.auxdata.obs_npos_7];
obspos = [center-problem.data.auxdata.obs_r_7 2*problem.data.auxdata.obs_r_7 2*problem.data.auxdata.obs_r_7];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_8 problem.data.auxdata.obs_npos_8];
obspos = [center-problem.data.auxdata.obs_r_8 2*problem.data.auxdata.obs_r_8 2*problem.data.auxdata.obs_r_8];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_9 problem.data.auxdata.obs_npos_9];
obspos = [center-problem.data.auxdata.obs_r_9 2*problem.data.auxdata.obs_r_9 2*problem.data.auxdata.obs_r_9];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
center=[problem.data.auxdata.obs_epos_10 problem.data.auxdata.obs_npos_10];
obspos = [center-problem.data.auxdata.obs_r_10 2*problem.data.auxdata.obs_r_10 2*problem.data.auxdata.obs_r_10];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');

plot(speval(solution,'X',3,xx),speval(solution,'X',2,xx),'b-','LineWidth',2)
xlabel('East Position [m]')
ylabel('North Position [m]')
% grid on
plot(problem.states.x0(3)-1000, problem.states.x0(2)-1000, '.k', 'MarkerSize',20)
plot(problem.states.xfl(3)-1000, problem.states.xfl(2)-1000, '.k', 'MarkerSize',20)
text(problem.states.x0(3)+20000,problem.states.x0(2),'ORG')
text(problem.states.xfl(3)+10000,problem.states.xfl(2)+10000,'DES')
xlim([-1 10]*10^5)
ylim([-0.5 9]*10^5)

figure
center=[problem.data.auxdata.obs_epos_1 problem.data.auxdata.obs_npos_1];
obspos = [center-problem.data.auxdata.obs_r_1 2*problem.data.auxdata.obs_r_1 2*problem.data.auxdata.obs_r_1];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
hold on
plot(speval(solution,'X',3,xx),speval(solution,'X',2,xx),'b-','LineWidth',2)
xlabel('East Position [m]')
ylabel('North Position [m]')
grid on
plot(problem.states.x0(3)-1000, problem.states.x0(2)-1000, '.k', 'MarkerSize',20)
plot(problem.states.xfl(3)-1000, problem.states.xfl(2)-1000, '.k', 'MarkerSize',20)
text(problem.states.x0(3)+10000,problem.states.x0(2),'ORG')
text(problem.states.xfl(3)-10000,problem.states.xfl(2)-5000,'DES')
text(8.52e05, 7.8e05,'NO FLIGHT ZONE','Color','white','FontSize',14)
xlim([8.5 9.1]*10^5)
ylim([7.5 8.1]*10^5)

%%
figure
subplot(2,1,1)
hold on
plot([solution.T(1,1); solution.tf],[problem.states.xu(1), problem.states.xu(1)],'r-' )
plot(xx,speval(solution,'X',1,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

subplot(2,1,2)
plot([solution.T(1,1); solution.tf],[problem.states.xl(4), problem.states.xl(4)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(4), problem.states.xu(4)],'r-' )
hold on
plot(xx,speval(solution,'X',4,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('True Airspeed [m/s]')
grid on

figure
subplot(2,1,1)
hold on
plot(xx,speval(solution,'X',5,xx)*180/pi,'b-' )
hold on
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on


subplot(2,1,2)
hold on
plot(xx,speval(solution,'X',6,xx)*180/pi,'b-' )
hold on
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Tracking Angle [deg]')
grid on

figure
subplot(2,1,1)
hold on
plot(xx,speval(solution,'X',7,xx)/9.81,'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on

subplot(2,1,2)
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)]*180/pi,'r-' )
plot(xx,speval(solution,'U',1,xx)*180/pi,'b-' )

xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Angle of attack (Control) [deg]')
grid on


figure
subplot(2,1,1)
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(2), problem.inputs.ul(2)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(2), problem.inputs.uu(2)]*180/pi,'r-' )
hold on
plot(xx,speval(solution,'U',2,xx)*180/pi,'b-' )

xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Roll angle (Control) [deg]')
grid on

subplot(2,1,2)
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(3), problem.inputs.ul(3)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(3), problem.inputs.uu(3)],'r-' )
hold on
plot(xx,speval(solution,'U',3,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Throttle Setting (Control) [-]')
grid on


