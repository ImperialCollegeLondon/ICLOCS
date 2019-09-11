% Main script to solve the Optimal Control Problem
%
% Supersonic Aircraft Minimum Time-to-climb 
%
% The problem was adapted from the supersonic aircraft minimum time-to-climb problem originally presented by
% A. E. Bryson, M. N. Desai, and W. C. Hoffman, "Energy-State Approximation in Performance Optimization of Supersonic Aircraft," Journal of Aircraft, Vol. 6, No. 6, November-December, 1969, pp. 481-488. 
% This example was formulated originally by:
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
% With aerodynamic data and modifications to thrust data by:
% M.A. Patterson and A.V. Rao, User's Manual, "GPOPS-II: A General Purpose MATLAB Software for Solving Multiple-Phase Optimal Control Problems, 2.3 edition", 2016
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

[problem,guess]=LowThrustOrbitTransfer;          % Fetch the problem definition
options= problem.settings(150);                  % h method
% options= problem.settings(100,4);                  % hp method
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.1 );

%% figure
xx=linspace(solution.T(1,1),solution.tf,10000);
% 
figure
p=speval(solution,'X',1,xx);
plot(xx/3600,p/1e06,'b-')
hold on
plot(tv/3600,xv(:,1)/1e06,'k-.')
xlabel('Time [hrs]')
ylabel('p [10^6 ft]')
grid on

figure
f=speval(solution,'X',2,xx);
plot(xx/3600,f,'b-')
hold on
plot(tv/3600,xv(:,2),'k-.')
xlabel('Time [hrs]')
ylabel('f [-]')
grid on

figure
g=speval(solution,'X',3,xx);
plot(xx/3600,g,'b-')
hold on
plot(tv/3600,xv(:,3),'k-.')
xlabel('Time [hrs]')
ylabel('g [-]')
grid on

figure
h=speval(solution,'X',4,xx);
plot(xx/3600,h,'b-')
hold on
plot(tv/3600,xv(:,4),'k-.')
xlabel('Time [hrs]')
ylabel('h [-]')
grid on

figure
k=speval(solution,'X',5,xx);
plot(xx/3600,k,'b-')
hold on
plot(tv/3600,xv(:,5),'k-.')
xlabel('Time [hrs]')
ylabel('k [-]')
grid on

figure
L=speval(solution,'X',6,xx);
plot(xx/3600,L,'b-')
hold on
plot(tv/3600,xv(:,6),'k-.')
xlabel('Time [hrs]')
ylabel('L [rev]')
grid on

figure
u1=speval(solution,'U',1,xx);
plot(xx/3600,u1,'b-')
hold on
plot(tv/3600,uv(:,1),'k-.')
xlabel('Time [hrs]')
ylabel('u_{radial} [-]')
grid on

figure
u2=speval(solution,'U',2,xx);
plot(xx/3600,u2,'b-')
hold on
plot(tv/3600,uv(:,2),'k-.')
xlabel('Time [hrs]')
ylabel('u_{tangential} [-]')
grid on

figure
u3=speval(solution,'U',3,xx);
plot(xx/3600,u3,'b-')
hold on
plot(tv/3600,uv(:,3),'k-.')
xlabel('Time [hrs]')
ylabel('u_{normal} [-]')
grid on


%%
q=1+f.*cos(L)+g.*sin(L);
r=p./q;
alpha_sq=h.^2-k.^2;
chi_sq=h.^2+k.^2;
s_sq=1+chi_sq;

r1=r./s_sq.*(cos(L)+alpha_sq.*cos(L)+2.*h.*k.*sin(L));
r2=r./s_sq.*(sin(L)-alpha_sq.*sin(L)+2.*h.*k.*cos(L));
r3=2*r./s_sq.*(h.*sin(L)-k.*cos(L));

PlotEarth;
plot3(r1,r2,r3,'linewidth',2)