function [solution]=output(problem,solution,options,data,plotid)
%output - Format and display the solution for h collocation methods. Estimate the discretization error and constraint violation.
%
% Syntax:  output(problem,solution,options,data,plotid)
%
% Inputs:
%    problem - Optimal control problem definition
%    solutions - Structure containing the solution
%    options - Display options
%    data - structure containing matrices to format data 
%    plotid - flag for figure generation
%       * 0: No plotting
%       * 1: Plot all figures
%       * 2: Plot the optimal state trajectory only
%       * 3: Plot the multipliers
%       * 4: Plot the discretization error and constraint violation only
%
% Subfunctions: estimateError, estimateConstraintViolation, calcNumActiveConstraint 
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

% if isfield(problem.inputs,'singular_arc_lift')
%     [ problem ] = avoidSingularArc( problem );
% end


if (strcmp(options.transcription,'globalLGR')) || (strcmp(options.transcription,'hpLGR'))
    % Define some variables
    [nt,np,n,m,~,~,M,~,ns,npd,~,npduidx,~,~,~,~,~]=deal(data.sizes{:});
    vdat=data.data;
    t0=data.t0;

    z=solution.z;
    f=data.functions_unscaled{3};

    % Generate time vector
    if nt
        tf=z(end);
    else
        tf=data.tf(1);
    end

    T=solution.T;
    Tdec=T(1:ns:end);

    % Time segment barrier
    if data.options.adaptseg==1 
        TSeg_Bar=solution.z((end-nt+1):end);
    else
        TSeg_Bar=(tf-t0)/2.*data.tau_segment'+(tf+t0)/2; %Time at start/end of each segmen
    end


    % Extract design parameters
    if np
        p=solution.p;
        P=repmat(p',M,1);
    else
        P=[];
        p=[];
    end

    % Extract states and controls

    X_Np1=solution.X;
    X=X_Np1(1:M,:);
    U=solution.U;
    U_Np1=U;
    U_Np1(end+1,:)=U_Np1(end,:);

    % Define piecewise polynomials
    F=f(X,U,P,T,vdat);
    F_Np1=f(X_Np1,U_Np1,P,[T;solution.tf],vdat);
    Xp=cell(length(npd),n);dXp=cell(length(npd),n);Up=cell(1,m);

    if strcmp(data.options.resultRep,'direct') % Barycentric Lagrange Interpolation
        for i=1:n 
             for j=1:length(npd)
                 idxst=sum(npd(1:j-1))+1;
                 idxed=sum(npd(1:j))+1;
                 Xp{j,i}=[[data.map.LGR.points{npduidx(j)};1],X_Np1(idxst:idxed,i)];
                 dXp{j,i}=[[data.map.LGR.points{npduidx(j)};1],F_Np1(idxst:idxed,i)];
             end
        end
        for i=1:m 
         for j=1:length(npd)
             idxst=sum(npd(1:j-1))+1;
             idxed=sum(npd(1:j))+1;
             Up{j,i}=[[data.map.LGR.points{npduidx(j)};1],U_Np1(idxst:idxed,i)];
         end
        end
    else
        % State trajectory
        if strcmp(data.options.stateRep,'pchip') || strcmp(data.options.resultRep,'default')% pchip interpolation
            for i=1:n
                Xp{i}=pchip([T;solution.tf],X_Np1(:,i)');
                dXp{i}=pchip([T;solution.tf],F_Np1(:,i)');
            end
        elseif strcmp(data.options.stateRep,'Legendre') % Legendre interpolation
            for i=1:n 
                 for j=1:length(npd)
                     idxst=sum(npd(1:j-1))+1;
                     idxed=sum(npd(1:j))+1;
                     Xp{j,i}=legendrefit([data.map.LGR.points{npduidx(j)};1],X_Np1(idxst:idxed,i), npd(j)+1, 'qr');
                     dXp{j,i}=legendrefit(data.map.LGR.points{npduidx(j)},F(idxst:idxed-1,i), npd(j), 'qr');
                 end
            end
        elseif strcmp(data.options.stateRep,'Barycentric') %Barycentric Lagrange Interpolation
            for i=1:n 
                 for j=1:length(npd)
                     idxst=sum(npd(1:j-1))+1;
                     idxed=sum(npd(1:j))+1;
                     Xp{j,i}=[[data.map.LGR.points{npduidx(j)};1],X_Np1(idxst:idxed,i)];
                     dXp{j,i}=[[data.map.LGR.points{npduidx(j)};1],F_Np1(idxst:idxed,i)];
                 end
            end
        else
            error('State representation method invalid or not supported by the selected transcription method');
        end

        % Input trajectory
        if strcmp(data.options.inputRep,'pchip') || strcmp(data.options.resultRep,'default') % pchip interpolation
            for i=1:m 
                Up{i}=pchip([T;solution.tf],U_Np1(:,i)');
            end
        elseif strcmp(data.options.inputRep,'constant') % Piecewise constant polynomials
            for i=1:m 
                Up{i}=mkpp([T;solution.tf],U(:,i)');
            end
        elseif strcmp(data.options.inputRep,'linear') % Piecewise linear polynomials
            for i=1:m 
                Up{i}=Linsplines([T;solution.tf],U_Np1(:,i),diff(U_Np1(:,i)));
            end
        elseif strcmp(data.options.inputRep,'Legendre') % Legendre fitting
            for i=1:m 
             for j=1:length(npd)
                 idxst=sum(npd(1:j-1))+1;
                 idxed=sum(npd(1:j))+1;
                 Up{j,i}=legendrefit([data.map.LGR.points{npduidx(j)};1],U_Np1(idxst:idxed,i), npd(j)+1, 'qr');
             end
            end
        elseif strcmp(data.options.inputRep,'Barycentric') % Barycentric Lagrange Interpolation
            for i=1:m 
             for j=1:length(npd)
                 idxst=sum(npd(1:j-1))+1;
                 idxed=sum(npd(1:j))+1;
                 Up{j,i}=[[data.map.LGR.points{npduidx(j)};1],U_Np1(idxst:idxed,i)];
             end
            end
        else
            error('Input representation method invalid or not supported by the selected transcription method');
        end
    end

    % Display computation time
    if (options.print.time)
        disp('computation time:');disp(solution.computation_time);
    end

    % Display minimized cost
    if (options.print.cost)
        if (strcmp(data.options.transcription,'globalLGR') || strcmp(data.options.transcription,'hpLGR')) && data.options.reorderLGR
            disp('minimized cost:');disp(costFunction(solution.z_org,data));
        else
            disp('minimized cost:');disp(costFunction(z,data));
        end
    end

    % Display discretization error, constaint violation and number of active
    % constraints
        
        [ solution ] = getVariableRate( solution,data );
    
        [Error, ErrorRelative]=estimateError_LGR(Xp,Up,dXp,p,T,TSeg_Bar,n,m,f,M,data);
        
        [ConstraintError,T_ConstraintError,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation_LGR(Xp,Up,p,t0,tf,TSeg_Bar,n,m,solution,problem,data,1);
            
    if options.print.relative_local_error
        disp('Maximum absolute local error:');disp(max(Error));
        disp('Maximum relative local error:');disp(max(ErrorRelative));
        disp('Maximum absolute constraint violation:');disp(max(ConstraintError));
        disp('Number of active constraints:');disp(NumActiveConstraint);
    end

    % Figure generation
    if plotid==1 || plotid==2
        % Plot state trajectory
        if (options.plot.states)
            figure(1)
            for i=1:n
                hold on
                plot([T;solution.tf],ppvalL( Xp,i,TSeg_Bar,[T;solution.tf]));
            end
            title('Optimal state trajectories')
            xlabel('Time')
            grid on; axis tight;
        end

        % Plot inputs trajectory
        if (options.plot.inputs)
            figure(2)
            for i=1:m
                hold on
                plot([T;solution.tf],ppvalL( Up,i,TSeg_Bar,[T;solution.tf]));
            end
            title('Optimal input sequences')
            xlabel('Time')
            grid on; axis tight;
        end
    end

    if plotid==1 || plotid==3
        % Plot multipliers
        if (options.plot.multipliers==1)

          switch(data.options.NLPsolver)
           case{'ipopt'} 
             % Estimate the adjoint variables
               lambda_midpoint=solution.multipliers.lambda;
           case{'worhp'} 
             % Estimate the adjoint variables
               lambda_midpoint=solution.multipliers.lambda;
           case{'fmincon'}  
              lambda_midpoint=reshape(solution.multipliers.eqnonlin(n+1:n*M),n,M-1)';
          end    
          lambda_midpoint=lambda_midpoint(ns:ns:M-1,:);

          if M<=2
              lambda=lambda_midpoint;
          else
              lambda=interp1((1.5:(M-1)/ns+1)',lambda_midpoint,(1:(M-1)/ns+1)','linear','extrap');
          end     

          figure(3)
          for i=1:n
              hold on
              plot(Tdec,lambda(:,i))
          end
          title('Adjoint variables')
          xlabel('Time')
          grid on; axis tight;
        end

    end

    % Plot discretization error, constaint violation and number of active
    % constraints
    if plotid==1 || plotid==4
        if options.print.relative_local_error

          figure(4)
            hold on
            plot(T,max(Error,[],2))
            xlabel('Time')
            title('Maximum absolute Local error') 
            grid on; axis tight; 

           figure(5)
            hold on
            plot(T,max(ErrorRelative,[],2))
            xlabel('Time')
            title('Maximum relative Local error') 
            grid on; axis tight; 

            figure(6)
            hold on
            plot(T_ConstraintError,max(ConstraintError,[],2))
            xlabel('Time')
            title('Maximum absolute constraint violation') 
            grid on; axis tight; 
        end
    end

    % Save results
    solution.Xp=Xp;
    solution.dXp=dXp;
    solution.Up=Up;
    solution.Error=Error;
    solution.ErrorRelative=ErrorRelative;
    solution.ConstraintError=ConstraintError;
    solution.T_ConstraintError=T_ConstraintError;
%     if isfield(options,'AutoDirect')
         solution.ActiveConstraint=ActiveConstraint;
         solution.NumActiveConstraint=NumActiveConstraint;
%     end
    solution.TSeg_Bar=TSeg_Bar;


    
else
    % Define some variables
    [nt,np,n,m,~,~,M,N,ns]=deal(data.sizes{1:9});
    vdat=data.data;
    t0=data.t0;

    z=solution.z;
    Vx=data.map.Vx;
    Vu=data.map.Vu;
    f=data.functions_unscaled{3};

    % Generate time vector
    if nt
        tf=solution.tf;
    else
        tf=data.tf(1);
    end

    T=(tf-t0)*[0;cumsum(data.tau)]*data.Nm/ns+data.k0;
    Tdec=T(1:ns:end);

    % Extract design parameters
    if np
        p=solution.p;
        P=repmat(p',M,1);
    else
        P=[];
        p=[];
    end

    % Extract states and controls

    if strcmp(data.options.transcription,'multiple_shooting')
        X=reshape(Vx*z,n,M)';
        usp=reshape(data.map.Vu*z,m,N)';
        U=kron(usp,ones((M-1)/N,1));
        if data.options.scaling
            X=scale_variables_back( X, data.data.Xscale, data.data.Xshift );
            U=scale_variables_back( U, data.data.Uscale, data.data.Ushift );
        end
        if isfield(vdat,'singular_arc_lift')
            U(:,vdat.singular_arc_lift)=(U(:,vdat.singular_arc_lift)+1).^2+vdat.singular_arc_lift_shift(vdat.singular_arc_lift);
        end
        solution.z_orgscale=[tf*ones(nt);p;zeros(n*M+m*N,1)]+data.map.xV*reshape(X',M*n,1)+data.map.uV*reshape(U',m*N,1);
        U=[U;U(end,:)];

    else
        X=reshape(Vx*z,n,M)';
        U=reshape(Vu*z,m,N)';
        if data.options.scaling
            X=scale_variables_back( X, data.data.Xscale, data.data.Xshift );
            U=scale_variables_back( U, data.data.Uscale, data.data.Ushift );
        end
        if isfield(vdat,'singular_arc_lift')
            U(:,vdat.singular_arc_lift)=(U(:,vdat.singular_arc_lift)+1).^2+vdat.singular_arc_lift_shift(vdat.singular_arc_lift);
        end
        solution.z_orgscale=[tf*ones(nt);p;zeros(n*M+m*N,1)]+data.map.xV*reshape(X',M*n,1)+data.map.uV*reshape(U',m*N,1);
    end





    % Define piecewise polynomials
    F=f(X,U,P,T,vdat);
    Xp=cell(n,1);dXp=cell(n,1);Up=cell(m,1);

    % Solution reconstruction
    if strcmp(data.options.resultRep,'default') || strcmp(data.options.resultRep,'direct')
        if strcmp(data.options.transcription,'hermite')
            for i=1:n % Cubic Hermite interpolation
                [Xp{i}, dXp{i}]=Hsplines(T,X(:,i),F(:,i));
            end
            for i=1:m % Piecewise quadratic polynomials
                [Up{i}]=QuadsplinesU(T,U(:,i));
            end
        elseif strcmp(data.options.transcription,'trapezoidal')
            for i=1:n % Quadratic interpolation
                [Xp{i}, dXp{i}]=Quadsplines(T,X(:,i),F(:,i));
            end
            for i=1:m % Piecewise linear polynomials
                Up{i}=Linsplines(T,U(:,i),diff(U(:,i)));
            end
        else
            for i=1:n % Linear interpolation
                [Xp{i}, dXp{i}]=Linsplines(T,X(1:end,i),F(1:end,i));
            end
            for i=1:m % Piecewise constant polynomials
                Up{i}=mkpp(T,U(1:end-1,i)');
            end
        end
    else
        % State trajectory
        if strcmp(data.options.stateRep,'linear') && strcmp(data.options.transcription,'euler')
            for i=1:n % Linear interpolation
                [Xp{i}, dXp{i}]=Linsplines(T,X(1:end,i),F(1:end,i));
            end
        elseif strcmp(data.options.stateRep,'quadratic') && (strcmp(data.options.transcription,'euler') || strcmp(data.options.transcription,'trapezoidal'))
            for i=1:n % Quadratic interpolation
                [Xp{i}, dXp{i}]=Quadsplines(T,X(:,i),F(:,i));
            end
        elseif strcmp(data.options.stateRep,'cubic') && strcmp(data.options.transcription,'hermite')
            for i=1:n % Cubic Hermite interpolation
                [Xp{i}, dXp{i}]=Hsplines(T,X(:,i),F(:,i));
            end
        elseif strcmp(data.options.stateRep,'pchip')
            for i=1:n
                Xp{i}=pchip(T,X(:,i)');
                dXp{i}=pchip(T,F(:,i)');
            end
        else
            error('State representation method invalid or not supported by the selected transcription method');
        end

        % Input trajectory
        if strcmp(data.options.inputRep,'constant')
            for i=1:m % Piecewise constant polynomials
                Up{i}=mkpp(T,U(1:end-1,i)');
            end
        elseif strcmp(data.options.inputRep,'linear')
            for i=1:m % Piecewise linear polynomials
                Up{i}=Linsplines(T,U(:,i),diff(U(:,i)));
            end
        elseif strcmp(data.options.inputRep,'quadratic') && strcmp(data.options.transcription,'hermite')
            for i=1:m % Piecewise quadratic polynomials
                [Up{i}]=QuadsplinesU(T,U(:,i));
            end
        elseif strcmp(data.options.inputRep,'pchip')
            for i=1:m 
                Up{i}=pchip(T,U(:,i)');
            end
        else
            error('Input representation method invalid or not supported by the selected transcription method');
        end
    end

    % Display computation time
    if (options.print.time)
        disp('computation time:');disp(solution.computation_time);
    end

    % Display minimized cost
    if (options.print.cost)
        disp('minimized cost:');disp(costFunction(z,data));
    end

    % Calculate variable rates
    if ~strcmp(data.options.transcription,'multiple_shooting')
        [ solution ] = getVariableRate( solution,data );
    end
    
    % Display discretization error, constaint violation and number of active
    % constraints
    [Error, ErrorRelative,T_error]=estimateError(Xp,Up,dXp,p,tf,n,m,f,M,ns,data);
    disp('Maximum absolute local error:');disp(max(Error));
    disp('Maximum relative local error:');disp(max(ErrorRelative));
    
    
    [ConstraintError,T_ConstraintError,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation(Xp,Up,p,tf,n,m,solution,problem,data,1);
    disp('Maximum absolute constraint violation:');disp(max(ConstraintError));
%     if isfield(options,'AutoDirect')
    disp('Number of active constraints:');disp(NumActiveConstraint);
%     end

    % Obtain multipliers
    switch(data.options.NLPsolver)
    case{'ipopt'} 
     % Estimate the adjoint variables
       lambda_midpoint=reshape(solution.multipliers.lambda(n+1:n*M),n,M-1)';
    case{'worhp'} 
     % Estimate the adjoint variables
       lambda_midpoint=reshape(solution.multipliers.lambda(n+1:n*M),n,M-1)';
    case{'builtinSQP'} 
     % Estimate the adjoint variables
       lambda_midpoint=reshape(solution.multipliers.lambda(n+1:n*M),n,M-1)';
    case{'fmincon'}  
       lambda_midpoint=reshape(solution.multipliers.lambda(n+1:n*M),n,M-1)';
    end    
    lambda_midpoint=lambda_midpoint(ns:ns:M-1,:);
    if M<=2
      lambda=lambda_midpoint;
    else
      lambda=interp1((1.5:(M-1)/ns+1)',lambda_midpoint,(1:(M-1)/ns+1)','linear','extrap');
    end     

    % Figure generation
    if plotid==1 || plotid==2
        % Plot states
        if (options.plot.states)
            figure(1)
            for i=1:n
                hold on
                plot(T,speval(Xp,i,T),'r')

            end
            title('Optimal state trajectories')
            xlabel('Time')
            grid on; axis tight;
        end


        % Plot inputs
        if (options.plot.inputs)
            figure(2)
            for i=1:m
                hold on
                plot(T,speval(Up,i,T),'r')
            end
            title('Optimal input sequences')
            xlabel('Time')
            grid on; axis tight;
        end
    end

    if plotid==1 || plotid==3
        % Plot multipliers
        if (options.plot.multipliers==1);
            figure(3)
            for i=1:n
                hold on
                plot(Tdec,lambda(:,i),'r')
            end
            title('Adjoint variables')
            xlabel('Time')
            grid on; axis tight;
        end
    end

    % Plot discretization error, constaint violation and number of active
    % constraints
    if plotid==1 || plotid==4
        if options.print.relative_local_error

          figure(4)
            hold on
            plot(T_error(2:end),max(Error,[],2))
            xlabel('Time')
            title('Maximum absolute local error') 
            grid on; axis tight; 

            figure(5)
            hold on
            plot(T_error(2:end),max(ErrorRelative,[],2))
            xlabel('Time')
            title('Maximum relative local error') 
            grid on; axis tight; 

            figure(6)
            hold on
            plot(T_ConstraintError,max(ConstraintError,[],2))
            xlabel('Time')
            title('Maximum absolute constraint violation') 
            grid on; axis tight; 
        end
    end

    % Save results
    solution.org.X=X;
    solution.org.U=U;
    solution.org.T=T;
    solution.org.F=F;
    solution.Xp=Xp;
    solution.dXp=dXp;
    solution.Up=Up;
    solution.Error=Error;
    solution.ErrorRelative=ErrorRelative;
    solution.multipliers.lambda_1toN=lambda;
    solution.ConstraintError=ConstraintError;
    solution.T_ConstraintError=T_ConstraintError;
%     if isfield(options,'AutoDirect')
         solution.ActiveConstraint=ActiveConstraint;
         solution.NumActiveConstraint=NumActiveConstraint;
%     end

end
%------------- END OF CODE --------------

end


function [Error,ErrorRelative,T]=estimateError(Xp,Up,dXp,p,tf,n,m,f,M,ns,data)
%ESTIMATEERROR - Estimate the discretization error
%Represent the solution by the chosen method and intergrate dp/dt-f(p) using Romberg quadrature
%
% Syntax:  [Error,ErrorRelative,T]=estimateError(Xp,Up,dXp,p,tf,n,m,f,M,ns,data)
%
% Inputs:
%    Xp  - Structure containing polynomial coefficients for states.
%    Up  - Structure containing polynomial coefficients for inputs.
%    dXp - Structure containing polynomial coefficients for state derivatives.
%    p  - parameter
%    tf - terminal time
%    n - number of state variables
%    m - number of control variables
%    f - function handle
%    M - total number of mesh nodes (scalar)
%    data - structure containing matrices to format data 
%    
% Output:
%    Error - absolute local error
%    ErrorRelative - relative local error
%    T - time vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%------------- BEGIN CODE --------------


q=11;% estimated relative error 10^(-q);
T=(tf-data.t0)*[0;cumsum(data.tau)]*data.Nm/ns+data.k0;
vdat=data.data;

%Computation of scale weights for the discretization error
wi=zeros(n,1);
for i=1:n
     wi(i)=max(max(abs(speval(Xp,i,T)),abs(speval(dXp,i,T))));
end

Error=zeros(M-1,n);         %Pre-allocation
ErrorRelative=zeros(M-1,n);         %Pre-allocation   

for k=1:M-1


a=T(k); b=T(k+1);
h = 2.^((1:q)-1)*(b-a)/2^(q-1);             % These are the intervals used.
k1 = 2.^((q-2):-1:-1)*2+1;                  % Index into the intervals.
tq=a:h(1):b;



Xq=zeros(length(tq),n);
dXq=zeros(length(tq),n);
for i=1:n
    xp=speval(Xp,i,tq);
    Xq(:,i)=xp(:);
    dxp=speval(dXp,i,tq);
    dXq(:,i)=dxp(:); 
end

Uq=zeros(length(tq),m);
for i=1:m
    up=speval(Up,i,tq); 
    Uq(:,i)=up(:);
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p,k1(1),1);
else
  P=[];  
end

F=dXq-f(Xq,Uq,P,tq',vdat);           % Function evaluations.


% computation of the absolute local error on step k

romberg=zeros(n,1);
for kk=1:n
  R = zeros(1,q);                                           % Pre-allocation.
  % Define the starting vector:
  for ii = 1:q
	R(ii) = 0.5*h(ii)*(F(1,kk)+2*...
                       sum(F(k1(end-ii+1):k1(end-ii+1)-1:end-1,kk))+F(end,kk));
  end
% Interpolations:
for jj = 2:q
    jpower = (4^(jj-1)-1);
    for ii = 1:(q-jj+1)
        R(ii) = R(ii)+(R(ii)-R(ii+1))/jpower; 
    end 
end
romberg(kk) = R(1);
end

ErrorRelative(k,:)=abs(transpose(romberg./(wi+1)));
Error(k,:)=abs(romberg');
end
end


function [MaxConstraintError,tg,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation(Xp,Up,p,tf,n,m,solution,problem,data, id)

%estimateConstraintViolation - Estimate the absolute local constraint violation
%
% Syntax:  [MaxConstraintError,tg]=estimateConstraintViolation(Xp,Up,p,tf,n,m,problem,data)
%
% Inputs:
%    Xp  - Structure containing polynomial coefficients for states.
%    Up  - Structure containing polynomial coefficients for inputs.
%    p  - parameter
%    tf - terminal time
%    n - number of state variables
%    m - number of control variables
%    problem - Optimal control problem definition
%    data - structure containing matrices to format data 
%    
% Output:
%    MaxConstraintError - absolute local constraint violation
%    tg - time vector corresponds to the absolute local constraint violation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%------------- BEGIN CODE --------------

ntg=1000;
tg=linspace(0,tf,ntg);

Xg=zeros(ntg,n);

for i=1:n
    Xg(:,i)=speval(Xp,i,tg);
end

Ug=zeros(ntg,m);
for i=1:m
    Ug(:,i)=speval(Up,i,tg); 
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p,ntg,1);
else
    P=[];
end


g=data.functions_unscaled{4};

G=g(Xg,Ug,P,tg',data.data);           % Function evaluations.

Gl=G-problem.constraints.gl;
Gl(Gl>0)=0;

Gu=problem.constraints.gu-G;
Gu(Gu>0)=0;

Xl=Xg-problem.states.xl;
Xl(Xl>0)=0;

Xu=problem.states.xu-Xg;
Xu(Xu>0)=0;

Ul=Ug-problem.inputs.ul;
Ul(Ul>0)=0;

Uu=problem.inputs.uu-Ug;
Uu(Uu>0)=0;

MaxConstraintError=[abs(Gu),abs(Gl),abs(Xu),abs(Xl),abs(Uu),abs(Ul)];


%% 

if id==1
    Gl=G-problem.constraints.gl;
    Glpactive=min(Gl,problem.constraints.gTol)-problem.constraints.gTol;
    Glpactive=Glpactive<0;
    Glabs=min(min(Gl,[],1),problem.constraints.gTol);
    Glactive=Glabs<problem.constraints.gTol;

    Gu=problem.constraints.gu-G;
    Gupactive=min(Gu,problem.constraints.gTol)-problem.constraints.gTol;
    Gupactive=Gupactive<0;
    Guabs=min(min(Gu,[],1),problem.constraints.gTol);
    Guactive=Guabs<problem.constraints.gTol;

    Xl=solution.X-problem.states.xl;
    Xlabs=min(min(Xl,[],1),problem.states.xConstraintTol);
    Xlactive=Xlabs<problem.states.xConstraintTol;

    Xu=problem.states.xu-solution.X;
    Xuabs=min(min(Xu,[],1),problem.states.xConstraintTol);
    Xuactive=Xuabs<problem.states.xConstraintTol;

    Ul=solution.U-problem.inputs.ul;
    Ulabs=min(min(Ul,[],1),problem.inputs.uConstraintTol);
    Ulactive=Ulabs<problem.inputs.uConstraintTol;

    Uu=problem.inputs.uu-solution.U;
    Uuabs=min(min(Uu,[],1),problem.inputs.uConstraintTol);
    Uuactive=Uuabs<problem.inputs.uConstraintTol;

    NumActiveRateConstraint=0;
    if isfield(solution,'dX') && isfield(problem.states,'xrl')
        Xrl=solution.dX-problem.states.xrl;
        Xrlabs=min(min(Xrl,[],1),problem.states.xrConstraintTol);
        Xrlactive=Xrlabs<problem.states.xrConstraintTol;

        Xru=problem.states.xru-solution.dX;
        Xruabs=min(min(Xru,[],1),problem.states.xrConstraintTol);
        Xruactive=Xruabs<problem.states.xrConstraintTol;

        Url=solution.dU-problem.inputs.url;
        Urlabs=min(min(Url,[],1),problem.inputs.urConstraintTol);
        Urlactive=Urlabs<problem.inputs.urConstraintTol;

        Uru=problem.inputs.uru-solution.dU;
        Uruabs=min(min(Uru,[],1),problem.inputs.urConstraintTol);
        Uruactive=Uruabs<problem.inputs.urConstraintTol;

        ActiveConstraint.Xrlactive=Xrlactive;
        ActiveConstraint.Xruactive=Xruactive;
        ActiveConstraint.Urlactive=Urlactive;
        ActiveConstraint.Uruactive=Uruactive;
        
        NumActiveRateConstraint=nnz(Xrlactive)+nnz(Xruactive)+nnz(Urlactive)+nnz(Uruactive);
    end

    ActiveConstraint.Glactive=Glactive;
    ActiveConstraint.Guactive=Guactive;
    ActiveConstraint.Xlactive=Xlactive;
    ActiveConstraint.Xuactive=Xuactive;
    ActiveConstraint.Ulactive=Ulactive;
    ActiveConstraint.Uuactive=Uuactive;
    ActiveConstraint.Glpactive=Glpactive;
    ActiveConstraint.Gupactive=Gupactive;
    
    NumActiveConstraint=nnz(Glactive)+nnz(Guactive)+nnz(Xlactive)+nnz(Xuactive)+nnz(Ulactive)+nnz(Uuactive)+NumActiveRateConstraint;
else
    ActiveConstraint=[];
    NumActiveConstraint=[];
end



end
%------------- END OF CODE --------------
% 
% function [ActiveConstraint,NumActiveConstraint]=calcNumActiveConstraint(solution,problem,data)
% 
% %calcNumActiveConstraint - calculate the number of active constraints
% %
% % Syntax:  NumActiveConstraint=calcNumActiveConstraint(solution,X,problem,data)
% %
% % Inputs:
% %    solutions - Structure containing the solution
% %    problem - Optimal control problem definition
% %    data - structure containing matrices to format data 
% %
% % Output:
% %    NumActiveConstraint - number of active constraints
% %
% % Other m-files required: none
% % Subfunctions: none
% % MAT-files required: none
% %
% %------------- BEGIN CODE --------------
% 
% 
% g=data.functions_unscaled{4};
% 
% G=g(solution.X,solution.U,solution.p,solution.T,data.data);           % Function evaluations.
% 
% Gl=G-problem.constraints.gl;
% Glabs=min(min(Gl,[],1),problem.constraints.gTol);
% Glactive=Glabs<problem.constraints.gTol;
% 
% Gu=problem.constraints.gu-G;
% Guabs=min(min(Gu,[],1),problem.constraints.gTol);
% Guactive=Guabs<problem.constraints.gTol;
% 
% Xl=solution.X-problem.states.xl;
% Xlabs=min(min(Xl,[],1),problem.states.xConstraintTol);
% Xlactive=Xlabs<problem.states.xConstraintTol;
% 
% Xu=problem.states.xu-solution.X;
% Xuabs=min(min(Xu,[],1),problem.states.xConstraintTol);
% Xuactive=Xuabs<problem.states.xConstraintTol;
% 
% Ul=solution.U-problem.inputs.ul;
% Ulabs=min(min(Ul,[],1),problem.inputs.uConstraintTol);
% Ulactive=Ulabs<problem.inputs.uConstraintTol;
% 
% Uu=problem.inputs.uu-solution.U;
% Uuabs=min(min(Uu,[],1),problem.inputs.uConstraintTol);
% Uuactive=Uuabs<problem.inputs.uConstraintTol;
% 
% if isfield(solution,'dX')
%     Xrl=solution.dX-problem.states.xrl;
%     Xrlabs=min(min(Xrl,[],1),problem.states.xrConstraintTol);
%     Xrlactive=Xrlabs<problem.states.xrConstraintTol;
% 
%     Xru=problem.states.xru-solution.dX;
%     Xruabs=min(min(Xru,[],1),problem.states.xrConstraintTol);
%     Xruactive=Xruabs<problem.states.xrConstraintTol;
% 
%     Url=solution.dU-problem.inputs.url;
%     Urlabs=min(min(Url,[],1),problem.inputs.urConstraintTol);
%     Urlactive=Urlabs<problem.inputs.urConstraintTol;
% 
%     Uru=problem.inputs.uru-solution.dU;
%     Uruabs=min(min(Uru,[],1),problem.inputs.urConstraintTol);
%     Uruactive=Uruabs<problem.inputs.urConstraintTol;
%     
%     ActiveConstraint.Xrlactive=Xrlactive;
%     ActiveConstraint.Xruactive=Xruactive;
%     ActiveConstraint.Urlactive=Urlactive;
%     ActiveConstraint.Uruactive=Uruactive;
% end
% 
% ActiveConstraint.Glactive=Glactive;
% ActiveConstraint.Guactive=Guactive;
% ActiveConstraint.Xlactive=Xlactive;
% ActiveConstraint.Xuactive=Xuactive;
% ActiveConstraint.Ulactive=Ulactive;
% ActiveConstraint.Uuactive=Uuactive;
% 
% NumActiveConstraint=nnz(Glactive)+nnz(Guactive)+nnz(Xlactive)+nnz(Xuactive)+nnz(Ulactive)+nnz(Uuactive)+nnz(Xrlactive)+nnz(Xruactive)+nnz(Urlactive)+nnz(Uruactive);
% 
% 
% end
% %------------- END OF CODE --------------

function [Error, ErrorRelative]=estimateError_LGR(Xp,Up,dXp,p,T,TSeg_Bar,n,m,f,M,data)
%ESTIMATEERROR - Estimate the discretization error
%Represent the solution by the chosen method and intergrate dp/dt-f(p) using Romberg quadrature
%
% Syntax:  [Error, ErrorRelative]=estimateError(Xp,Up,dXp,p,T,TSeg_Bar,n,m,f,M,ns,data)
%
% Inputs:
%    Xp  - Structure containing polynomial coefficients for states.
%    Up  - Structure containing polynomial coefficients for inputs.
%    dXp - Structure containing polynomial coefficients for state derivatives.
%    p  - parameter
%    T - time vector
%    TSeg_Bar - time interval barrier
%    n - number of state variables
%    m - number of control variables
%    M - total number of mesh nodes (scalar)
%    data - structure containing matrices to format data 
%    
% Output:
%    Error - absolute local error
%    ErrorRelative - relative local error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%------------- BEGIN CODE --------------


q=11;% estimated relative error 10^(-q);
tf=TSeg_Bar(end);
T=[T;tf];

vdat=data.data;
%Computation of scale weights for the discretization error

wi=zeros(n,1);
% tau = normalizeT( [T;tf], t0,tf );
for i=1:n
     wi(i)=max(max(abs(speval( Xp,i,TSeg_Bar,T)),abs(speval( dXp,i,TSeg_Bar,T))));
end

Error=zeros(M,n);         %Pre-allocation   
ErrorRelative=zeros(M,n);         %Pre-allocation   
for k=1:M


a=T(k); b=T(k+1);
h = 2.^((1:q)-1)*(b-a)/2^(q-1);             % These are the intervals used.
k1 = 2.^((q-2):-1:-1)*2+1;                  % Index into the intervals.
tq=a:h(1):b;

if ~isempty(tq)
    Xq=zeros(length(tq),n);
    dXq=zeros(length(tq),n);
    for i=1:n
        xp=speval( Xp,i,TSeg_Bar,tq);
        Xq(:,i)=xp(:);
        dxp=speval( dXp,i,TSeg_Bar,tq);
        dXq(:,i)=dxp(:); 
    end

    Uq=zeros(length(tq),m);
    for i=1:m
        up=speval( Up,i,TSeg_Bar,tq);
        Uq(:,i)=up(:);
    end

    % Extract design parameters if specified and convert to cells
    if ~isempty(p) 
        P=repmat(p,k1(1),1);
    else
      P=[];  
    end

    F=dXq-f(Xq,Uq,P,tq',vdat);           % Function evaluations.


    % computation of the absolute local error on step k

    romberg=zeros(n,1);
    for kk=1:n
      R = zeros(1,q);                                           % Pre-allocation.
      % Define the starting vector:
      for ii = 1:q
        R(ii) = 0.5*h(ii)*(F(1,kk)+2*...
                           sum(F(k1(end-ii+1):k1(end-ii+1)-1:end-1,kk))+F(end,kk));
      end
    % Interpolations:
    for jj = 2:q
        jpower = (4^(jj-1)-1);
        for ii = 1:(q-jj+1)
            R(ii) = R(ii)+(R(ii)-R(ii+1))/jpower; 
        end 
     end
    romberg(kk) = R(1);
    end

    ErrorRelative(k,:)=abs(transpose(romberg./(wi+1)));
    Error(k,:)=abs(romberg');
else
    ErrorRelative(k,:)=0;
    Error(k,:)=0;
end

end


end

function [MaxConstraintError,tg,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation_LGR(Xp,Up,p,t0,tf,TSeg_Bar,n,m,solution,problem,data,id)
%estimateConstraintViolation - Estimate the absolute local constraint violation
%
% Syntax:  [MaxConstraintError,tg]=estimateConstraintViolation(Xp,Up,p,tf,TSeg_Bar,n,m,problem,data)
%
% Inputs:
%    Xp  - Structure containing polynomial coefficients for states.
%    Up  - Structure containing polynomial coefficients for inputs.
%    p  - parameter
%    tf - terminal time
%    TSeg_Bar - time interval barrier
%    n - number of state variables
%    m - number of control variables
%    problem - Optimal control problem definition
%    data - structure containing matrices to format data 
%    
% Output:
%    MaxConstraintError - absolute local constraint violation
%    tg - time vector corresponds to the absolute local constraint violation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%------------- BEGIN CODE --------------

ntg=1000;
tg=linspace(t0,tf,ntg);

Xg=zeros(ntg,n);

for i=1:n
    Xg(:,i)=speval(Xp,i,TSeg_Bar,tg);
end

Ug=zeros(ntg,m);
for i=1:m
    Ug(:,i)=speval(Up,i,TSeg_Bar,tg);
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p,ntg,1);
else
    P=[];
end


g=data.functions_unscaled{4};

G=g(Xg,Ug,P,tg',data.data);           % Function evaluations.

Gl=G-problem.constraints.gl;
Gl(Gl>0)=0;

Gu=problem.constraints.gu-G;
Gu(Gu>0)=0;

Xl=Xg-problem.states.xl;
Xl(Xl>0)=0;

Xu=problem.states.xu-Xg;
Xu(Xu>0)=0;

Ul=Ug-problem.inputs.ul;
Ul(Ul>0)=0;

Uu=problem.inputs.uu-Ug;
Uu(Uu>0)=0;

MaxConstraintError=[abs(Gu),abs(Gl),abs(Xu),abs(Xl),abs(Uu),abs(Ul)];


%% 

if id==1
    Gl=G-problem.constraints.gl;
    Glpactive=min(Gl,problem.constraints.gTol)-problem.constraints.gTol;
    Glpactive=Glpactive<0;
    Glabs=min(min(Gl,[],1),problem.constraints.gTol);
    Glactive=Glabs<problem.constraints.gTol;

    Gu=problem.constraints.gu-G;
    Gupactive=min(Gu,problem.constraints.gTol)-problem.constraints.gTol;
    Gupactive=Gupactive<0;
    Guabs=min(min(Gu,[],1),problem.constraints.gTol);
    Guactive=Guabs<problem.constraints.gTol;

    Xl=solution.X-problem.states.xl;
    Xlabs=min(min(Xl,[],1),problem.states.xConstraintTol);
    Xlactive=Xlabs<problem.states.xConstraintTol;

    Xu=problem.states.xu-solution.X;
    Xuabs=min(min(Xu,[],1),problem.states.xConstraintTol);
    Xuactive=Xuabs<problem.states.xConstraintTol;

    Ul=solution.U-problem.inputs.ul;
    Ulabs=min(min(Ul,[],1),problem.inputs.uConstraintTol);
    Ulactive=Ulabs<problem.inputs.uConstraintTol;

    Uu=problem.inputs.uu-solution.U;
    Uuabs=min(min(Uu,[],1),problem.inputs.uConstraintTol);
    Uuactive=Uuabs<problem.inputs.uConstraintTol;

    NumActiveRateConstraint=0;
    if isfield(solution,'dX') && isfield(problem.states,'xrl')
        Xrl=solution.dX-problem.states.xrl;
        Xrlabs=min(min(Xrl,[],1),problem.states.xrConstraintTol);
        Xrlactive=Xrlabs<problem.states.xrConstraintTol;

        Xru=problem.states.xru-solution.dX;
        Xruabs=min(min(Xru,[],1),problem.states.xrConstraintTol);
        Xruactive=Xruabs<problem.states.xrConstraintTol;

        Url=solution.dU-problem.inputs.url;
        Urlabs=min(min(Url,[],1),problem.inputs.urConstraintTol);
        Urlactive=Urlabs<problem.inputs.urConstraintTol;

        Uru=problem.inputs.uru-solution.dU;
        Uruabs=min(min(Uru,[],1),problem.inputs.urConstraintTol);
        Uruactive=Uruabs<problem.inputs.urConstraintTol;

        ActiveConstraint.Xrlactive=Xrlactive;
        ActiveConstraint.Xruactive=Xruactive;
        ActiveConstraint.Urlactive=Urlactive;
        ActiveConstraint.Uruactive=Uruactive;
        
        NumActiveRateConstraint=nnz(Xrlactive)+nnz(Xruactive)+nnz(Urlactive)+nnz(Uruactive);
    end

    ActiveConstraint.Glactive=Glactive;
    ActiveConstraint.Guactive=Guactive;
    ActiveConstraint.Xlactive=Xlactive;
    ActiveConstraint.Xuactive=Xuactive;
    ActiveConstraint.Ulactive=Ulactive;
    ActiveConstraint.Uuactive=Uuactive;
    ActiveConstraint.Glpactive=Glpactive;
    ActiveConstraint.Gupactive=Gupactive;

    NumActiveConstraint=nnz(Glactive)+nnz(Guactive)+nnz(Xlactive)+nnz(Xuactive)+nnz(Ulactive)+nnz(Uuactive)+NumActiveRateConstraint;
else
    ActiveConstraint=[];
    NumActiveConstraint=[];
end

end

% function [NumActiveConstraint]=calcNumActiveConstraint_LGR(solution,problem,data)
% %calcNumActiveConstraint - calculate the number of active constraints
% %
% % Syntax:  NumActiveConstraint=calcNumActiveConstraint(solution,X,problem,data)
% %
% % Inputs:
% %    solutions - Structure containing the solution
% %    problem - Optimal control problem definition
% %    data - structure containing matrices to format data 
% %
% % Output:
% %    NumActiveConstraint - number of active constraints
% %
% % Other m-files required: none
% % Subfunctions: none
% % MAT-files required: none
% %
% %------------- BEGIN CODE --------------
% X=solution.X(1:end-1,:);
% 
% g=data.functions_unscaled{4};
% 
% G=g(X,solution.U,solution.p,solution.T,data.data);           % Function evaluations.
% 
% Gl=G-problem.constraints.gl;
% Glabs=min(min(Gl,[],1),0);
% Glactive=length(Glabs(Glabs<1e-04));
% 
% Gu=problem.constraints.gu-G;
% Guabs=min(min(Gu,[],1),0);
% Guactive=length(Guabs(Guabs<1e-04));
% 
% Xl=solution.X-problem.states.xl;
% Xlabs=min(min(Xl,[],1),0);
% Xlactive=length(Xlabs(Xlabs<1e-04));
% 
% Xu=problem.states.xu-solution.X;
% Xuabs=min(min(Xu,[],1),0);
% Xuactive=length(Xuabs(Xuabs<1e-04));
% 
% Ul=solution.U-problem.inputs.ul;
% Ulabs=min(min(Ul,[],1),0);
% Ulactive=length(Ulabs(Ulabs<1e-04));
% 
% Uu=problem.inputs.uu-solution.U;
% Uuabs=min(min(Uu,[],1),0);
% Uuactive=length(Uuabs(Uuabs<1e-04));
% 
% NumActiveConstraint=Glactive+Guactive+Xlactive+Xuactive+Ulactive+Uuactive;
% 
% end
% 
% %------------- END OF CODE --------------