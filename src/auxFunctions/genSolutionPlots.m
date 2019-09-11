function genSolutionPlots(options, solution)
%genSolutionPlots - generate plots for the obtained solution
%
% Syntax:  [ X,Xr,U,Ur,P,x0,xf,u0,uf,p,data ] = batchScaleBack(X,Xr,U,Ur,P,x0,xf,u0,uf,p,data)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
% 
%------------- BEGIN CODE --------------

if isfield(solution,'mp')
    
    plotid=options.mp.plot;

    for pidx=1:length(solution.phaseSol)

            if (strcmp(options.phaseoptions{pidx}.discretization,'globalLGR')) || (strcmp(options.phaseoptions{pidx}.discretization,'hpLGR'))
               % Figure generation
                if plotid==1 || plotid==2
                    % Plot state trajectory
                        figure(1)
                        for i=1:size(solution.phaseSol{pidx}.X,2)
                            hold on
                            plot([solution.phaseSol{pidx}.T;solution.phaseSol{pidx}.tf],speval(solution.phaseSol{pidx},'X',i,[solution.phaseSol{pidx}.T;solution.phaseSol{pidx}.tf]));
                        end
                        title('Optimal state trajectories')
                        xlabel('Time')
                        grid on; axis tight;

                    % Plot inputs trajectory
                        figure(2)
                        for i=1:size(solution.phaseSol{pidx}.U,2)
                            hold on
                            plot([solution.phaseSol{pidx}.T;solution.phaseSol{pidx}.tf],speval(solution.phaseSol{pidx},'U',i,[solution.phaseSol{pidx}.T;solution.phaseSol{pidx}.tf]));
                        end
                        title('Optimal input sequences')
                        xlabel('Time')
                        grid on; axis tight;
                end

                if (plotid==1 || plotid==3) && (strcmp(options.mp.transcription,'direct_collocation'))
                    % Plot multipliers

                      figure(3)
                      for i=1:size(solution.phaseSol{pidx}.X,2)
                          hold on
                          plot(solution.phaseSol{pidx}.Tdec,solution.phaseSol{pidx}.multipliers.lambda_1toN(:,i))
                      end
                      title('Adjoint variables')
                      xlabel('Time')
                      grid on; axis tight;

                end

                % Plot discretization error, constaint violation and number of active
                % constraints
                if plotid==1 || plotid==4

                      figure(4)
                        hold on
                        plot(solution.phaseSol{pidx}.T,max(solution.phaseSol{pidx}.Error,[],2))
                        xlabel('Time')
                        title('Maximum absolute Local error') 
                        grid on; axis tight; 

                       figure(5)
                        hold on
                        plot(solution.phaseSol{pidx}.T,max(solution.phaseSol{pidx}.ErrorRelative,[],2))
                        xlabel('Time')
                        title('Maximum relative Local error') 
                        grid on; axis tight; 

                        figure(6)
                        hold on
                        plot(solution.phaseSol{pidx}.T_ConstraintError,max(solution.phaseSol{pidx}.ConstraintError,[],2))
                        xlabel('Time')
                        title('Maximum absolute constraint violation') 
                        grid on; axis tight; 
                end

            else


                % Figure generation
                 if plotid==1 || plotid==2
                    % Plot states
                        figure(1)
                        for i=1:size(solution.phaseSol{pidx}.X,2)
                            hold on
                            plot(solution.phaseSol{pidx}.T,speval(solution.phaseSol{pidx},'X',i,solution.phaseSol{pidx}.T))
                        end
                        title('Optimal state trajectories')
                        xlabel('Time')
                        grid on; axis tight;


                    % Plot inputs
                        figure(2)
                        for i=1:size(solution.phaseSol{pidx}.U,2)
                            hold on
                            plot(solution.phaseSol{pidx}.T,speval(solution.phaseSol{pidx},'U',i,solution.phaseSol{pidx}.T))
                        end
                        title('Optimal input sequences')
                        xlabel('Time')
                        grid on; axis tight;
                end

                if (plotid==1 || plotid==3) && (strcmp(options.mp.transcription,'direct_collocation'))
                    % Plot multipliers
                        figure(3)
                        for i=1:size(solution.phaseSol{pidx}.X,2)
                            hold on
                            plot(solution.phaseSol{pidx}.Tdec,solution.phaseSol{pidx}.multipliers.lambda_1toN(:,i))
                        end
                        title('Adjoint variables')
                        xlabel('Time')
                        grid on; axis tight;
                end

                % Plot discretization error, constaint violation and number of active
                % constraints
                if plotid==1 || plotid==4
                      figure(4)
                        hold on
                        plot(solution.phaseSol{pidx}.T_error(2:end),max(solution.phaseSol{pidx}.Error,[],2))
                        xlabel('Time')
                        title('Maximum absolute local error') 
                        grid on; axis tight; 

                        figure(5)
                        hold on
                        plot(solution.phaseSol{pidx}.T_error(2:end),max(solution.phaseSol{pidx}.ErrorRelative,[],2))
                        xlabel('Time')
                        title('Maximum relative local error') 
                        grid on; axis tight; 

                        figure(6)
                        hold on
                        plot(solution.phaseSol{pidx}.T_ConstraintError,max(solution.phaseSol{pidx}.ConstraintError,[],2))
                        xlabel('Time')
                        title('Maximum absolute constraint violation') 
                        grid on; axis tight; 
                end

            end


        
    end
    
    
else
    
    plotid=options.plot;

    if (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
       % Figure generation
        if plotid==1 || plotid==2
            % Plot state trajectory
                figure(1)
                for i=1:size(solution.X,2)
                    hold on
                    plot([solution.T;solution.tf],speval(solution,'X',i,[solution.T;solution.tf]));
                end
                title('Optimal state trajectories')
                xlabel('Time')
                grid on; axis tight;

            % Plot inputs trajectory
                figure(2)
                for i=1:size(solution.U,2)
                    hold on
                    plot([solution.T;solution.tf],speval(solution,'U',i,[solution.T;solution.tf]));
                end
                title('Optimal input sequences')
                xlabel('Time')
                grid on; axis tight;
        end

        if (plotid==1 || plotid==3) && (strcmp(options.transcription,'direct_collocation'))
            % Plot multipliers

              figure(3)
              for i=1:size(solution.X,2)
                  hold on
                  plot(solution.Tdec,solution.multipliers.lambda_1toN(:,i))
              end
              title('Adjoint variables')
              xlabel('Time')
              grid on; axis tight;

        end

        % Plot discretization error, constaint violation and number of active
        % constraints
        if plotid==1 || plotid==4

              figure(4)
                hold on
                plot(solution.T,max(solution.Error,[],2))
                xlabel('Time')
                title('Maximum absolute Local error') 
                grid on; axis tight; 

               figure(5)
                hold on
                plot(solution.T,max(solution.ErrorRelative,[],2))
                xlabel('Time')
                title('Maximum relative Local error') 
                grid on; axis tight; 

                figure(6)
                hold on
                plot(solution.T_ConstraintError,max(solution.ConstraintError,[],2))
                xlabel('Time')
                title('Maximum absolute constraint violation') 
                grid on; axis tight; 
        end

    else


        % Figure generation
         if plotid==1 || plotid==2
            % Plot states
                figure(1)
                for i=1:size(solution.X,2)
                    hold on
                    plot(solution.T,speval(solution,'X',i,solution.T))
                end
                title('Optimal state trajectories')
                xlabel('Time')
                grid on; axis tight;


            % Plot inputs
                figure(2)
                for i=1:size(solution.U,2)
                    hold on
                    plot(solution.T,speval(solution,'U',i,solution.T))
                end
                title('Optimal input sequences')
                xlabel('Time')
                grid on; axis tight;
        end

        if (plotid==1 || plotid==3) && (strcmp(options.transcription,'direct_collocation'))
            % Plot multipliers
                figure(3)
                for i=1:size(solution.X,2)
                    hold on
                    plot(solution.Tdec,solution.multipliers.lambda_1toN(:,i))
                end
                title('Adjoint variables')
                xlabel('Time')
                grid on; axis tight;
        end

        % Plot discretization error, constaint violation and number of active
        % constraints
        if plotid==1 || plotid==4
              figure(4)
                hold on
                plot(solution.T_error(2:end),max(solution.Error,[],2))
                xlabel('Time')
                title('Maximum absolute local error') 
                grid on; axis tight; 

                figure(5)
                hold on
                plot(solution.T_error(2:end),max(solution.ErrorRelative,[],2))
                xlabel('Time')
                title('Maximum relative local error') 
                grid on; axis tight; 

                figure(6)
                hold on
                plot(solution.T_ConstraintError,max(solution.ConstraintError,[],2))
                xlabel('Time')
                title('Maximum absolute constraint violation') 
                grid on; axis tight; 
        end

    end
end
