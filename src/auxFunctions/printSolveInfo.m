function printSolveInfo(solution,options)
%printSolveInfo - print out the solution information
%
% Syntax:  printSolveInfo(solution,options)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

    if isfield(solution,'mp')
        
        if (options.mp.print.time)
            disp('computation time:');disp(solution.mp.computation_time);
        end
        
        % Display minimized cost
        if (options.mp.print.cost)
            disp('minimized cost:');disp(solution.mp.cost);
        end

        if options.mp.print.relative_local_error
            for i=1:length(solution.phaseSol)
                disp(['Phase' num2str(i) ':']);
                disp('Maximum absolute local error:');disp(solution.phaseSol{i}.MaxAbsError);
                disp('Maximum relative local error:');disp(solution.phaseSol{i}.MaxRelError);
                if isfield(options.mp.print,'residual_error') && options.mp.print.residual_error
                    disp('Integrated Squared ODE residual:');disp(solution.phaseSol{i}.residuals.r');
                end
                disp('Maximum absolute constraint violation:');disp(solution.phaseSol{i}.MaxConstVioError);
                disp('Number of active constraints:');disp(solution.phaseSol{i}.NumActiveConstraint);
            end
        end
        
        
    else
        if (options.print.time)
            disp('computation time:');disp(solution.computation_time);
        end

    
        % Display minimized cost
        if (options.print.cost)
            disp('minimized cost:');disp(solution.cost);
        end

        if options.print.relative_local_error
            disp('Maximum absolute local error:');disp(solution.MaxAbsError);
            disp('Maximum relative local error:');disp(solution.MaxRelError);
                if isfield(options.print,'residual_error') && options.print.residual_error
                    disp('Integrated Squared ODE residual:');disp(solution.residuals.r');
                end
            disp('Maximum absolute constraint violation:');disp(solution.MaxConstVioError);
            disp('Number of active constraints:');disp(solution.NumActiveConstraint);
        end
    end
end

