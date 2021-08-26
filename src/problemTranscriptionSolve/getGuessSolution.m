function [ x_guess, u_guess, data ] = getGuessSolution( options, guess, problem,  data )
%getGuessSolution - generate initial guess for state and input variables,
%as well as multipliers
%
% Syntax:  [ x_guess, u_guess, data ] = getGuessSolution( options, guess, problem,  data )
%
% Inputs:
%    options, guess, problem,  data - Defined in transcribeOCP
%
% Outputs:
%    x_guess,u_guess - initial guess for state and input variables
%    data - data structure, updated with the guess values of the
%    multipliers
%     
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
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

    
if strcmp(options.discretization,'globalLGR') || strcmp(options.discretization,'hpLGR')
    [nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:17});
    nrc=nrcl+nrcu+nrce;
    x_guess=zeros(M+1,n);u_guess=zeros(M,m);
    if (strcmp(options.start,'Cold'))
        if isfield(guess,'time') && ~isempty(guess.time)
            Tx=linspace(guess.time(1),guess.time(end), M+1);
            Tu=linspace(guess.time(1),guess.time(end), M);
            x_guess=interp1(guess.time, guess.states,Tx,'linear','extrap');
            x_guess=x_guess(:);
            u_guess=interp1(guess.time, guess.inputs,Tu,'linear','extrap');
            u_guess=u_guess(:);
        else
            for i=1:n
              x_guess(:,i)=linspace(guess.states(1,i),guess.states(2,i),M+1);
            end
            x_guess=reshape(x_guess,(M+1)*n,1);
            for i=1:m
             if M>1  
               u_guess(:,i)=linspace(guess.inputs(1,i),guess.inputs(2,i),M)';
              else
               u_guess(:,i)=guess.inputs(1,i);
              end
            end
            u_guess=reshape(u_guess,M*m,1);
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Warm'))
        T=(guess.tf-guess.t0)/2*data.tau_inc+(guess.tf+guess.t0)/2;
        if isfield(guess,'Xp')
            if isfield(guess,'TSeg_Bar')
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,[T;guess.tf]);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,T);
                end
            else
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,[T;guess.tf]);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,T);
                end
            end
            if options.scaling
                x_guess=scale_variables(x_guess, problem.states.scales, problem.states.shifts );
                u_guess=scale_variables(u_guess, problem.inputs.scales, problem.inputs.shifts );
            end
            x_guess=x_guess(:);
            u_guess=u_guess(:);
        else
            if size(guess.time,1)==size(guess.states,1)
                x_guess=interp1(guess.time, guess.states,[T;guess.tf],'pchip');
                x_guess=x_guess(:);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=u_guess(:);
            else
                x_guess=interp1([guess.time;guess.tf], guess.states,[T;guess.tf],'pchip');
                x_guess=x_guess(:);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=u_guess(:);
            end
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Hot'))
        T=(guess.tf-guess.t0)/2*data.tau_inc+(guess.tf+guess.t0)/2;
        if isfield(guess,'Xp')
            if isfield(guess,'TSeg_Bar')
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,[T;guess.tf]);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,T);
                end
            else
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,[T;guess.tf]);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,T);
                end
            end
            if options.scaling
                x_guess=scale_variables(x_guess, problem.states.scales, problem.states.shifts );
                u_guess=scale_variables(u_guess, problem.inputs.scales, problem.inputs.shifts );
            end
            x_guess=x_guess(:);
            u_guess=u_guess(:);
        else
            if size(guess.time,1)==size(guess.states,1)
                x_guess=interp1(guess.time, guess.states,[T;guess.tf],'pchip');
                x_guess=x_guess(:);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=u_guess(:);
            else
                x_guess=interp1([guess.time;guess.tf], guess.states,[T;guess.tf],'pchip');
                x_guess=x_guess(:);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=u_guess(:);
            end
        end
        if isfield(guess.multipliers,'lambda') && ~isempty(guess.multipliers.lambda) && size(guess.multipliers.lambda,2)==1 && size(guess.multipliers.lambda,1)==M*n+M*ng+nb+nrc
            data.multipliers.lambda=guess.multipliers.lambda;
        else
            if isfield(guess.multipliers,'lambda') && ~isempty(guess.multipliers.lambda)
                if isfield(guess,'time_org')
                    data.multipliers.lambda=interp1(guess.time_org, guess.multipliers.lambda,T,'linear','extrap');
                else
                    data.multipliers.lambda=interp1(guess.time, guess.multipliers.lambda,T,'linear','extrap');
                end
                data.multipliers.lambda=data.multipliers.lambda(:);
            else
                data.multipliers.lambda=[];
            end

            if ng
                lambda_g=interp1(guess.timeFull, guess.multipliers.lambda_g,T,'linear','extrap');
                data.multipliers.lambda=[data.multipliers.lambda;lambda_g(data.gAllidx)];
            end
            if nb
                data.multipliers.lambda=[data.multipliers.lambda;guess.multipliers.lambda_b(1:nb)];
            end
            if nrc
                data.multipliers.lambda=[data.multipliers.lambda;zeros(nrc,1)];
            end
            if isfield(guess.multipliers,'lambda_resconst') %&& strcmp(options.transcription,'direct_collocation')
%                 data.multipliers.lambda=[0;data.multipliers.lambda;guess.multipliers.lambda_resconst];
                data.multipliers.lambda=[data.multipliers.lambda;0;guess.multipliers.lambda_resconst];
            end
        end
    end
else
    [nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:13});
    nrc=nrcl+nrcu+nrce;
    x_guess=zeros(M,n);u_guess=zeros(N,m);
    if (strcmp(options.start,'Cold'))
        if isfield(guess,'time') && ~isempty(guess.time)
            Tx=linspace(guess.time(1),guess.time(end), M);
            Tu=linspace(guess.time(1),guess.time(end-1), N);
            x_guess=interp1(guess.time, guess.states,Tx,'linear','extrap');
            x_guess=reshape(x_guess',M*n,1);
            u_guess=interp1(guess.time, guess.inputs,Tu,'linear','extrap');
            u_guess=reshape(u_guess',N*m,1);
        else
            for i=1:n
              x_guess(:,i)=linspace(guess.states(1,i),guess.states(2,i),M);
            end
            for i=1:m
             if N>1  
               u_guess(:,i)=linspace(guess.inputs(1,i),guess.inputs(2,i),N);
              else
               u_guess(:,i)=guess.inputs(1,i);
              end
            end
            x_guess=reshape(x_guess',M*n,1);
            u_guess=reshape(u_guess',N*m,1);
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Warm'))
        T=(guess.tf-guess.t0)*[0;cumsum(data.tau)]*data.Nm/ns+guess.t0;
        if isfield(guess,'Xp')
            if isfield(guess,'TSeg_Bar')
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,T);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,T);
                end
            else
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,T);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,T);
                end
            end
            if options.scaling
                x_guess=scale_variables(x_guess, problem.states.scales, problem.states.shifts );
                u_guess=scale_variables(u_guess, problem.inputs.scales, problem.inputs.shifts );
            end
            x_guess=reshape(x_guess',M*n,1);
            u_guess=reshape(u_guess',M*m,1);
        else
            if size(guess.time,1)==size(guess.states,1)
                x_guess=interp1(guess.time, guess.states,T,'pchip');
                x_guess=reshape(x_guess',M*n,1);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=reshape(u_guess',M*m,1);
            else
                x_guess=interp1([guess.time;guess.tf], guess.states,T,'pchip');
                x_guess=reshape(x_guess',M*n,1);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=reshape(u_guess',M*m,1);
            end
        end
        data.multipliers=[];
    elseif (strcmp(options.start,'Hot'))
        T=(guess.tf-guess.t0)*[0;cumsum(data.tau)]*data.Nm/ns+guess.t0;
        if isfield(guess,'Xp')
            if isfield(guess,'TSeg_Bar')
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,T);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,guess.TSeg_Bar,T);
                end
            else
                for i=1:n
                    x_guess(:,i)=speval(guess,'X',i,T);
                end
                for i=1:m
                    u_guess(:,i)=speval(guess,'U',i,T);
                end
            end
            if options.scaling
                x_guess=scale_variables(x_guess, problem.states.scales, problem.states.shifts );
                u_guess=scale_variables(u_guess, problem.inputs.scales, problem.inputs.shifts );
            end
            x_guess=reshape(x_guess',M*n,1);
            u_guess=reshape(u_guess',M*m,1);
        else
            if size(guess.time,1)==size(guess.states,1)
                x_guess=interp1(guess.time, guess.states,T,'pchip');
                x_guess=reshape(x_guess',M*n,1);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=reshape(u_guess',M*m,1);
            else
                x_guess=interp1([guess.time;guess.tf], guess.states,T,'pchip');
                x_guess=reshape(x_guess',M*n,1);
                u_guess=interp1(guess.time, guess.inputs,T,'pchip');
                u_guess=reshape(u_guess',M*m,1);
            end
        end
        if isfield(guess.multipliers,'lambda') && ~isempty(guess.multipliers.lambda) && size(guess.multipliers.lambda,2)==1 && size(guess.multipliers.lambda,1)==M*n+M*ng+nb+nrc
            data.multipliers.lambda=guess.multipliers.lambda;
        else
            if isfield(guess.multipliers,'lambda') && ~isempty(guess.multipliers.lambda)
                if isfield(guess,'time_org')
                    data.multipliers.lambda=interp1(guess.time_org, guess.multipliers.lambda,T,'linear','extrap');
                else
                    data.multipliers.lambda=interp1(guess.time, guess.multipliers.lambda,T,'linear','extrap');
                end
                data.multipliers.lambda=reshape(data.multipliers.lambda',M*n,1);
            else
                data.multipliers.lambda=[];
            end
            if ng 
                if ~isempty(guess.multipliers.lambda_g)
                    lambda_g=interp1(guess.timeFull, guess.multipliers.lambda_g,T,'linear','extrap');
                else
                    lambda_g=interp1(guess.timeFull, zeros(length(guess.multipliers.lambda_g),1),T,'linear','extrap');
                end
                lambda_g=reshape(lambda_g',M*ng,1);
                data.multipliers.lambda=[data.multipliers.lambda;lambda_g(data.gAllidx)];
            end
            if nb
                data.multipliers.lambda=[data.multipliers.lambda;guess.multipliers.lambda_b(1:nb)];
            end
            if nrc
                data.multipliers.lambda=[data.multipliers.lambda;zeros(nrc,1)];
            end
            if isfield(guess.multipliers,'lambda_resconst') %&& strcmp(options.transcription,'direct_collocation')
%                 data.multipliers.lambda=[0;data.multipliers.lambda;guess.multipliers.lambda_resconst];
                data.multipliers.lambda=[data.multipliers.lambda;0;guess.multipliers.lambda_resconst];
            end 
        end
        
    end
end

end

