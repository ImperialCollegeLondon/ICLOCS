function [solution]=postSolveAnalysis(problem,solution,options,data)
%output - Format and display the solution for h collocation methods. Estimate the discretization error and constraint violation.
%
% Syntax:  output(problem,solution,options,data,plotid)
%
% Inputs:
%    problem - Optimal control problem definition
%    solutions - Structure containing the solution
%    options - Display options
%    data - structure containing matrices to format data 
%
% Subfunctions: estimateError, estimateConstraintViolation, calcNumActiveConstraint 
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

if isfield(data,'dataNLP')
    dataNLP=data.dataNLP;
else
    dataNLP=data;
end

if isfield(data,'mode') && data.mode==0
    data.mode=1;
end

if (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
    [nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive,ng_eq,ng_neq]=deal(dataNLP.sizes{1:19});
else
    [nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive,nps,ng_eq,ng_neq]=deal(dataNLP.sizes{1:16});
end

if strcmp(options.transcription,'integral_res_min')
    problem.constraints.gTol=[problem.constraints.gTol_neq];
else
    problem.constraints.gl=[zeros(1,ng_eq) problem.constraints.gl];
    problem.constraints.gu=[zeros(1,ng_eq) problem.constraints.gu];
    problem.constraints.gTol=[problem.constraints.gTol_eq problem.constraints.gTol_neq];
end


if (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
    vdat=dataNLP.data;
    z=solution.z;
    
    if (options.print.cost)
        if strcmp(dataNLP.options.transcription,'integral_res_min')
            if (strcmp(dataNLP.options.discretization,'globalLGR') || strcmp(dataNLP.options.discretization,'hpLGR')) && dataNLP.options.reorderLGR
                J=costFunction(solution.z_org,data);
            else
                J=costFunction(z,data);
            end
        else
            if (strcmp(dataNLP.options.discretization,'globalLGR') || strcmp(dataNLP.options.discretization,'hpLGR')) && dataNLP.options.reorderLGR
                J=costFunction(solution.z_org,dataNLP);
            else
                J=costFunction(z,dataNLP);
            end
        end
    end
    solution.cost=J;
    f=dataNLP.functions_unscaled{3};

    % Generate time vector
    if nt==1
        t0=dataNLP.t0;
        tf=z(end);
    elseif nt>=2
        t0=z(end-nt+1);
        tf=z(end);
    else
        tf=dataNLP.tf(1);
    end

    T=solution.T;
    Tdec=T(1:ns:end);
    solution.Tdec=Tdec;
    
    % Time segment barrier
    if dataNLP.options.adaptseg==1 
        TSeg_Bar=solution.z((end-nt+1):end);
    else
        TSeg_Bar=(tf-t0)/2.*dataNLP.tau_segment'+(tf+t0)/2; %Time at start/end of each segmen
    end


    % Extract design parameters
    if np
        p=solution.p;
        P=repmat(p',M,1);
        P_Np1=repmat(p',M+1,1);
    else
        P=[];
        P_Np1=[];
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
    F_Np1=f(X_Np1,U_Np1,P_Np1,[T;tf],vdat);
    Xp=cell(length(npd),n);dXp=cell(length(npd),n);Up=cell(1,m);
    
    if isfield(data,'dataNLP')
        tf=z(end);
        t0=z(end-nt+1);
        
        X_Np1=solution.coll.X; %States including the end point
        X=X_Np1(1:M,:); %States exlduing the end point
        U=solution.coll.U; %input exluding the end point
        U_Np1=[U;U(end,:)];
        
        if np
            p=z(n*(M+1)+M*m+1:n*(M+1)+M*m+np);
        else
            p=[];
        end

        
        X_quad=(data.InterpH*X_Np1)./data.sumInterpH;
        U_quad=(data.InterpH*U_Np1)./data.sumInterpH;
        X_quad(data.interp_fixi,:)=X_Np1(data.interp_fixj,:);
        U_quad(data.interp_fixi,:)=U_Np1(data.interp_fixj,:);

        DT=(data.tf_list-data.t0_list)'*(tf-t0)/2;
        DT_seg_quad=repelem(DT',data.npd_quad+1,1);
        dX_quad=data.D_mat*X_quad./DT_seg_quad;

        LGR_points_quad=[data.LGR.points{1,data.idx_quad};1];
        
        npd_quad_p1=data.npd_quad+1;
        for i=1:n 
             for j=1:length(npd_quad_p1)
                 idxst=sum(npd_quad_p1(1:j-1))+1;
                 idxed=sum(npd_quad_p1(1:j));
                 Xp{j,i}=[LGR_points_quad,X_quad(idxst:idxed,i)];
                 dXp{j,i}=[LGR_points_quad,dX_quad(idxst:idxed,i)];
             end
        end
        for i=1:m 
         for j=1:length(data.npd_quad)
             idxst=sum(npd_quad_p1(1:j-1))+1;
             idxed=sum(npd_quad_p1(1:j));
             Up{j,i}=[LGR_points_quad,U_quad(idxst:idxed,i)];
         end
        end

        solution.tf=tf;
        solution.p=p;
        solution.resmin.X=X_quad;
        solution.resmin.U=U_quad;
        solution.resmin.T=data.tau_quad*(tf-t0);
         
        solution.Xp=Xp;
        solution.dXp=dXp;
        solution.Up=Up; 
        solution.TSeg_Bar=TSeg_Bar;
        solution.x0=speval(solution,'X',1:n,t0);
        solution.X=speval(solution,'X',1:n,[T;tf]);
        solution.U=speval(solution,'U',1:m,T);
        solution.T=T;

        dataNLP.minres_auxdata.D_mat_quad=data.D_mat;
        dataNLP.minres_auxdata.npd_quad=data.npd_quad;
        

        
    else
        if strcmp(dataNLP.options.resultRep,'res_min') % Direct Min-Residual Representation
            [ X_quad, dX_quad, U_quad, t_quad, LGR_points_quad, npd_quad, t0_minres, tf_minres, p_minres, scaled_solution, D_mat_quad ]  = solInterpolationMinResidual( problem,solution,TSeg_Bar,options,dataNLP );
            npd_quad_p1=npd_quad+1;
            TSeg_Bar=TSeg_Bar/TSeg_Bar(end)*tf_minres;
            T=(T-t0_minres)/solution.tf*tf_minres+t0_minres;
            for i=1:n 
                 for j=1:length(npd_quad_p1)
                     idxst=sum(npd_quad_p1(1:j-1))+1;
                     idxed=sum(npd_quad_p1(1:j));
                     Xp{j,i}=[LGR_points_quad,X_quad(idxst:idxed,i)];
                     dXp{j,i}=[LGR_points_quad,dX_quad(idxst:idxed,i)];
                 end
            end
            for i=1:m 
             for j=1:length(npd_quad)
                 idxst=sum(npd_quad_p1(1:j-1))+1;
                 idxed=sum(npd_quad_p1(1:j));
                 Up{j,i}=[LGR_points_quad,U_quad(idxst:idxed,i)];
             end
            end

            if dataNLP.options.scaling
                solution.scaledVariables.resmin.X=scaled_solution.X;
                solution.scaledVariables.resmin.U=scaled_solution.U;
                solution.scaledVariables.X=scaled_solution.X;
                solution.scaledVariables.U=scaled_solution.U(1:end-1,:);
                if isfield(dataNLP.data,'Pscale')
                    solution.scaledVariables.p=scaled_solution.p;
                end
                solution.scaledVariables.x0=scaled_solution.X(1,:);
            end
            solution.coll.tf=solution.tf;
            solution.tf=tf_minres;
            
            solution.Xp=Xp;
            solution.dXp=dXp;
            solution.Up=Up; 
            solution.TSeg_Bar=TSeg_Bar;
            
            tf=tf_minres;
            solution.p=p_minres;
            solution.resmin.X=X_quad;
            solution.resmin.U=U_quad;
            solution.resmin.T=t_quad*(tf-t0)+t0;
            solution.x0=speval(solution,'X',1:n,t0);
            solution.X=speval(solution,'X',1:n,[T;tf]);
            solution.U=speval(solution,'U',1:m,T);
            solution.T=T;

            dataNLP.minres_auxdata.D_mat_quad=D_mat_quad;
            dataNLP.minres_auxdata.npd_quad=npd_quad;

        else
            % State trajectory
            if strcmp(dataNLP.options.resultRep,'default')  || strcmp(dataNLP.options.resultRep,'res_min_final_default') || strcmp(dataNLP.options.stateRep,'Barycentric')%Barycentric Lagrange Interpolation
                for i=1:n 
                     for j=1:length(npd)
                         idxst=sum(npd(1:j-1))+1;
                         idxed=sum(npd(1:j))+1;
                         Xp{j,i}=[[dataNLP.map.LGR.points{npduidx(j)};1],X_Np1(idxst:idxed,i)];
                         dXp{j,i}=[[dataNLP.map.LGR.points{npduidx(j)};1],F_Np1(idxst:idxed,i)];
                     end
                end
            elseif strcmp(dataNLP.options.stateRep,'pchip')  % pchip interpolation
                for i=1:n
                    Xp{i}=pchip([T;solution.tf],X_Np1(:,i)');
                    dXp{i}=pchip([T;solution.tf],F_Np1(:,i)');
                end
            elseif strcmp(dataNLP.options.stateRep,'Legendre') % Legendre interpolation
                for i=1:n 
                     for j=1:length(npd)
                         idxst=sum(npd(1:j-1))+1;
                         idxed=sum(npd(1:j))+1;
                         Xp{j,i}=legendrefit([dataNLP.map.LGR.points{npduidx(j)};1],X_Np1(idxst:idxed,i), npd(j)+1, 'qr');
                         dXp{j,i}=legendrefit(dataNLP.map.LGR.points{npduidx(j)},F(idxst:idxed-1,i), npd(j), 'qr');
                     end
                end

            else
                error('State representation method invalid or not supported by the selected discretization method');
            end

            % Input trajectory
            if strcmp(dataNLP.options.resultRep,'default') || strcmp(dataNLP.options.resultRep,'res_min_final_default') || strcmp(dataNLP.options.inputRep,'Barycentric')  % Barycentric Lagrange Interpolation
                for i=1:m 
                 for j=1:length(npd)
                     idxst=sum(npd(1:j-1))+1;
                     idxed=sum(npd(1:j))+1;
                     Up{j,i}=[[dataNLP.map.LGR.points{npduidx(j)};1],U_Np1(idxst:idxed,i)];
                 end
                end
            elseif strcmp(dataNLP.options.inputRep,'pchip')   % pchip interpolation
                for i=1:m 
                    Up{i}=pchip([T;solution.tf],U_Np1(:,i)');
                end
            elseif strcmp(dataNLP.options.inputRep,'constant') % Piecewise constant polynomials
                for i=1:m 
                    Up{i}=mkpp([T;solution.tf],U(:,i)');
                end
            elseif strcmp(dataNLP.options.inputRep,'linear') % Piecewise linear polynomials
                for i=1:m 
                    Up{i}=Linsplines([T;solution.tf],U_Np1(:,i),diff(U_Np1(:,i)));
                end
            elseif strcmp(dataNLP.options.inputRep,'Legendre') % Legendre fitting
                for i=1:m 
                 for j=1:length(npd)
                     idxst=sum(npd(1:j-1))+1;
                     idxed=sum(npd(1:j))+1;
                     Up{j,i}=legendrefit([dataNLP.map.LGR.points{npduidx(j)};1],U_Np1(idxst:idxed,i), npd(j)+1, 'qr');
                 end
                end
            else
                error('Input representation method invalid or not supported by the selected discretization method');
            end
            
            solution.Xp=Xp;
            solution.dXp=dXp;
            solution.Up=Up; 
            solution.TSeg_Bar=TSeg_Bar;
        end
    end


    % Display discretization error, constaint violation and number of active
    % constraints
    [ solution ] = getVariableRate( solution,dataNLP );
    [Error, ErrorRelative]=estimateError_LGR(solution,p,T,n,m,f,M,dataNLP);
    [ConstraintError,T_ConstraintError,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation_LGR(p,t0,tf,n,m,solution,problem,dataNLP,1);
            
    % Save results
    solution.Error=Error;
    solution.ErrorRelative=ErrorRelative;
    solution.ConstraintError=ConstraintError;
    solution.T_ConstraintError=T_ConstraintError;
    solution.ActiveConstraint=ActiveConstraint;
    solution.NumActiveConstraint=NumActiveConstraint;

    solution.MaxAbsError=max(solution.Error);
    solution.MaxRelError=max(solution.ErrorRelative);
    solution.MaxConstVioError=max(solution.ConstraintError);

    if isfield(dataNLP.options.print,'residual_error') && dataNLP.options.print.residual_error
        [r,r_seg]=estimateResidual_LGR(solution,p,t0,tf,n,m,data);
        solution.residuals.r=r;
        solution.residuals.r_seg=r_seg;
    end

    
else
    % Define some variables
    vdat=dataNLP.data;
    z=solution.z;
    J=costFunction(z,data);
    solution.cost=J;
    Vx=dataNLP.map.Vx;
    Vu=dataNLP.map.Vu;
    f=dataNLP.functions_unscaled{3};

    % Generate time vector
    if nt==1
        tf=solution.tf;
        t0=dataNLP.t0;
    elseif nt==2
        tf=solution.tf;
        t0=solution.t0;
    else
        tf=dataNLP.tf(1);
        t0=dataNLP.t0;
    end

    T=(tf-t0)*[0;cumsum(dataNLP.tau)]*dataNLP.Nm/ns+t0;
    Tdec=T(1:ns:end);
    solution.Tdec=Tdec;
    
    % Extract design parameters
    if np
        p=solution.p;
        P=repmat(p',M,1);
    else
        P=[];
        p=[];
    end

    % Extract states and controls

    if strcmp(dataNLP.options.transcription,'multiple_shooting')
        X=reshape(Vx*z,n,M)';
        usp=reshape(dataNLP.map.Vu*z,m,N)';
        U=kron(usp,ones((M-1)/N,1));
        if dataNLP.options.scaling
            X=scale_variables_back( X, dataNLP.data.Xscale, dataNLP.data.Xshift );
            U=scale_variables_back( U, dataNLP.data.Uscale, dataNLP.data.Ushift );
        end
        solution.z_orgscale=[tf*ones(nt);p;zeros(n*M+m*N,1)]+dataNLP.map.xV*reshape(X',M*n,1)+dataNLP.map.uV*reshape(U',m*N,1);
        U=[U;U(end,:)];

    else
        X=reshape(Vx*z,n,M)';
        U=reshape(Vu*z,m,N)';
        if dataNLP.options.scaling
            X=scale_variables_back( X, dataNLP.data.Xscale, dataNLP.data.Xshift );
            U=scale_variables_back( U, dataNLP.data.Uscale, dataNLP.data.Ushift );
        end
        if nt==1
            solution.z_orgscale=[tf;p;zeros(n*M+m*N,1)]+dataNLP.map.xV*reshape(X',M*n,1)+dataNLP.map.uV*reshape(U',m*N,1);
        elseif nt==2
            solution.z_orgscale=[t0;tf;p;zeros(n*M+m*N,1)]+dataNLP.map.xV*reshape(X',M*n,1)+dataNLP.map.uV*reshape(U',m*N,1);
        else
            solution.z_orgscale=[p;zeros(n*M+m*N,1)]+dataNLP.map.xV*reshape(X',M*n,1)+dataNLP.map.uV*reshape(U',m*N,1);
        end
    end





    
    F=f(X,U,P,T,vdat);
    solution.coll.F=F;
    Fi=repelem(F,2,1);Fi(4:4:end,:)=[];Fi(1,:)=[];Fi(end,:)=[];
    
    % Define piecewise polynomials
    Xp=cell(n,1);dXp=cell(n,1);Up=cell(m,1);

    % Solution reconstruction
    if isfield(data,'dataNLP')
        
        if nt==1
            t0=data.dataNLP.t0;
            tf=z(1);
        elseif nt==2
            t0=z(1);
            tf=z(2);
        else
            t0=data.dataNLP.t0;
            tf=data.dataNLP.tf;
        end
        
        p=solution.p;
        F_k=f(X(1:2:end,:),U(1:2:end,:),repmat(p',length(data.tau),1),data.tau(1:2:end)*(tf-t0),dataNLP.data);
        F_kph=data.DxHS_hf*X/(tf-t0)-F_k(1:end-1,:)/2;
        F_kp1=data.DxHS_p1*X/(tf-t0)+F_k(1:end-1,:);
        F=[F_k(1:end-1,:) F_kph F_kp1]';
        F=reshape(F(:),n,3*data.nps)';

        T=data.tau*(tf-t0)+t0;
        
        for i=1:n % Cubic Hermite interpolation
            [Xp{i}, dXp{i}]=HSInterpolation(T,X(:,i),F(:,i));
        end
        for i=1:m % Piecewise quadratic polynomials
            [Up{i}]=HSInterpolationU(T,U(:,i));
        end       
    else

        if strcmp(dataNLP.options.resultRep,'res_min') % Direct Min-Residual Construction
            [ Xp, dXp, Up, t0_minres, tf_minres, p_minres, minres_solution, scaled_solution ]  = solInterpolationMinResidual( problem,solution,solution.T,options,dataNLP );

            if dataNLP.options.scaling
                solution.scaledVariables.resmin.X=scaled_solution.X;
                solution.scaledVariables.resmin.U=scaled_solution.U;
                solution.scaledVariables.X=scaled_solution.X(1:ns:end,:);
                solution.scaledVariables.U=scaled_solution.U(1:ns:end,:);
                if isfield(dataNLP.data,'Pscale')
                    solution.scaledVariables.p=scaled_solution.p;
                end
                solution.scaledVariables.x0=scaled_solution.X(1,:);
            end
            solution.coll.t0=solution.t0;
            solution.coll.tf=solution.tf;
            solution.t0=t0_minres;
            solution.tf=tf_minres;
            solution.p=p_minres;
            solution.resmin.X=minres_solution.X;
            solution.resmin.U=minres_solution.U;
            solution.resmin.T=minres_solution.T;
            solution.X=minres_solution.X(1:ns:end,:);
            solution.U=minres_solution.U(1:ns:end,:);
            solution.T=minres_solution.T(1:ns:end,:);
            solution.x0=minres_solution.X(1,:);
            t0=t0_minres;
            tf=tf_minres;
        elseif strcmp(dataNLP.options.resultRep,'default') || strcmp(dataNLP.options.resultRep,'res_min_final_default')
            if strcmp(dataNLP.options.discretization,'hermite')
                for i=1:n % Cubic Hermite interpolation
                    [Xp{i}, dXp{i}]=HSInterpolation(T,X(:,i),Fi(:,i));
                end
                for i=1:m % Piecewise quadratic polynomials
                    [Up{i}]=HSInterpolationU(T,U(:,i));
                end
            elseif strcmp(dataNLP.options.discretization,'trapezoidal')
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
            if strcmp(dataNLP.options.stateRep,'linear') && strcmp(dataNLP.options.discretization,'euler')
                for i=1:n % Linear interpolation
                    [Xp{i}, dXp{i}]=Linsplines(T,X(1:end,i),F(1:end,i));
                end
            elseif strcmp(dataNLP.options.stateRep,'quadratic') && (strcmp(dataNLP.options.discretization,'euler') || strcmp(dataNLP.options.discretization,'trapezoidal'))
                for i=1:n % Quadratic interpolation
                    [Xp{i}, dXp{i}]=Quadsplines(T,X(:,i),F(:,i));
                end
            elseif strcmp(dataNLP.options.stateRep,'cubic') && strcmp(dataNLP.options.discretization,'hermite')
                for i=1:n % Cubic Hermite interpolation
    %                 [Xp{i}, dXp{i}]=Hsplines(T,X(:,i),F(:,i));
                    [Xp{i}, dXp{i}]=HSInterpolation(T,X(:,i),Fi(:,i));
                end
            elseif strcmp(dataNLP.options.stateRep,'pchip')
                for i=1:n
                    Xp{i}=pchip(T,X(:,i)');
                    dXp{i}=pchip(T,F(:,i)');
                end
            else
                error('State representation method invalid or not supported by the selected discretization method');
            end

            % Input trajectory
            if strcmp(dataNLP.options.inputRep,'constant')
                for i=1:m % Piecewise constant polynomials
                    Up{i}=mkpp(T,U(1:end-1,i)');
                end
            elseif strcmp(dataNLP.options.inputRep,'linear')
                for i=1:m % Piecewise linear polynomials
                    Up{i}=Linsplines(T,U(:,i),diff(U(:,i)));
                end
            elseif strcmp(dataNLP.options.inputRep,'quadratic') && strcmp(dataNLP.options.discretization,'hermite')
                for i=1:m % Piecewise quadratic polynomials
                    [Up{i}]=HSInterpolationU(T,U(:,i));
                end
            elseif strcmp(dataNLP.options.inputRep,'pchip')
                for i=1:m 
                    Up{i}=pchip(T,U(:,i)');
                end
            else
                error('Input representation method invalid or not supported by the selected discretization method');
            end
        end
    end

    % Calculate variable rates
    if ~strcmp(dataNLP.options.transcription,'multiple_shooting')
        [ solution ] = getVariableRate( solution,dataNLP );
    end
    
    solution.Xp=Xp;
    solution.dXp=dXp;
    solution.Up=Up;
    
    % Display discretization error, constaint violation and number of active
    % constraints
    [Error, ErrorRelative,T_error]=estimateError(solution,p,t0,tf,n,m,f,M,ns,dataNLP);
    
    [ConstraintError,T_ConstraintError,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation(p,t0,tf,n,m,solution,problem,dataNLP,1);
    


    % Obtain multipliers
    if ~strcmp(dataNLP.options.transcription,'integral_res_min') && ~strcmp(dataNLP.options.NLPsolver,'NOMAD')
        switch(dataNLP.options.NLPsolver)
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
        solution.multipliers.lambda_1toN=lambda;
    else
        solution.multipliers.lambda_1toN=[];
    end
    


    % Save results
    solution.Error=Error;
    solution.T_error=T_error;
    solution.ErrorRelative=ErrorRelative;

    solution.ConstraintError=ConstraintError;
    solution.T_ConstraintError=T_ConstraintError;

    solution.ActiveConstraint=ActiveConstraint;
    solution.NumActiveConstraint=NumActiveConstraint;
         
    solution.MaxAbsError=max(solution.Error);
    solution.MaxRelError=max(solution.ErrorRelative);
    solution.MaxConstVioError=max(solution.ConstraintError);

    if isfield(dataNLP.options.print,'residual_error') && dataNLP.options.print.residual_error
        if strcmp(dataNLP.options.discretization,'hermite')
            [r,r_seg]=estimateResidual(solution,p,t0,tf,n,m,data);
        else
            r_seg=Error';
            r=sum(r_seg,2);
        end
        solution.residuals.r=r;
        solution.residuals.r_seg=r_seg;
    end
end
%------------- END OF CODE --------------

end


function [Error,ErrorRelative,T]=estimateError(solution,p,t0,tf,n,m,f,M,ns,data)
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
T=(tf-t0)*[0;cumsum(data.tau)]*data.Nm/ns+t0;
vdat=data.data;

%Computation of scale weights for the discretization error
wi=zeros(n,1);
for i=1:n
     wi(i)=max(max(abs(speval(solution,'X',i,T)),abs(speval(solution,'dX',i,T))));
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
    xp=speval(solution,'X',i,tq);
    Xq(:,i)=xp(:);
    dxp=speval(solution,'dX',i,tq);
    dXq(:,i)=dxp(:); 
end

Uq=zeros(length(tq),m);
for i=1:m
    up=speval(solution,'U',i,tq); 
    Uq(:,i)=up(:);
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p',k1(1),1);
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


function [MaxConstraintError,tg,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation(p,t0,tf,n,m,solution,problem,data,id)

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

ntg=length(solution.T)*10;
tg=linspace(t0,tf,ntg);
Xg=zeros(ntg,n);

for i=1:n
    Xg(:,i)=speval(solution,'X',i,tg);
end

Ug=zeros(ntg,m);
for i=1:m
    Ug(:,i)=speval(solution,'U',i,tg); 
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p',ntg,1);
else
    P=[];
end


g=data.functions_unscaled{4};

G = g(Xg,Ug,P,tg',data.data);           % Function evaluations.


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

function [r,r_seg]=estimateResidual(solution,p,t0,tf,n,m,data)

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


if isfield(data,'resmin')
    ntg=length(data.resmin.tau_quad);
    tg=data.resmin.tau_quad*(tf-t0)+t0;
    tcoll=data.resmin.tau*(tf-t0)+t0;
    vdat=data.data;
    sum_nps_quad=data.resmin.sum_nps_quad;
    idx_quad=data.resmin.idx_quad;
    LGRweights=data.resmin.LGR.weights;
    nps=data.resmin.nps;
    tau=data.resmin.tau;
    DxHS_hf=data.resmin.DxHS_hf;
    DxHS_p1=data.resmin.DxHS_p1;
    AfHS=data.resmin.AfHS;
    f=data.functions_unscaled{3};
elseif isfield(data,'dataNLP')
    ntg=length(data.tau_quad);
    tg=data.tau_quad*(tf-t0)+t0;
    tcoll=data.tau*(tf-t0)+t0;
    vdat=data.dataNLP.data;
    sum_nps_quad=data.sum_nps_quad;
    idx_quad=data.idx_quad;
    LGRweights=data.LGR.weights;
    nps=data.nps;
    tau=data.tau;
    DxHS_hf=data.DxHS_hf;
    DxHS_p1=data.DxHS_p1;
    AfHS=data.AfHS;
    f=data.dataNLP.functions_unscaled{3};
end

Xg=zeros(ntg,n);
Xcoll=zeros(length(tcoll),n);

for i=1:n
    Xg(:,i)=speval(solution,'X',i,tg);
    Xcoll(:,i)=speval(solution,'X',i,tcoll);
end

Ug=zeros(ntg,m);
Ucoll=zeros(length(tcoll),m);
for i=1:m
    Ug(:,i)=speval(solution,'U',i,tg); 
    Ucoll(:,i)=speval(solution,'U',i,tcoll); 
end


if ~isfield(data,'resmin')
    F_k=f(Xcoll(1:2:end,:),Ucoll(1:2:end,:),repmat(p',length(tau),1),tau(1:2:end)*(tf-t0)+t0,vdat);
    F_kph=DxHS_hf*Xcoll/(tf-t0)-F_k(1:end-1,:)/2;
    F_kp1=DxHS_p1*Xcoll/(tf-t0)+F_k(1:end-1,:);
    F=[F_k(1:end-1,:) F_kph F_kp1]';
    F=reshape(F(:),n,3*nps)';
    dXg=AfHS*F;
else
    dXg=zeros(ntg,n);
    for i=1:n
        dXg(:,i)=speval(solution,'dX',i,tg);
    end
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p',ntg,1);
else
    P=[];
end




F=f(Xg,Ug,P,tg,vdat);           % Function evaluations.


Res=dXg-F;
r_seg=1/((tf-t0).^2).*transpose(sum_nps_quad*(repmat([LGRweights{1, idx_quad(1)};0],nps,1).*(Res.^2)));
r=sum(r_seg,2);
end
%------------- END OF CODE --------------

function [Error, ErrorRelative]=estimateError_LGR(solution,p,T,n,m,f,M,data)
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
TSeg_Bar=solution.TSeg_Bar;
tf=TSeg_Bar(end);
T=[T;tf];

vdat=data.data;
%Computation of scale weights for the discretization error

wi=zeros(n,1);
% tau = normalizeT( [T;tf], t0,tf );
for i=1:n
     wi(i)=max(max(abs(speval( solution,'X',i,T)),abs(speval( solution,'X',i,T))));
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
        xp=speval( solution,'X',i,tq);
        Xq(:,i)=xp(:);
        dxp=speval( solution,'dX',i,tq);
        dXq(:,i)=dxp(:); 
    end

    Uq=zeros(length(tq),m);
    for i=1:m
        up=speval( solution,'U',i,tq);
        Uq(:,i)=up(:);
    end

    % Extract design parameters if specified and convert to cells
    if ~isempty(p) 
        P=repmat(p',k1(1),1);
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

function [MaxConstraintError,tg,ActiveConstraint,NumActiveConstraint]=estimateConstraintViolation_LGR(p,t0,tf,n,m,solution,problem,data,id)
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
ntg=length(solution.T)*10;
tg=linspace(t0,tf,ntg);

Xg=zeros(ntg,n);

for i=1:n
    Xg(:,i)=speval(solution,'X',i,tg);
end

Ug=zeros(ntg,m);
for i=1:m
    Ug(:,i)=speval(solution,'U',i,tg);
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p',ntg,1);
else
    P=[];
end

g=data.functions_unscaled{4};
G = g(Xg,Ug,P,tg',data.data);           % Function evaluations.

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

function [r,r_seg]=estimateResidual_LGR(solution,p,t0,tf,n,m,data)
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

if isfield(data,'resmin')
    ntg=length(data.resmin.tau_quad);
    tg=data.resmin.tau_quad*(tf-t0)+t0;
    vdat=data.data;
    sum_nps_quad=data.resmin.sum_nps_quad;
    idx_quad=data.resmin.idx_quad;
    LGRweights=data.resmin.LGR.weights;
    nps=data.resmin.nps;
    f=data.functions_unscaled{3};
else
    ntg=length(data.tau_quad);
    tg=data.tau_quad*(tf-t0)+t0;
    vdat=data.dataNLP.data;
    sum_nps_quad=data.sum_nps_quad;
    idx_quad=data.idx_quad;
    LGRweights=data.LGR.weights;
    nps=data.nps;
    f=data.dataNLP.functions_unscaled{3};
end

Xg=zeros(ntg,n);
dXg=zeros(ntg,n);
for i=1:n
    Xg(:,i)=speval(solution,'X',i,tg);
end
for i=1:n
    dXg(:,i)=speval(solution,'dX',i,tg);
end

Ug=zeros(ntg,m);
for i=1:m
    Ug(:,i)=speval(solution,'U',i,tg);
end

% Extract design parameters if specified and convert to cells
if ~isempty(p) 
    P=repmat(p',ntg,1);
else
    P=[];
end


F=f(Xg,Ug,P,tg,vdat);           % Function evaluations.

Res=dXg-F;
r_seg=1/((tf-t0).^2).*transpose(sum_nps_quad*(repmat([LGRweights{1, idx_quad(1)};0],nps,1).*(Res.^2)));
r=sum(r_seg,2);

end
