function [ ] = genAdigator4ICLOCS_resmin( options, data, n, m, np, nt, M, ng_eq )
%genAdigator4ICLOCS_resmin - Initialize Adigator for use with integrated residual minimization method
%
% Syntax:  [ ] = genAdigator4ICLOCS_resmin( options, data, n, m, np, nt, M )
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

switch options.discretization
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        xidx=[data.dataNLP.FD.index.Ly(:,1:n);data.dataNLP.FD.index.Ly(end,1:n)+1];
        uidx=data.dataNLP.FD.index.Ly(:,(n+1):(n+m));

        [x_idx2,odr]=sort(xidx(:));
        x_idx1=(1:(M+1)*n).';
        x_idx1=x_idx1(odr);

        [u_idx2,odr]=sort(uidx(:));
        u_idx1=(1:(M)*m).';
        u_idx1=u_idx1(odr);
        
        nz=nt+np+(M+1)*n+(M)*m;
        
        gX = adigatorCreateDerivInput([M+1,n],struct('vodname','Y','vodsize',[nz,1],...
          'nzlocs',[x_idx1 x_idx2]));
        gU = adigatorCreateDerivInput([M,m],struct('vodname','Y','vodsize',[nz,1],...
          'nzlocs',[u_idx1 u_idx2]));
        if nt
            gT = adigatorCreateDerivInput([1,nt],struct('vodname','Y','vodsize',[nz,1],...
              'nzlocs',[(1:nt).' (((M+1)*n+(M)*m+np+1):((M+1)*n+(M)*m+np+nt)).']));
        else
            gT=[data.dataNLP.t0 data.dataNLP.tf];
        end
        if np
           gP = adigatorCreateDerivInput([1,np],struct('vodname','Y','vodsize',[nz,1],...
          'nzlocs',[(1:np).' (((M+1)*n+(M)*m+1):((M+1)*n+(M)*m+np)).']));
        else
           gP=sparse(1,np);
        end

        % Create 1st Derivative File
        if options.scaling
            if ng_eq
                adigator('costResidualMin_ModeMinRes_Adigator_Scaling_DAE',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            else
                adigator('costResidualMin_ModeMinRes_Adigator_Scaling',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            end
        else
            if ng_eq
                adigator('costResidualMin_ModeMinRes_Adigator_DAE',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            else
                adigator('costResidualMin_ModeMinRes_Adigator',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            end
        end


        % Create 2nd Derivative File
        if ~strcmp(options.ipopt.hessian_approximation,'limited-memory')
            gXs.f = gX; gUs.f = gU; gTs.f=gT; gPs.f=gP;
            gXs.dY = adigatorCreateAuxInput([(M+1)*n,1]);
            gUs.dY = adigatorCreateAuxInput([(M)*m,1]);
            gTs.dY = adigatorCreateAuxInput([nt,1]);
            if np
                gPs.dY = adigatorCreateAuxInput([np,1]);
            else
                gPs=sparse(1,np);
            end
            adigator('minresCost_Y',{gXs,gUs,gPs,gTs,data},'minresCost_YY',...
              adigatorOptions('overwrite',1,'comments',0));
        end

    otherwise % h Transcription Method
        xidx=data.dataNLP.FD.index.Ly(:,(nt+np+1):(nt+np+n));
        uidx=data.dataNLP.FD.index.Ly(:,(nt+np+n+1):(nt+np+n+m));

        [x_idx2,odr]=sort(xidx(:));
        x_idx1=(1:M*n).';
        x_idx1=x_idx1(odr);

        [u_idx2,odr]=sort(uidx(:));
        u_idx1=(1:M*m).';
        u_idx1=u_idx1(odr);
        gX = adigatorCreateDerivInput([M,n],struct('vodname','Y','vodsize',[M*(m+n)+np+nt,1],...
          'nzlocs',[x_idx1 x_idx2]));
        gU = adigatorCreateDerivInput([M,m],struct('vodname','Y','vodsize',[M*(m+n)+np+nt,1],...
          'nzlocs',[u_idx1 u_idx2]));
        if nt>=2
            gT = adigatorCreateDerivInput([1,nt],struct('vodname','Y','vodsize',[M*(m+n)+np+nt,1],...
              'nzlocs',[(1:nt).' (1:nt).']));
        else
            gT=[data.dataNLP.t0;data.dataNLP.tf];
        end
        if np
           gP = adigatorCreateDerivInput([1,np],struct('vodname','Y','vodsize',[M*(m+n)+np+nt,1],...
          'nzlocs',[(1:np).' (nt+1:nt+np).']));
        else
           gP=sparse(1,np);
        end

        % Create 1st Derivative File
        if options.scaling
            if ng_eq
                adigator('costResidualMin_ModeMinRes_Adigator_Scaling_DAE',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            else
                adigator('costResidualMin_ModeMinRes_Adigator_Scaling',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            end
        else
            if ng_eq
                adigator('costResidualMin_ModeMinRes_Adigator_DAE',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            else
                adigator('costResidualMin_ModeMinRes_Adigator',{gX,gU,gP,gT,data},'minresCost_Y',adigatorOptions('overwrite',1));
            end
        end


        % Create 2nd Derivative File
        if ~strcmp(options.ipopt.hessian_approximation,'limited-memory')
            gXs.f = gX; gUs.f = gU; gTs.f=gT; gPs.f=gP;
            gXs.dY = adigatorCreateAuxInput([M*n,1]);
            gUs.dY = adigatorCreateAuxInput([M*m,1]);
            gTs.dY = adigatorCreateAuxInput([nt,1]);
            if np
                gPs.dY = adigatorCreateAuxInput([np,1]);
            else
                gPs=sparse(1,np);
            end
            adigator('minresCost_Y',{gXs,gUs,gPs,gTs,data},'minresCost_YY',...
              adigatorOptions('overwrite',1,'comments',0));
        end
end



end

