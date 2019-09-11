function [ data ] = genAdigator4ICLOCS( options, data, n, m, np, nt, M )
%genAdigator4ICLOCS - Initialize Adigator for use with direct collocation method 
%
% Syntax:  [ data ] = genAdigator4ICLOCS( options, data, n, m, np, nt, M )
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
        ng_eq=data.sizes{18};
        ng_neq=data.sizes{19};
        
        xidx=[data.FD.index.Ly(:,1:n);data.FD.index.Ly(end,1:n)+1];
        uidx=data.FD.index.Ly(:,(n+1):(n+m));

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
            gT=[data.t0 data.tf];
        end
        if np
           gP = adigatorCreateDerivInput([1,np],struct('vodname','Y','vodsize',[nz,1],...
          'nzlocs',[(1:np).' (((M+1)*n+(M)*m+1):((M+1)*n+(M)*m+np)).']));
        else
           gP=sparse(1,np);
        end

        % Create dynamics and path constraint 1st Derivative File
        FuncName=['Func_Y_' data.data.plantmodel];
        data.data.AdigatorUserFunction_Y=str2func(FuncName);
    
        if ng_eq && ng_neq
            if options.scaling
                adigator('userFunction_Adigator_LGR_All_Scaling',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            else
                adigator('userFunction_Adigator_LGR_All',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            end
        elseif ng_eq || ng_neq
            if options.scaling
                adigator('userFunction_Adigator_LGR_oneConst_Scaling',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            else
                adigator('userFunction_Adigator_LGR_oneConst',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            end
        else
            if options.scaling
                adigator('userFunction_Adigator_LGR_noConst_Scaling',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            else
                adigator('userFunction_Adigator_LGR_noConst',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            end
        end

        % Create dynamics and path constraint 2nd st Derivative File
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
            FuncName1=['Func_Y_' data.data.plantmodel];
            FuncName2=['Func_YY_' data.data.plantmodel];
            data.data.AdigatorUserFunction_YY=str2func(FuncName2);
            adigator(FuncName1,{gXs,gUs,gPs,gTs,data},FuncName2,...
              adigatorOptions('overwrite',1,'comments',0));
        end

    otherwise % h Transcription Method
        ng_eq=data.sizes{15};
        ng_neq=data.sizes{16};

        xidx=data.FD.index.Ly(:,(nt+np+1):(nt+np+n));
        uidx=data.FD.index.Ly(:,(nt+np+n+1):(nt+np+n+m));

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
        if nt
            gT = adigatorCreateDerivInput([1,nt],struct('vodname','Y','vodsize',[M*(m+n)+np+nt,1],...
              'nzlocs',[(1:nt).' (1:nt).']));
        else
            gT=data.tf;
        end
        if np
           gP = adigatorCreateDerivInput([1,np],struct('vodname','Y','vodsize',[M*(m+n)+np+nt,1],...
          'nzlocs',[(1:np).' (nt+1:nt+np).']));
        else
           gP=sparse(1,np);
        end

        % Create dynamics and path constraint 1st Derivative File
        FuncName=['Func_Y_' data.data.plantmodel];
        data.data.AdigatorUserFunction_Y=str2func(FuncName);
        if ng_eq && ng_neq
            if options.scaling
                adigator('userFunction_Adigator_All_Scaling',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            else
                adigator('userFunction_Adigator_All',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            end
        elseif ng_eq || ng_neq
            if options.scaling
                adigator('userFunction_Adigator_oneConst_Scaling',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            else
                adigator('userFunction_Adigator_oneConst',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            end
        else
            if options.scaling
                adigator('userFunction_Adigator_noConst_Scaling',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            else
                adigator('userFunction_Adigator_noConst',{gX,gU,gP,gT,data},FuncName,adigatorOptions('overwrite',1));
            end
        end

        % Create dynamics and path constraint 2nd st Derivative File
        if ~strcmp(options.ipopt.hessian_approximation,'limited-memory')
            gXs.f = gX; gUs.f = gU; gTs.f=gT; gPs.f=gP;
            % Need to make aux inputs for derivative fields since we dont know their
            % vectorized dimension
            gXs.dY = adigatorCreateAuxInput([M*n,1]);
            gUs.dY = adigatorCreateAuxInput([M*m,1]);
            gTs.dY = adigatorCreateAuxInput([nt,1]);
            if np
                gPs.dY = adigatorCreateAuxInput([np,1]);
            else
                gPs=sparse(1,np);
            end
            FuncName1=['Func_Y_' data.data.plantmodel];
            FuncName2=['Func_YY_' data.data.plantmodel];
            data.data.AdigatorUserFunction_YY=str2func(FuncName2);
            adigator(FuncName1,{gXs,gUs,gPs,gTs,data},FuncName2,...
              adigatorOptions('overwrite',1,'comments',0));
        end
end



end