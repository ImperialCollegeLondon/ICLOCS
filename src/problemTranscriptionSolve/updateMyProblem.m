function [OCP_updated] = updateMyProblem(OCP_old,varargin)
% updateMyProblem - Update problem formulation without redo transcription
%
% Syntax:  
%          [OCP_updated] = updateMyProblem(OCP_old, 'options', values)	
% 
% Inputs:
%    OCP_old - Previous OCP formulation
%    options - field to be updated ('x0', 'z0', 'bu', 'bl', 'userdata')
%    values - values corresponding to options
% 
% Output:
%    OCP_updated - Updated OCP formulation
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk



   p = inputParser;
   validScalarPosNum = @(x) isnumeric(x);
   addOptional(p,'x0','none',validScalarPosNum);
   addOptional(p,'z0','none',validScalarPosNum);
   addOptional(p,'bu','none',validScalarPosNum);
   addOptional(p,'bl','none',validScalarPosNum);
   addOptional(p,'userdata','none',@isstruct);
   
   parse(p,varargin{:});
   
   OCP_updated=OCP_old;
   
   if ~strcmp(OCP_updated.data.options.transcription,'direct_collocation') || isfield(OCP_updated.data.options,'mp')
       error('Problem update currectly only supports single-phase direct collocation formulation')
   end
   
   if (strcmp(OCP_updated.data.options.discretization,'globalLGR')) || (strcmp(OCP_updated.data.options.discretization,'hpLGR'))
        [nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive,ng_eq,ng_neq]=deal(OCP_updated.data.sizes{1:19});
   else
        [nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive,nps,ng_eq,ng_neq]=deal(OCP_updated.data.sizes{1:16});
   end

   if isnumeric(p.Results.z0)
      if size(OCP_updated.infoNLP.z0,1)~=size(p.Results.z0,1)
          warning('The dimension of z0 changes for the update');
      end
      OCP_updated.infoNLP.z0=p.Results.z0;   % Update the initial condition from new system information
   end

   if isnumeric(p.Results.x0)
       
       if OCP_updated.data.options.scaling
          x0=scale_variables(p.Results.x0,OCP_updated.data.data.Xscale,OCP_updated.data.data.Xshift);  
       else
          x0=p.Results.x0';
       end     
     
       if (strcmp(OCP_updated.data.options.discretization,'globalLGR')) || (strcmp(OCP_updated.data.options.discretization,'hpLGR')) 
           
           OCP_updated.infoNLP.z0(1:M+1:n*(M+1),1)=x0;   
           OCP_updated.infoNLP.zl(1:M+1:n*(M+1),1)=x0;
           OCP_updated.infoNLP.zu(1:M+1:n*(M+1),1)=x0;
           OCP_updated.data.x0t(:,1)=x0;
           OCP_updated.data.x0=x0';
       else
           OCP_updated.infoNLP.z0(nt+np+1:nt+np+n,1)=x0;   
           OCP_updated.infoNLP.zl(nt+np+1:nt+np+n,1)=x0;
           OCP_updated.infoNLP.zu(nt+np+1:nt+np+n,1)=x0;
           OCP_updated.data.x0t(:,1)=x0;
           OCP_updated.data.x0=x0';
       end
   end
   

   
   if isnumeric(p.Results.bu)
      OCP_updated.infoNLP.cu(end-nb+1:end)=p.Results.bu';   % Update the initial condition from new system information
      OCP_updated.data.map.bu=p.Results.bu;
   end
   
   if isnumeric(p.Results.bl)
      OCP_updated.infoNLP.cl(end-nb+1:end)=p.Results.bl';   % Update the initial condition from new system information
      OCP_updated.data.map.bl=p.Results.bl;
   end
   
   if isstruct(p.Results.userdata)
      OCP_updated.data.data=catstruct(OCP_updated.data.data,p.Results.userdata);
   end
   
   OCP_updated.data.infoNLP=OCP_updated.infoNLP;
   
end

