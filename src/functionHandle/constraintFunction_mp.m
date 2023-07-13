function constraints_mp=constraintFunction_mp(z_mp,auxdata)
%constraintFunction_mp - Evaluate the constraint functions for a multiphase
%problem
%
% Syntax:  constraints_mp=constraintFunction_mp(z_mp,auxdata)
%
% Inputs:
%    z    - NLP variable 
%    data - Structure with data needed to evaluate ConstFn
%
% Outputs:
%    constraints - Vector of values of the constraint function at x
%
% Other m-files required: multipleShooting.m or directTranscription.m
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
nbl=auxdata.mpdata.mpsizes.nbl;
% xu0f_mp=zeros(1,auxdata.mpdata.mpsizes.nxu0f_inc(end));

constraints_mp=zeros(auxdata.mpdata.mpsizes.nConst-nbl,1);

x0=cell(auxdata.mpdata.mpsizes.nphase,1);xf=cell(auxdata.mpdata.mpsizes.nphase,1);
u0=cell(auxdata.mpdata.mpsizes.nphase,1);uf=cell(auxdata.mpdata.mpsizes.nphase,1);

for i=1:length(auxdata.phasedata)
    data=auxdata.phasedata{i};
    z=zeros(data.nz,1);
    z(data.zidx.org.z,1)=z_mp(data.zidx.mp.z,1);
    
    switch data.options.transcription

        case {'multiple_shooting'} 

            constraints=multipleShooting('const',z,data);

        % Direct Collocation
        case {'direct_collocation','direct_collocation_intres_reg'}

            switch data.options.discretization

                case{'discrete','euler','trapezoidal','hermite'}

                    [constraints,XU0f]=directCollocation('const',z,data,i);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;
                     
                case{'globalLGR','hpLGR'}

                    [constraints,XU0f]=directCollocationLGR('const',z,data,i);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;

                case {'resMinConstructionForSolution'}

                    constraints=resMinConstructionForSolution('const',z,data);

                case {'resMinInterpolationForSolution'}

                    [constraints,XU0f]=resMinOptimization('const',z,data);
                     x0{i}=XU0f.x0;
                     xf{i}=XU0f.xf;
                     u0{i}=XU0f.u0;
                     uf{i}=XU0f.uf;
            end

        case {'integral_res_min'}

            [constraints,XU0f]=resMinOptimization('const',z,data);
             x0{i}=XU0f.x0;
             xf{i}=XU0f.xf;
             u0{i}=XU0f.u0;
             uf{i}=XU0f.uf;
    end
    
    constraints_mp(data.zidx.mp.const,1)=constraints_mp(data.zidx.mp.const,1)+constraints(data.zidx.org.const,1);

    
end

linkfunctions=auxdata.mpdata.linkfunctions;
t0=z_mp(end-auxdata.mpdata.mpsizes.nt+1:end-1);
tf=z_mp(end-auxdata.mpdata.mpsizes.nt+2:end);
p=z_mp(end-auxdata.mpdata.mpsizes.nt-auxdata.mpdata.mpsizes.np+1:end-auxdata.mpdata.mpsizes.nt);
if isfield(data.data,'Pscale')
   p=scale_variables_back( p', data.data.Pscale, data.data.Pshift )';
end

[blink_l,blink_nl]=linkfunctions(x0,xf,u0,uf,p,t0,tf,auxdata.mpdata.data);

constraints_mp=[constraints_mp;blink_l;blink_nl];

end

%------------- END OF CODE --------------