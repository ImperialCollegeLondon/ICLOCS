function [  Xout] = speval( solution,solType,nidx,T)
%speval - evaluate the spline functions
%
% Syntax:   Xout = speval( Xp,solType,nidx,T)
%
% Inputs:
%    solution  - solution to the OCP after post-processing.
%    solType - 'X' for states and 'U' for input variable
%    nidx - Return result for the nth state/input variable
%    T - Time vector
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


Xout=zeros(length(T),length(nidx));
for k=1:length(nidx)
    n=nidx(k);
    if ~isfield(solution,'TSeg_Bar')
        switch solType
            case{'x','X','state','states'}
                Pp=solution.Xp;
            case{'u','U','input','inputs'}
                Pp=solution.Up;
            case{'dx','dX','dynamics'}
                Pp=solution.dXp;
            case{'du','dU'}
                Pp=solution.dUp;
            otherwise
                error('requested solution type not defined!!');
        end
        Sp=ppval(Pp{n},T);
    else
        TSeg_Bar=solution.TSeg_Bar;
        switch solType
            case{'x','X','state','states'}
                Pp=solution.Xp;
                diffT_seg_Bar=ones(length(TSeg_Bar)-1,1);
            case{'u','U','input','inputs'}
                Pp=solution.Up;
                diffT_seg_Bar=ones(length(TSeg_Bar)-1,1);
            case{'dx','dX','dynamics'}
                Pp=solution.dXp;
                diffT_seg_Bar=diff((TSeg_Bar-TSeg_Bar(1))./(TSeg_Bar(end)-TSeg_Bar(1)))*(solution.tf-solution.t0);
            case{'du','dU'}
                Pp=solution.dUp;
                diffT_seg_Bar=diff((TSeg_Bar-TSeg_Bar(1))./(TSeg_Bar(end)-TSeg_Bar(1)))*(solution.tf-solution.t0);
            otherwise
                error('requested solution type not defined!!');
        end
        if isstruct(Pp{n})
            Sp=ppval(Pp{n},T);
        else
            Sp=zeros(length(T),1);Ppi=Pp(:,n);
            
            for i=1:length(TSeg_Bar)-1 
                Tau_Seg=normalizeT(T(T-TSeg_Bar(i)>=-1e-10 & T<TSeg_Bar(i+1)),TSeg_Bar(i),TSeg_Bar(i+1));
                if size(Ppi{i},1)==1
                    Sp((T-TSeg_Bar(i)>=-1e-10 & T<TSeg_Bar(i+1)),:)=polyval(Ppi{i},Tau_Seg/2+0.5)./diffT_seg_Bar(i);
                elseif size(Ppi{i},2)==1
                    Sp((T-TSeg_Bar(i)>=-1e-10 & T<TSeg_Bar(i+1)),:)=legendreEval(Ppi{i},Tau_Seg);
                else
                    Sp((T-TSeg_Bar(i)>=-1e-10 & T<TSeg_Bar(i+1)),:)=BaryEval(Ppi{i},Tau_Seg);
                end
            end
            if TSeg_Bar(end)==T(end)
                if size(Ppi{i},1)==1
                    Sp(end)=polyval(Ppi{i},1)./diffT_seg_Bar(i);
                elseif size(Ppi{i},2)==1
                    Sp(end)=legendreEval(Ppi{i},1);
                else
                    Sp(end)=BaryEval(Ppi{i},1);
                end
            end
        end
    end
    Xout(:,k)=Sp;
end

end

