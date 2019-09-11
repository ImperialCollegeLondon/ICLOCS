function [LGR]=generateLGR_All(npd)
%generateLGR_All - Generate the points, weights and differentiation
%metrices for all segments with various degrees (LGR direct transcription
%formulation)
%
% Syntax:  LGR=directCollocation(npd)
%
% Inputs:
%    npd - A vector contains different polynomial degrees used in the mesh
%
% Outputs:
%    LGR - Data structure containing the LGR points, weights and differentiation
%          metrices
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
for inpd=1:length(npd)
    
    % Solve root finding problem with eigenvalue computation of the Jacobi matrix
    N=npd(inpd);
    indices=0:N-2;                  
    maindiagonal=diag(1./((2*indices+1).*(2*indices+3))); %Diagonal element
    indices=1:N-2;
    Jacobimat=maindiagonal+diag(sqrt(indices.*(indices+1))./(2*indices+1),1) +diag(sqrt(indices.*(indices+1))./(2*indices+1),-1); %Jacobi matrix
    LGR_points= [-1; sort(eig(sparse(Jacobimat)))]; %LGR points
    Legendrepoly=legendreP(N-1,LGR_points);
    LGR_weights=(1-LGR_points)./(N^2*Legendrepoly.^2); %LGR weights
    LGR_points_include_noncol=[LGR_points; 1]; 

    % Symbolic Check
    % syms PRs
    % expr1 = legendreP(N-1,PRs);
    % expr2 = legendreP(N,PRs);
    % rts=vpasolve(expr1+expr2 == 0);
    % difexp1=diff(expr1,PRs,1);
    % PRs=LGR_Points;
    % Leg_Poly_Dot=eval(difexp1);
    % LGR_Weights = [2/N^2; (1./(1-LGR_Points(2:end)))./Leg_Poly_Dot(2:end).^2];

    % Calculate the LGR differentiation matrix
    N=N+1;
    tau=LGR_points_include_noncol;
    Diff_Matrix_LGR = zeros(N,N);
    for id=1:N
        z=tau(id);
        for j=1:N
            temp=0;
            for i=1:N
                if not(i==j)
                    k = 1/(tau(j)-tau(i));
                    for m=1:N
                        if not(m==j) && not(m==i)
                            k = k*(z-tau(m))/(tau(j)-tau(m));
                        end
                    end
                    temp = temp + k;
                end
            end
        Diff_Matrix_LGR(id,j)=temp;
        end
    end
    Diff_Matrix_LGR_square=Diff_Matrix_LGR;
    Diff_Matrix_LGR = Diff_Matrix_LGR(1:end-1,:);
    Diag_Diff_Matrix_LGR=zeros(size(Diff_Matrix_LGR));
    Diag_Diff_Matrix_LGR(:,1:end-1) = diag(diag(Diff_Matrix_LGR));

    %Save variables in a structure
    LGR.points{inpd}=LGR_points;
    LGR.weights{inpd}=LGR_weights;
    LGR.diff_matrix{inpd}=Diff_Matrix_LGR;
    LGR.diff_matrix_square{inpd}=Diff_Matrix_LGR_square;
    LGR.diag_diff_matrix{inpd}=Diag_Diff_Matrix_LGR;
    
end
