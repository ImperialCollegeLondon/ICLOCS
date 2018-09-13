function [ T_active ] = identifyConstActive( problem, options, solution, lambda_g,ng,ig )
%identifyConstActive - identify time intervals that the path constraints are
%potentially active
%
% Syntax:   [ T_active ] = identifyConstActive( options, solution, lambda_g )
%
% Inputs:
%    options, solution - Defined in transcribeOCP
%    lambda_g - the corresponding Lagrange multiplier values
%
% Outputs:
%    T_active - time intervals the path constraints are potentially active
%     
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

multi_g=mat2gray(abs(lambda_g(2:end-1)));

if (strcmp(options.transcription,'globalLGR')) || (strcmp(options.transcription,'hpLGR'))
    time_vec=[solution.T;solution.tf]';
elseif strcmp(options.transcription,'hermite')
    time_vec=solution.org.T';
else
    time_vec=solution.T';
end

[ipt,~] = findchangepts(multi_g(2:end-1),'Statistic','linear','MinThreshold',0.1);

ipt=unique([1;ipt+1;length(multi_g(:,1))+2]);
ipt=repelem(ipt,2);ipt(1)=[];ipt(end)=[];
ipt=reshape(ipt,2,length(ipt)/2)';

temp=1:length(multi_g(:,1));

temp1=cell(1,size(ipt,1));
avg_g=zeros(1,size(ipt,1));
for i=1:size(ipt,1)
    temp1{1,i}=temp(temp>=ipt(i,1) & temp<=ipt(i,2));
    avg_g(1,i)=mean(multi_g(temp1{1,i},1));
end
active_idx=avg_g>0.05;
ipt_active=ipt(active_idx,:);
T_active=time_vec(ipt_active);


time_vec=solution.T_ConstraintError;
temp=time_vec(solution.ActiveConstraint.Glpactive(:,ig));
if ~isempty(temp)
    temp=repelem(temp,2);temp(1)=[];temp(end)=[];
    temp=reshape(temp,2,length(temp)/2)';
    T_active=[T_active;temp];
end

temp=time_vec(solution.ActiveConstraint.Gupactive(:,ig));
if ~isempty(temp)
    temp=repelem(temp,2);temp(1)=[];temp(end)=[];
    temp=reshape(temp,2,length(temp)/2)';
    T_active=[T_active;temp];
end


% for j=1:2
%     time_vec=solution.T_ConstraintError;
%     [ipt,~] = findchangepts(solution.ConstraintError(:,ig+(j-1)*ng),'Statistic','linear','MinThreshold',0.1);
% 
%     if ~isempty(ipt)
%         ipt=repelem(ipt,2);ipt(1)=[];ipt(end)=[];
%         ipt=reshape(ipt,2,length(ipt)/2)';
% 
%         temp=1:length(solution.T_ConstraintError);
%         temp1=cell(1,size(ipt,1));
%         max_g=zeros(1,size(ipt,1));
%         for i=1:size(ipt,1)
%             temp1{1,i}=temp(temp>=ipt(i,1) & temp<=ipt(i,2));
%             max_g(1,i)=max(abs(solution.ConstraintError(temp1{1,i},ig+(j-1)*ng)));
%         end
%         active_idx=max_g>problem.constraintErrorTol(ig+(j-1)*ng);
%         ipt_active=ipt(active_idx,:);
%         T_active=[T_active;time_vec(ipt_active)];
%     end
% end



T_active(:,1)=T_active(:,1)-0.1*solution.tf;
T_active(:,2)=T_active(:,2)+0.1*solution.tf;
T_active(T_active<=0)=0;
T_active(T_active>=solution.tf)=solution.tf;
T_active=mergeTimeIntervals( T_active );

end

