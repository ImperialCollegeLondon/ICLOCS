
function H=overlapping(M)
%   OVERLAPPING - It  groups orthogonal colums of the matrix M together. 
%   Each group is associated to an integer number starting from 1,2 ... 
%
%   Syntax:  H=overlapping(M)
%
% Inputs: 
%    M- matrix  
%
% Outputs:
%    H - vector containing the set number for the column of M 
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


 [r,c]=size(M);
 H=zeros(c,1);
if ~isempty(M)
  H(1,1)=1;
  numsets=1;  % number of sets 
  elem_set=zeros(1,c);
  elem_set(1)=1; % number of elemets in the set
  sets{1}=M(:,1);
 
  for k=2:c  %it checks to which set the column k belongs to 
    %it checks if the column k belongs to the set i  
    for i=1:numsets  
     flag=1;  
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % it performs the check with any element belonging to the set i
 
    for j=1:elem_set(i) 
       test_ov=M(:,k).*sets{i}(:,j);
       set_belong=find(test_ov,1);
       if ~isempty(set_belong)
         flag=0;
         break
       end
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%  

   if flag==1
    sets{i}=[sets{i},M(:,k)];
    elem_set(i)= elem_set(i)+1;  
    H(k,1)=i;
    break
   else
     if (i==numsets)
       numsets=numsets+1;  
       sets{numsets}=M(:,k); 
       H(k,1)=numsets;
       elem_set(numsets)=1;
     end    
   end
 
    end
 end    

end