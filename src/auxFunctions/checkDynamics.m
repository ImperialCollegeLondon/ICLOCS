function checkDynamics( z0,data )
%checkDynamics - pre-testing of the user defined function, to catch errors more easily
%
% Syntax:  checkDynamics( z0,data )
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
% 
%------------- BEGIN CODE --------------


if isfield(data,'dataNLP')
    data=data.dataNLP;
end
if (strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))
    data.sizes{5}=data.sizes{18}+data.sizes{19};
else  
    data.sizes{5}=data.sizes{15}+data.sizes{16};
end
data.options.transcription='direct_collocation';
data.options.derivatives='numeric';
data.data.resmin=0;
constraintFunction(z0,data);
    

end

