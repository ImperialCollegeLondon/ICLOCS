function [ scaledvari ] = scale_variables_back( vari, vscales, vshifts )
%UNTITLED

scaledvari= bsxfun(@times, bsxfun(@minus, vari, vshifts),1./vscales);

end

