function [ scaledvari ] = scale_variables( vari, vscales, vshifts )
%UNTITLED

scaledvari= bsxfun(@plus, bsxfun(@times, vari, vscales),vshifts);

end

