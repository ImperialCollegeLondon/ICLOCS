function c=g_unscaled(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------
Dynamics=vdat.InternalDynamics;
ng_group=nargout(Dynamics);

if ng_group==1
    c=[];
elseif ng_group==2
    if vdat.resmin && vdat.ng_eq 
        c=[];
    else
        [~,c] = Dynamics(x,u,p,t,vdat);
    end
else
    [~,ceq,cneq] = Dynamics(x,u,p,t,vdat);
    if vdat.resmin && vdat.ng_eq 
        c=cneq;
    else
        c=[ceq cneq];
    end
end




%------------- END OF CODE ---------------------