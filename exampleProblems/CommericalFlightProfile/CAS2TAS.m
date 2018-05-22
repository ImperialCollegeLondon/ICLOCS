function [ V_tas ] = CAS2TAS( kappa, p, rho, p0, rho0, V_cas )

V_tas=sqrt(2*kappa*p/(kappa-1)./rho.*((1+p0./p.*((1+(kappa-1)/2/kappa*rho0/p0.*V_cas.^2).^(kappa/(kappa-1))-1)).^((kappa-1)/kappa)-1));

end

