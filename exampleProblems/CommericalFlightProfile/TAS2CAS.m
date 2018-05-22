function [ V_cas ] = TAS2CAS( kappa, p, rho, p0, rho0, V_tas )

V_cas=sqrt(p0/rho0*2*kappa/(kappa-1)*((p/p0*((V_tas^2*(kappa-1)/2/kappa*rho/p+1)^(kappa/(kappa-1))-1)+1)^((kappa-1)/kappa)-1));

end

