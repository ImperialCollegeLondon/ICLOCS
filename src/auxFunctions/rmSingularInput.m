function [newsolp] = rmSingularInput(Up)
%Remove singular arc ringing from linear problems based on aera rule

coefs=Up.coefs;
breaks=Up.breaks';
areafun=@(coef,val)coefs(:,1).*val.^3+coefs(:,2).*val.^2+coefs(:,3).*val;
area=areafun(coefs,breaks(2:end,1)-breaks(1:end-1,1));
newsol=area./diff(breaks);
newsolp=mkpp(breaks,newsol);


