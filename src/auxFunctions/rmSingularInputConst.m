function [newsolp] = rmSingularInputConst(Up)
%Remove singular arc ringing from linear problems based on aera rule

coefs=Up.coefs;
breaks=Up.breaks';
areafun=@(coef,val)coefs(:,1).*val.^3/3+coefs(:,2).*val.^2/2+coefs(:,3).*val;
DT=breaks(2:end,1)-breaks(1:end-1,1);
area=areafun(coefs,DT);
newsol=area./DT;
newsolp=mkpp(breaks,newsol);


