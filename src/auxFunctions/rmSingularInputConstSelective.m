function [newTraj] = rmSingularInputConstSelective(Up,idxIntSing)
%Remove singular arc ringing from linear problems based on aera rule

if size(Up.coefs,2)==3
    areafun=@(coef,val)coef(:,1).*val.^3/3+coef(:,2).*val.^2/2+coef(:,3).*val;
else
    error('the method currently not supported')
end

newTraj=Up;
newcoef=Up.coefs;
% newcoef=cell(size(idxIntSing));
for i=1:length(idxIntSing)
    coefs=Up.coefs(idxIntSing{i},:);
    idx_bks=union(idxIntSing{i},idxIntSing{i}+1);
    breaks=Up.breaks(idx_bks)';
    DT=breaks(2:end,1)-breaks(1:end-1,1);
    area=areafun(coefs,DT);
    if size(Up.coefs,2)==3
        newcoef(idxIntSing{i},:)=[zeros(length(area),2) area./DT];
    end
end

newTraj.coefs=newcoef;
