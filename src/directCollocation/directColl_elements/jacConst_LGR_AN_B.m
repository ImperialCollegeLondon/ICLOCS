function [ bz ] = jacConst_LGR_AN_B( bz, db, nt, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if data.options.adaptseg==1 
    bz(data.costStruct.B==1)=[db.dx0 db.dxf db.du0 db.duf db.dp db.dt];
else
    bz(data.costStruct.B==1)=[db.dx0 db.dxf db.du0 db.duf db.dp db.dt0 db.dtf];
end

end

