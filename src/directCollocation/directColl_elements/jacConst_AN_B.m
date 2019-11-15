function [ bz ] = jacConst_AN_B( bz, db, nt, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if nt==1
        bz(data.costStruct.B==1)=[db.dtf db.dp db.dx0 db.du0 db.duf db.dxf];
    elseif nt==2
        bz(data.costStruct.B==1)=[db.dt0 db.dtf db.dp db.dx0 db.du0 db.duf db.dxf];
    else
        bz(data.costStruct.B==1)=[db.dp db.dx0 db.du0 db.duf db.dxf];
    end


end

