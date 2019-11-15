function [ Ez ] = gradientCost_LGR_AN_E( dE,Ez,data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ez(data.costStruct.E==1)=[dE.dx0 dE.dxf dE.du0 dE.duf dE.dp dE.dt0 dE.dtf ];
    
end

