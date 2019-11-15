function [ Ez ] = gradientCost_AN_E( dE,Ez,data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ez(data.costStruct.E==1)=[dE.dt0 dE.dtf dE.dp dE.dx0 dE.du0 dE.dxf dE.duf];
    
end

