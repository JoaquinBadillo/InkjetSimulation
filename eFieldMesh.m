function [Ex,Ey] = eFieldMesh(rcx,rcy,rpx,rpy,q)
% Computes the electric field in multiple points due to a system of charged particles

% - rcx: horizontal coordinates of the charged particles
% - rcy: vertical coordinates of the charged particles
% - rpx: horizontal coordinates of the points
% - rpy: vertical coordinates of the points
% - q: charges

eps0=8.854e-12; % Epsiloon naught
k=1/(4*pi*eps0); % Coulomb's constant

Ex=0; Ey=0; % Initialize variables

for i=1:length(rcx)
    dx=rpx-rcx(i); % Array of horizontal distances
    dy=rpy-rcy(i); % Array of vertical distances
    r=sqrt(dx.^2+dy.^2);% Array of L2 Norms
    Ex=Ex+k*q(i).*dx./r.^3; % Horizontal components of the field
    Ey=Ey+k*q(i).*dy./r.^3; % Vertical components of the field
end
return
end