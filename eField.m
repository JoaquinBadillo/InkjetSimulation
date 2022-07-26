function [Ex,Ey] = eField(rcx,rcy,rpx,rpy,q)
% Computes the electric field at a given point.

% - rcx: horizontal coordinates of the charged particles
% - rcy: vertical coordinates of the charged particles
% - rpx: horizontal coordinate of the point
% - rpy: vertical coordinate of the point
% - q: charges

eps0=8.854e-12; % epsilon naught
k=1/(4*pi*eps0); % Coulomb's Constant

dx=rpx-rcx; % Horizontal distance
dy=rpy-rcy; % Vertical distance
r=sqrt(dx.^2+dy.^2); % L2 Norm

% Electric Field
Ex=k*q.*dx./r.^3;
Ey=k*q.*dy./r.^3;
Ex=sum(Ex);
Ey=sum(Ey);

return
end