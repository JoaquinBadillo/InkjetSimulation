clc; clear all; close all;

% Runge-Kutta step size
he = 7e-9; 

% Dimensions
L=0.5e-3;
s=0.2e-3;
d=0.5e-3;
h= 0.4e-3;

% Discrete model of the distribution
N=100; 
qtotal1=9e-9; 
qtotal2=-9e-9;
q1=ones(1,N)*qtotal1/N;
q2=ones(1,N)*qtotal2/N;
q=[q1,q2];

% Horizontal coordinates
rcx1=linspace(s,s+L,N); 
rcx=[rcx1,rcx1];

% Vertical coordinates
rcy1=ones(1,N)*(d/2);
rcy2=ones(1,N)*(-d/2);
rcy=[rcy1,rcy2];

% Parameters
m = 1.13e-10; % Mass
v0 = 20; % Initial speed

len = 500;
Q = linspace(-2.8e-14, 2.8e-14,len); % Posible charges
y = zeros(1,size(Q,2)); % Vertical position as a function of charge y(i) = f(Q(i))
j=1;

for i=1:(size(Q,2))
    v=[v0,0,0,0];
    ce=Q(i);
    while v(1,3)<(s+L+h)
        % Runge-Kutta 4 Integrator
        [Ex, Ey] = eField(rcx, rcy, v(1, 3), v(1, 4), q);
        k1 = he * [ce*Ex/m, ce*Ey/m, v(1,1), v(1,2)];

        [Ex, Ey] = eField(rcx, rcy, v(1, 3) + k1(3)/2, v(1, 4) + k1(4)/2, q);
        k2 = he * [ce*Ex/m, ce*Ey/m, v(1,1)+k1(1), v(1,2)+k1(2)];

        [Ex, Ey] = eField(rcx, rcy, v(1, 3) + k2(3)/2, v(1, 4) + k2(4)/2, q);
        k3 = he * [ce*Ex/m, ce*Ey/m, v(1,1)+k2(1), v(1,2)+k2(2)];

        [Ex, Ey] = eField(rcx, rcy, v(1, 3) + k3(3), v(1, 4) + k3(4), q);
        k4 = he * [ce*Ex/m, ce*Ey/m, v(1,1)+k3(1), v(1,2)+k3(2)];
        
        v(1,:) = v(1,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
    fprintf('%i / %i\n', j, len) % print status
    j = j+1;
    y(i)=v(1,4);
end

% Save results
writematrix(Q,'q.csv')
writematrix(y,'y.csv')