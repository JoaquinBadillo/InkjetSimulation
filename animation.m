clc; clear all; close all;

he = 8e-7; % Step size (Euler's Method)

% Dimensions (meters)
L=0.5e-3;
s=0.2e-3;
d=0.5e-3;
h= 0.4e-3;

% Discrete model of distributions
N=50; % Number of charged particles
qtotal1=9e-9; % Total charge of the top electrode (Coulombs)
qtotal2=-9e-9;% Total charge of the bottom electrode (Coulombs)
q1=ones(1,N)*qtotal1/N;
q2=ones(1,N)*qtotal2/N; 
q=[q1,q2]; % charged particles vector

% X coordinates
rcx1=linspace(s,s+L,N);
rcx=[rcx1,rcx1];

% Y Coordinates
rcy1=ones(1,N)*(d/2); 
rcy2=ones(1,N)*(-d/2);
rcy=[rcy1,rcy2];

% Parameters
m = 1.13e-10; % Mass
v0=20; % Initial speed

% Vector field
meshes=25;
x=linspace(0,s+L+h,meshes);
y=linspace(-d,d,meshes);
[xmesh, ymesh]=meshgrid(x,y);
[Exmesh, Eymesh]=eFieldMesh(rcx,rcy,xmesh,ymesh,q);
r=sqrt(Exmesh.^2+Eymesh.^2);

% Charge vs Final position arrays
Q = linspace(-2.8e-14, 2.8e-14,100);
y = zeros(1,size(Q,2));

% Euler's Method
for i=1:size(Q,2)
    v=[v0,0,0,0];
    k=1;
    ce=Q(i);
    while v(k,3)<(s+L+h)
        [Ex, Ey] = eField(rcx, rcy, v(k, 3), v(k, 4), q);
        v(k+1,:) = v(k,:) + he*[ce*Ex/m,ce*Ey/m, v(k, 1), v(k, 2)];
        k=k+1;
    end 
    y(i)=v(end,4);
end

% Graph vertical position vs charge
figure
hold on
plot(Q,y)
ylim([-d,d])
ylabel("Final vertical position (m)")
xlim([Q(1),Q(end)])
xlabel("Ink droplet charge (C)")
title("Vertical Position vs Charge")

% Figure
yimp = [0.985, 1.543,0.481,-0.123,1.876,-0.513, 1.945,-0.842, 1.854, -1.386 , 1.386, -0.842, 1.854,-0.513, 1.945, -0.123,1.876, 1.543,0.481,0.985]*1.5e-4;
ximp = [-1.73,-1.6,-1.6,-1.2,-1.2,-0.8,-0.8,-0.4,-0.4,0,0,0.4,0.4,0.8,0.8,1.2,1.2,1.6,1.6,1.73]*1.5 + 5;
colors = [[255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20];
          [255, 0, 20]];
colors = colors./255; % Normalize values

% Required Charges (Useful for a real printer)
qimp=zeros(1,size(yimp,2)); 
% Paths (array of ODE solutions)
vimp={size(qimp,2)};
% Final vertical positions
ycoords = ones(1,size(ximp,2)); 

% Simulation
for i=1:size(yimp,2)
    j=find(((yimp(1,i)-3e-5)<y)&(y<(yimp(1,i)+3e-5)));
    qimp(i)=Q(j(round((1+size(j,2))/2)));
    v=[v0,0,0,0];
    k=1;
    ce=qimp(i);
    while v(k,3)<(s+L+h)
        [Ex, Ey] = eField(rcx, rcy, v(k, 3), v(k, 4), q);
        v(k+1,:) = v(k,:) + he*[ce*Ex/m,ce*Ey/m, v(k, 1), v(k, 2)];
        k=k+1;
    end 
    vimp{i}=v;
    ycoords(1,i)=v(end,4);
end

% Animation
figure
for j=1:size(yimp,2)
    hold on
    % Plot field and electrodes
    quiver(xmesh,ymesh,Exmesh./r,Eymesh./r,"Color","blue")
    plot(rcx1,rcy1,".r","MarkerSize",50)
    plot(rcx1,rcy2,".b","MarkerSize",50)

    xlim([0, 2*s+L+h])
    ylim([-d, d])

    % Piece of paper
    xline(s+L+h,"LineWidth",2)

    particle=scatter(vimp{j}(:,3),vimp{j}(:,4),40,colors(j,:),"filled");
    
    % Ink Animation
    for i=1:size(vimp{j},1)
        delete(particle);
        particle=scatter(vimp{j}(i,3),vimp{j}(i,4),40,colors(j,:),"filled");
        drawnow('limitrate')
    end
    hold off
end 

% Show results
disp("Printing...")
ximp(2,:) = NaN; % Add row of NaNs (improves performance)
ycoords(2,:) = NaN; % Add row of NaNs (improves performance)

figure
colororder(colors) % Set the color of the ink droplets
plot(ximp, ycoords,'.','MarkerSize', 18);
ylim([-4e-4, 4e-4])
xlim([0, 10])