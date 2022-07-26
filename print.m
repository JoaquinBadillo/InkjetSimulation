clc; clear all; close all;

% Load data (computed previously with Runge-Kutta 4)
Q = readmatrix('./q.csv'); % Array of linearly spaced charges
y = readmatrix('./y.csv'); % Array of corresponding final vertical positions

% Plot the final vertical position vs charge
figure
plot(Q,y)
ylim([-6e-4,6e-4])
ylabel("Final vertical position (m)")
xlim([Q(1),Q(end)])
xlabel("Ink droplet charge (C)")
title("Vertical Position vs Charge")

% get RGB array and dimensions
[colors, width, height] = imgToMatrix('image.png');

ylin=linspace(4,-4,height);
xlin=linspace(0,10,width);
[x1, y1]=meshgrid(xlin,ylin);

qimp = zeros(size(x1));

ximp = reshape(x1.',1,[]);
yimp = reshape(y1.',1,[]);

yimp=yimp*1e-4; % Adjust the size
 
ycoords = zeros(1,size(ximp,2)); 

for i=1:size(yimp,2)
    j=find(((yimp(1,i)-3e-5)<y)&(y<(yimp(1,i)+3e-5)));
    j = j(round((1+size(j,2))/2)); % use midpoint
    qimp(i)=Q(j);
    ycoords(i)=y(j);
end

disp("Matrix of charges saved to 'charges.csv'.") % This is what a real printer might need
writematrix(qimp,'charges.csv') 

% Print
disp("Printing...")
ximp(2,:) = NaN; % Add row of NaNs
ycoords(2,:) = NaN; % ADd row of NaNs

figure
colororder(colors) % Set the color of the ink droplets
plot(ximp, ycoords,'.','MarkerSize', 18);
title('Print')
xticks([])
yticks([])
xlim([0, 10])