function [colors, width, height] = imgToMatrix(file)
% Transform png image to RGB matrix

width = imfinfo(file).Width;
height = imfinfo(file).Height;

[I, map] = imread(file, 'png');

RGB = ind2rgb(I,map);

colors = zeros(width*height, 3);

k=1;

for i=1:height
    for j=1:width
        colors(k, :) = RGB(i, j, :);
        k = k + 1;
    end
end
return
end