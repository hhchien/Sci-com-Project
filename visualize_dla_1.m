x = readmatrix('output_nutrients.txt');
y = readmatrix('output_bacteria.txt');

myColorMap = hsv(256);
myColorMap(1,:) = 0.9;
colormap(myColorMap);

colorbar;

axis equal;

subplot(1,2,1);
imagesc(x);

subplot(1,2,2);
imagesc(y);