xa = readmatrix('output_simulate.txt');
y = readmatrix('output_bacteria.txt');

myColorMap = hsv(256);
myColorMap(1,:) = 0.9;
colormap(myColorMap);

for i = 1:size(xa,1)/256
    xi = xa((i-1)*256+1:256*i,1:256);
    imagesc(xi);
    pause(0.001);
end

imagesc(y);