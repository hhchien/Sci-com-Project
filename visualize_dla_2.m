x = readmatrix('output_nutrients.txt');
y = readmatrix('output_bacteria.txt');

myColorMap = hsv(256);
myColorMap(1,:) = 0.9;
colormap(myColorMap);

colorbar;

axis equal;

for i = 1:256
    for j = 1:256
        if y(i,j) == 1
            x(i,j) = -0.01;
        end
    end
end

imagesc(x);
