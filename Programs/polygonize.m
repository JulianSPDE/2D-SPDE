%This function accepts an array that should be plotted,
%and cuts out a specified polygonal shape to the specified resolution

function result = polygonize(z,res,contour) %Contour: Array with 2 rows, x and y
x1 = linspace(-1,1,res);
x2 = x1;
result = z;
[xx,yy] = meshgrid(x1,x2);
for i =  1:res
  for j = 1:res
    x = xx(i,j);
    y = yy(i,j);
    [In,~] = inpolygon(x,y,contour(1,:),contour(2,:));
    if In == 0
     result(i,j) = NaN;  %Maybe later try to isnert white color
    end
  end
end
end
