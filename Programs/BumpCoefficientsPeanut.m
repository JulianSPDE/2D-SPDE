%This script computes coefficients for the initial condition of a bump
%function supported on [0.4,0.6]x[0.3,0.5] inside the unit circle.



function BumpCoefficientsPeanut(mainfolder)
    % Input: The main folder directory
    N = 400;
    reso = 301;
    ef = importdata(strcat(mainfolder,'/Data/PeanutEigenfunctions400.mat'));
    ef = real(ef);
    ef(isnan(ef)==1)=0;
    %bumpfun = @(r) exp(-1/(0.25-r*r));
    x = linspace(0,1,reso);
    y = x;
    [X,Y] = meshgrid(x,y);
    %[thetapol,rpol] = cart2pol(X,Y);
    bumpfunarray = zeros(reso);
    % We want the bump function supported on [0.4,0.6]x[0.3,0.5]
    
    for i1=1:reso
        for i2=1:reso
            bumpfunarray(i1,i2) = bumpfun(X(i1,i2),Y(i1,i2),0.4,0.6,0.3,0.5);
        end
    end
    
    fouriercoefficients = zeros(1,N);
    
    for i=1:N
        fouriercoefficients(i) = trapz(y,trapz(x,bumpfunarray.*ef(:,:,i),1),2);
        fprintf("Function %d of %d done.\n",i,N)
    end
    save(sprintf(strcat(mainfolder,'/Data/Bump Function coefficients N=%d, reso=%d.mat'),N,reso),'fouriercoefficients')
end
    
% General bump function supported on [x1,x2]x[y1,y2]
function z = bumpfun(x,y,x1,x2,y1,y2)
    lx = x2-x1;
    ly = y2-y1;
    midx = (x2+x1)/2;
    midy = (y2+y1)/2;
    rx = abs(x-midx);
    ry = abs(y-midy);
    r = sqrt((2/lx)^2*(x-midx)^2+(2/ly)^2*(y-midy)^2);
    %if rx>abs(lx)/2 || ry>abs(ly)/2
    if r>=1
        z = 0;
    else
        z = exp(-1/(1-r^2))/exp(-1);
    end
end
    
    
    
    