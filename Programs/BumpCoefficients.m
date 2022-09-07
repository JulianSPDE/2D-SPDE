%This script computes coefficients for the initial condition of a bump
%function supported on [-0.5,0.5]x[-0.5,0.5] inside the unit circle.



function BumpCoefficients(N,reso,datafolder)
    %N = 200;
    %reso = 200; %Depending on input
    ef = importdata(sprintf(strcat(datafolder,'/DiskEigenfunctions N=%d reso=%d.mat'),N,reso));
    %bumpfun = @(r) exp(-1/(0.25-r*r));
    x = linspace(-1,1,reso);
    y = x;
    [X,Y] = meshgrid(x,y);
    [~,rpol] = cart2pol(X,Y);
    bumpfunarray = zeros(reso);
    for i1=1:reso
        for i2=1:reso
            bumpfunarray(i1,i2) = bumpfun(rpol(i1,i2));
        end
    end
    
    fouriercoefficients = zeros(1,N);
    
    for i=1:N
        fouriercoefficients(i) = trapz(y,trapz(x,bumpfunarray.*ef(:,:,i),1),2);
    end
    
    save(sprintf(strcat(datafolder,'/Bump Function coefficients N=%d, reso=%d.mat'),N,reso),'fouriercoefficients')
    
end
    

function z = bumpfun(r)
    if r<0.5
        z = exp(-1/(0.25-r*r))/exp(-4);
    else
        z = 0;
    end
end
    
    
    
    