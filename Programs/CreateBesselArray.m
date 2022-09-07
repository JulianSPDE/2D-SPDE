% This script creates a three-dimensional array containing arrays of the
% first eigenfunctions of the disk. 
% This will be a reso x reso x N array.
% For simplicity, we will do the unit disk, so the domain is [-1,1]x[-1,1].

function CreateBesselArray(N,reso,datafolder)
    % We start by fixing the order of the eigenfunctions by creating a Nx3
    % array containing all the indices in order:
  
  Besselarray = zeros(reso,reso,N);  
    
  % Creating templates for the contour
  param = 2*pi*linspace(0,1,reso);
  circ = [cos(param);sin(param)];    %Circle contour
  Z = zeros(reso);
  Ztemp = polygonize(Z,reso,circ);
  nanarray = isnan(Ztemp);              %Array with 0 inside the disk and 1 outside the disk
  %Z = NaN*zeros(reso);
  %fun = @(r,theta) Vkcircle(V_start1,V_start2,integralnormarray1,integralnormarray2,r,theta,besselzeroarray);
  x = linspace(-1,1,reso);
  y = x;
  [X,Y] = meshgrid(x,y);
  [thetapol,rpol] = cart2pol(X,Y);   
   
    wavenumarray = besselzero((0:150)',150,1);
    wavenumvector = reshape(wavenumarray,1,22650);  % All wavenumbers
    wavenumvector = sort(wavenumvector);
    currentwavenum = 0;
    %wavenumwithmultiplicities = zeros(1,2000);
    %index = 1;
    i = 1;
    while i<=N
        currentwavenum = currentwavenum + 1;
        [index1,index2] = find(wavenumarray==wavenumvector(currentwavenum));
        if index1 == 1
            %wavenumwithmultiplicities(index) = wavenumvector(currentwavenum);
            %index = index + 1;
            fun = @(r,theta) besselj(index1-1,wavenumarray(index1,index2).*r);
            % Radially symmetric function: Simple eigenvalue, one
            % eigenfunction
            for i1=1:reso
                for i2=1:reso
                    if nanarray(i1,i2)==0
                        Besselarray(i1,i2,i) = fun(rpol(i1,i2),thetapol(i1,i2));
                    end    
                end
            end
            fprintf("Function %d done.\n",i)
            i = i+1;
        else   % Non-radially symmetric function: Two eigenfunctions
            %wavenumwithmultiplicities(index) = wavenumvector(currentwavenum);
            %wavenumwithmultiplicities(index+1) = wavenumvector(currentwavenum);
            %index = index + 2;
            fun1 = @(r,theta) besselj(index1-1,wavenumarray(index1,index2).*r).*cos((index1-1)*theta);
            for i1=1:reso
                for i2=1:reso
                    if nanarray(i1,i2)==0
                        Besselarray(i1,i2,i) = fun1(rpol(i1,i2),thetapol(i1,i2));
                    end
                end
            end
            fprintf("Function %d done.\n",i)
            i = i+1;
            if i<N
            fun2 = @(r,theta) besselj(index1-1,wavenumarray(index1,index2).*r).*sin((index1-1)*theta);
            for i1=1:reso
                for i2=1:reso
                    if nanarray(i1,i2)==0
                        Besselarray(i1,i2,i) = fun2(rpol(i1,i2),thetapol(i1,i2));
                    end
                end
            end
            fprintf("Function %d done.\n",i)
            i = i+1;
            end
        end 
    end
    
    % Norming the plot arrays
    
    for i=1:N
        Z = Besselarray(:,:,i);
        norm = sqrt(trapz(y,trapz(x,Z.*Z,1),2));
        Besselarray(:,:,i) = Z/norm;
    end
    
   
    %save('DiskWavenumbers2000new.mat','wavenumwithmultiplicities');
    
    save(sprintf(strcat(datafolder,'/DiskEigenfunctions N=%d reso=%d.mat'),N,reso),'Besselarray')

    % We go through the eigenfunctions u_mnl by alternating l=0 and l=1
    % whenever possible 
    

end

