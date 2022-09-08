% This script checks, for a given resolution and eigenfunction number, if
% the relevant basis is already in the designated folder. If it is, then it
% is imported straightaway, if it is not, then it is created in that folder
% and imported afterwards.

function [ef,bumpcoefficients] = importdiskfunctions(reso,N,datafolder)

    effilename = sprintf(strcat(datafolder,"/DiskEigenfunctions N=%d reso=%d.mat"),N,reso);   
    bumpfilename = sprintf(strcat(datafolder,"/Bump Function coefficients N=%d, reso=%d.mat"),N,reso);
    efinfo = dir(effilename);
    isinfolder = length(efinfo);
    bumpinfo = dir(bumpfilename);

    if isinfolder == 0
        fprintf("Creating Base functions.\n")
        CreateBesselArray(N,reso,datafolder);
    end
    
    isinfolder = length(bumpinfo);    
    if isinfolder == 0
        fprintf("Creating Bump function coefficients.\n")
        BumpCoefficients(N,reso,datafolder);
    end
    ef = importdata(effilename);
    bumpcoefficients = importdata(bumpfilename);
    
    
end