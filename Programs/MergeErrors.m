% This script takes all error vectors of the different eigenfunctions and
% combines them into one single array.

%% Generating/Importing eigenfunctions Bump function coefficients 
mainfolder = fileparts(fileparts(fileparts(pwd)));
datafolder = strcat(mainfolder,'/Data');
programsfolder = strcat(mainfolder,'/Programs');

dinfoev = dir('*Peanut ev*.mat');
dinfoef = dir('*Peanut L2*.mat');
len = length(dinfoev);

filenametemp = dinfoev(1).name;
datatemp = importdata(filenametemp);
faces = datatemp.facevector;
facelen = length(faces);

save(strcat(datafolder,'/Facevector.mat'),'faces');


functionnum = 200;
evcomplete = zeros(functionnum,facelen);
efcomplete = zeros(functionnum,facelen);
for i=1:functionnum
    filenameev = dinfoev(i).name;
    dataev = importdata(filenameev);
    
    
    filenameef = dinfoef(i).name;
    dataef = importdata(filenameef);
    
    evs = dataev.eigenvaluedifferences;
    efs = dataef.normdifferences;
    
    evcomplete(i,:) = evs;
    efcomplete(i,:) = efs;
end
    
save(strcat(datafolder,sprintf('/ConvergenceErrorsPeanutEf%d.mat',functionnum)),'evcomplete')
save(strcat(datafolder,sprintf('/ConvergenceErrorsPeanutWavenum%d.mat',functionnum)),'efcomplete')


