% this is a script which reads all documents in the folder with convergence
% information and plots how many boundary elements are needed to stay below
% an input error tolerance.

function plotnecessaryboundaryelements(evtol,eftol)
    %% Generating/Importing eigenfunctions Bump function coefficients 
    mainfolder = fileparts(fileparts(pwd));
    datafolder = strcat(mainfolder,'/Data');
    programsfolder = strcat(mainfolder,'/Programs');
    addpath(genpath(programsfolder)); % Adding the programs folder

    everrors = importdata(strcat(datafolder,'/ConvergenceErrorsPeanutWavenum200.mat'));
    eferrors = importdata(strcat(datafolder,'/ConvergenceErrorsPeanutEf200.mat'));
    wavenum = importdata(strcat(datafolder,'/PeanutWavenumbers.mat'));
    wavenum = real(wavenum);
    faces = importdata(strcat(datafolder,'/Facevector.mat'));
    lenx = size(everrors,1);
    leny = size(everrors,2);
    for i=1:lenx
        for j=1:leny
        everrors(i,j) = 2*wavenum(i)*everrors(i,j)+everrors(i,j)^2; % Considering eigenvalue error, not wavenumber error
        end    
    end
    ev = wavenum.^2;

    len = size(everrors,2);
    xvector = zeros(1,len);
    yvector = xvector;

    evnumber = size(everrors,1);
    c = ones(evnumber,3);  % Color vector
    for i=1:evnumber
        c(i,:) = [0,1,0];  % Green color
    end
    for i=1:evnumber
        everrorvector = everrors(i,:);
        eferrorvector = eferrors(i,:);
        [indexev,boolev] = largestunderthreshold(everrorvector,evtol);
        [indexef,boolef] = largestunderthreshold(eferrorvector,eftol);
        index = min(indexev,indexef);
        fprintf("indexev: %d\n",indexev)
        fprintf("indexef: %d\n",indexef)
        % xvector(i) = i;      % Choice: Either the plain number of the
        % wavenumber, or the wavenumber itself
        xvector(i) = ev(i);
        if index == 0
            yvector(i) = 0;
        else
            yvector(i) = faces(index);
        end
        if boolev == true || boolef == true
            c(i,:) = [1,0,0];   % Color red; Needed accuracy cannot be achieved here
        end
    end

    xvector = sqrt(xvector);  % We want to plot dependence on wavenumber
    figure(1)
    scatter(xvector,yvector,[],c)
    axis([0,80,100,500])
    xlabel('Wavenumber')
    ylabel('Needed boundary elements')
    title(sprintf('Evtol: %.1e, eftol: %.1e',evtol,eftol))
    saveas(gca,sprintf('Needed boundary elements evtol %.1e eftol %.1e.png',evtol,eftol))
    save(sprintf('Needed BoundaryElements evtol=%.1d eftol=%.1d.mat',evtol,eftol),'xvector','yvector','c')
end


% This function accepts a vector and a threshold and returns the last index
% to be under that threshold
function [index,toolarge] = largestunderthreshold(v,threshold)
    len = length(v);
    index = 0;
    toolarge = false;
    for i=2:len
        if v(i)>threshold % At this point, value is too large
            index = i-1;  % => Pick previous index
            break
        end
    end
    if (v<threshold) == ones(1,len)
        index = len;
    end
    if (v>threshold) == ones(1,len)
        index = 1;
        toolarge = true;  % We don't manage to reach this threshold; not enough boundary elements
    end
end
    
    
    
    


