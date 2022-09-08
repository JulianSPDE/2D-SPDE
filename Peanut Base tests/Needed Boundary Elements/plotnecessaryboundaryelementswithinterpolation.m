% This is a script which reads all documents in the folder with convergence
% information and plots how many boundary elements are needed to stay below
% an input error tolerance.

function plotnecessaryboundaryelementswithinterpolation(evtol,eftol)
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
    xvector2 = zeros(1,len);
    yvector = xvector;

    evnumber = size(everrors,1);
    c = ones(evnumber,3);  % Color vector
    for i=1:evnumber
        c(i,:) = [0,1,0];  % Green color
    end
    for i=1:evnumber
        everrorvector = everrors(i,:);
        eferrorvector = eferrors(i,:);
        [indexev,evtoolarge,~] = largestunderthreshold(everrorvector,evtol);
        [indexef,eftoolarge,~] = largestunderthreshold(eferrorvector,eftol);
        xvector(i) = ev(i);
        xvector2(i) = i;
        if indexev == 1 || indexef == 1
            yvector(i) = faces(1);
        end
        if indexef<len && indexev<len && evtoolarge == false && eftoolarge == false
            evinterpolation = ((evtol-everrorvector(indexev))/(everrorvector(indexev+1)-everrorvector(indexev)))*(faces(indexev+1)-faces(indexev))+faces(indexev);
            efinterpolation = ((eftol-eferrorvector(indexef))/(eferrorvector(indexef+1)-eferrorvector(indexef)))*(faces(indexef+1)-faces(indexef))+faces(indexef);
            yvector(i) = max(evinterpolation,efinterpolation);
        end
        if evtoolarge == true || eftoolarge == true
            yvector(i) = faces(1);
            c(i,:) = [1,0,0];   % Color red; Needed accuracy cannot be achieved here
        end
        if i==2
            disp(indexev)
            disp(indexef)
            fprintf("Entry:\n")
            disp(yvector(i))
        end
    end
    find(yvector==max(yvector))

    yes = 1; % Excluding the first value in plots
    if yes
        xvector(1) = [];
        yvector(1) = [];
        c(1,:) = [];
    end

    xvector = sqrt(xvector); % Wavenumbers
    figure(1)
    size(xvector)
    size(yvector)
    size(c)
    scatter(xvector,yvector,[],c)
    xlabel('Wavenumber')
    ylabel('Needed boundary elements')
    title(sprintf('Evtol: %.1e, eftol: %.1e',evtol,eftol))
    saveas(gca,sprintf('Needed boundary elements interpolated evtol %.1e eftol %.1e.png',evtol,eftol))
    save(sprintf('Needed boundary elements interpolated evtol %.1e eftol %.1e.mat',evtol,eftol),'xvector','yvector','c')
end


% This function accepts a vector and a threshold and returns the last index
% to be under that threshold
function [index,toolarge,toosmall] = largestunderthreshold(v,threshold)
    len = length(v);
    index = 0;
    toolarge = false;
    toosmall = false;
    for i=2:len
        if v(i)>threshold % At this point, value is too large
            index = i-1;  % => Pick previous index
            break
        end
    end
    if (v<threshold) == ones(1,len)
        index = len;
        toosmall = true;
    end
    if (v>threshold) == ones(1,len)
        index = 1;
        toolarge = true;  % We don't manage to reach this threshold; not enough boundary elements
    end
end
    
    
    
    


