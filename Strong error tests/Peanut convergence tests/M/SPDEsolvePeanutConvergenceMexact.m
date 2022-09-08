
function SPDEsolvePeanutConvergenceMexact(N,T,samples,samplestart,sampleend)
    % Inputs:
    % M: Number of timesteps
    % T: Final time
    % samples: The number of realizations of the solutions which are
    % computed and over which the mean is taken for the error
    % The program computes samples samplestart to sampleend
    
    global reso
    global ef
    reso = 301;
    identifier = 'Peanut M';  % Appears in file names related to these tests    
    %% Generating/Importing eigenfunctions Bump function coefficients 
    mainfolder = fileparts(fileparts(fileparts(pwd)));
    datafolder = strcat(mainfolder,'/Data');
    programsfolder = strcat(mainfolder,'/Programs');
    addpath(genpath(programsfolder)); % Adding the programs folder
    ef = importdata(strcat(datafolder,'/PeanutEigenfunctions400.mat')); % Generating necessary data in the data folder
    ef = real(ef);
    ef(isnan(ef)==1)=0;    
    wavenum = importdata(strcat(datafolder,'/PeanutWavenumbers.mat'));
    wavenum = real(wavenum);
    ev = wavenum.^2;    
    Bumpcoefficients = importdata(sprintf(strcat(datafolder,'/Bump Function coefficients N=400, reso=301.mat')));
    
    
    %% Initializing variables

    Mvector = [1000,500,200,100,50,20,10];
    len = length(Mvector);
    totalsamples = sampleend-samplestart+1;
    errorarray = zeros(samples,len);
    V_start = Bumpcoefficients;  % Startingcoefficients
    V_startexact = Bumpcoefficients;
    totaltimesteps =sum(Mvector)*totalsamples;
    progress = 0;
    qvector = 1./(linspace(1,N,N).^4);
    %% Solving SPDEs
    
    for sample = samplestart:sampleend
        for index = 1:len
            V_new = zeros(1,N);
            V_newexact = zeros(1,N);
            V_temp = V_start;
            V_tempexact = V_startexact;
            dt = T/Mvector(index);
            for i=1:Mvector(index)
                tic
                for j=1:N
                c1 = qvector(j)*(1-exp(-2*dt*(ev(j)+0.5)))/(2*(ev(j)+0.5));
                c3 = qvector(j)*(1-exp(-2*dt*ev(j)))/(2*ev(j));
                R = [sqrt(c1);sqrt(c3)]*randn;
                V_new(j) = exp(-(ev(j))*dt)*V_temp(j) + (1-exp(-ev(j)*dt))/ev(j)*(-0.5)*V_temp(j) ...  
                         + R(1);
                V_newexact(j) = exp(-(ev(j)+0.5)*dt)*V_tempexact(j) ...
                             + R(2);
                end
            V_temp = V_new;
            V_tempexact = V_newexact;
            progress = progress + 1;
            fprintf("Time step %d of %d done. Elapsed time: %d seconds.\n",progress,totaltimesteps,toc)
            end
        errorarray(sample,index) = errorarray(sample,index) + norm(V_new-V_newexact)^2;
        end
    end
    errorvector = zeros(1,len);
    for sample = 1:samples
        errorvector = errorvector + errorarray(sample,:);
    end
    errorvector = sqrt(errorvector/samples);
    errorvector(1)=[];
    Mvector(1)=[];
    figure(1)
    plot(Mvector,errorvector)
    set(gca,'Yscale','log')
    set(gca,'Xscale','log')
    save(strcat(identifier,' errors ',num2str(samples),sprintf('Samples %d to %d.mat',samplestart,sampleend)),'errorvector')
    save('xvector.mat','Mvector')
end
