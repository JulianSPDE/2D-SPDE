
function SPDEsolvePeanutConvergenceNexact(M,T,samples,samplestart,sampleend)
    % Inputs:
    % M: Number of timesteps
    % T: Final time
    % samples: The number of realizations of the solutions which are
    % computed and over which the mean is taken for the error
    % The program will carry out realizations samplestart to sampleend    

    global reso
    global ef
    reso = 301;
    identifier = 'Peanut N';  % Appears in file names related to these tests    
    %% Generating/Importing eigenfunctions Bump function coefficients 
    mainfolder = fileparts(fileparts(fileparts(pwd)));
    datafolder = strcat(mainfolder,'/Data');
    programsfolder = strcat(mainfolder,'/Programs');
    addpath(genpath(programsfolder)); % Adding the programs folder  
    wavenum = importdata(strcat(datafolder,'/PeanutWavenumbers.mat'));
    wavenum = real(wavenum);
    ev = wavenum.^2;    
    Bumpcoefficients = importdata(sprintf(strcat(datafolder,'/Bump Function coefficients N=400, reso=301.mat')));
    perturbation = 0;
    evperturbed = ev + perturbation;
    %% Initializing variables
    
    Nvector = [400,300,200,150,100,80,50,20,10];
    Nmax = Nvector(1);
    len = length(Nvector);
    totalsamples = sampleend-samplestart+1;
    allexact = zeros(samples,Nmax);
    allapproximate = zeros(samples,Nmax);
    errorarray = zeros(samples,len);
    V_start = Bumpcoefficients;  % Startingcoefficients
    V_startexact = Bumpcoefficients;
    totaltimesteps = M*totalsamples;
    progress = 0;
    qvector = 1./(linspace(1,Nmax,Nmax).^2);
    %% Solving SPDEs
    
    for sample = samplestart:sampleend
        dt = T/M;
        V_new = zeros(1,Nmax);
        V_newexact = zeros(1,Nmax);
        V_temp = V_start;
        V_tempexact = V_startexact;       
            for i=1:M
             tic
                for j=1:Nmax
                c1 = qvector(j)*(1-exp(-2*dt*(ev(j)+0.5)))/(2*(ev(j)+0.5));
                c3 = qvector(j)*(1-exp(-2*dt*ev(j)))/(2*ev(j));
                R = [sqrt(c1);sqrt(c3)]*randn;
                V_new(j) = exp(-(evperturbed(j))*dt)*V_temp(j) + (1-exp(-(evperturbed(j))*dt))/evperturbed(j)*(-0.5)*V_temp(j) ...  
                         + R(1);
                V_newexact(j) = exp(-(ev(j)+0.5)*dt)*V_tempexact(j) ...
                             + R(2);
                end
            V_temp = V_new;
            V_tempexact = V_newexact;
            progress = progress + 1;
            fprintf("Time step %d of %d done. Elapsed time: %d seconds.\n",progress,totaltimesteps,toc)
            end
        errorarray(sample,:) = errorarray(sample,:) + norm(V_new-V_newexact)^2;
        allexact(sample,:) = V_newexact;
        allapproximate(sample,:) = V_new;
    end
    for sample = 1:samples
        for index = 1:len
            comparevector = [allapproximate(sample,1:Nvector(index)),zeros(1,Nmax-Nvector(index))];
            errorarray(sample,index) = norm(allexact(sample,:)-comparevector)^2;
            if sample==1 && index==3
                disp([comparevector',allexact(sample,:)'])
            end
        end
    end
    errorvector = zeros(1,len);
    for sample = 1:samples
        errorvector = errorvector + errorarray(sample,:);
    end
    errorvector = sqrt(errorvector/samples);
    errorvector(1)=[];
    Nvector(1)=[];
    figure(1)
    plot(Nvector,errorvector)
    set(gca,'Yscale','log')
    set(gca,'Xscale','log')
    save(strcat(identifier,' errors ',num2str(samples),sprintf('Samples %d to %d.mat',samplestart,sampleend)),'errorvector')
    save('xvector.mat','Nvector')
end


    
    