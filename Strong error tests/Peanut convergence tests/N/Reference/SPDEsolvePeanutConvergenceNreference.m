% this program is for comparing computed solutions to a reference function,
% with a reference N=400


function SPDEsolvePeanutConvergenceNreference(M,T,funstring,samples,samplestart,sampleend)
    % Inputs:
    % M: Number of timesteps
    % T: Final time
    % noise: Parameter q_i, i=1,...,N, in the series expansion of the noise
    % funstring: The nonlinear function entered in the form '@(x) f(x)'. For
    % example: '@(x) 0', '@(x) 1./(1+x.^2)'
    % samples: The number of realizations of the solutions which are
    % computed and over which the mean is taken for the error
    % The program will carry out realizations samplestart to sampleend
    
    % Resolution, Number of Eigenfucntions
    global reso
    global ef
    reso = 301;
    identifier = 'Peanut N';  % Appears in file names related to these tests    
    %% Generating/Importing eigenfunctions Bump function coefficients 
    mainfolder = fileparts(fileparts(fileparts(fileparts(pwd))));
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
    
    fN = str2func(funstring);
    x = linspace(0,1,reso);
    y = x;
    Nvector = [400,300,200,150,100,80,50,20,10];
    Nmax = Nvector(1);
    len = length(Nvector);
    totalsamples = sampleend-samplestart+1;
    allexact = zeros(samples,Nmax);
    allapproximate = zeros(samples,Nmax,len);
    errorarray = zeros(samples,len);
    V_start = Bumpcoefficients;  % Startingcoefficients
    totaltimesteps = M*totalsamples*len;
    progress = 0;
    qvector = 1./(linspace(1,Nmax,Nmax).^2);
    %% Solving SPDEs
    V_newreference = zeros(1,Nmax);
    randomnums = randn(M,Nmax,totalsamples);
    for sample = samplestart:sampleend
        dt = T/M;
        for index = 1:len
        V_new = zeros(1,Nmax);
        V_temp = V_start;
            for i=1:M
             tic
                for j=1:Nvector(index)
                c1 = qvector(j)*(1-exp(-2*dt*ev(j)))/(2*(ev(j)));
                R = sqrt(c1)*randomnums(i,j,sample);
                integralarray = fN(linearcombination(V_temp)).*ef(:,:,j);
                integralarray(isnan(ef(:,:,j))==1)=0;
                V_new(j) = exp(-(ev(j))*dt)*V_temp(j) + (1-exp(-ev(j)*dt))/ev(j)*simps(y,simps(x,integralarray,1),2) ...  
                         + R(1);
                if index==1 % Creating the reference solution
                    V_newreference(j) = V_new(j);
                end
                end
                disp(V_new(1:20))
                disp(V_newreference(1:20))
            V_temp = V_new;
            progress = progress + 1;
            fprintf("Time step %d of %d done. Elapsed time: %d seconds.\n",progress,totaltimesteps,toc)
            end
            allapproximate(sample,:,index) = V_new;
        end
        allexact(sample,:) = V_newreference;
    end
    for sample = 1:samples
        for index = 1:len
            comparevector = [allapproximate(sample,1:Nvector(index),index),zeros(1,Nmax-Nvector(index))];
            errorarray(sample,index) = norm(allexact(sample,:)-comparevector)^2;
            if sample==1 && index==4
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

function result = linearcombination(V)  % accepts a coefficient vector V and returns the function array
    global reso    
    global ef
    len = length(V);
    result = zeros(reso);
    for i=1:len
        result = result + V(i)*ef(:,:,i);
    end
end




    
    