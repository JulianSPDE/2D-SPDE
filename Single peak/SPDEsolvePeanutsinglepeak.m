
function SPDEsolvePeanutsinglepeak(M,N,T,peak,fignum)
    % Inputs:
    % M: Number of timesteps
    % T: Final time
    % peak: the peak of the nonlinear bump function
    
    global reso
    global ef
    reso = 301;
    %% Generating/Importing eigenfunctions Bump function coefficients 
    mainfolder = fileparts(pwd);
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
  
    fN = @(x) exp(-10*(x-peak).^2);  % Bump function
    x = linspace(0,1,reso);
    y = x;
    ti = linspace(0,2*pi,500);
    funx = @(x) 0.06*((cos(x)+2).*(cos(x+0.6)+2).*(0.1*cos(3*x)+2))-0.1;
    funy = @(x) 0.06*(sin(x)+2).*((sin(x-0.5)+2).*(0.4*cos(2*x)+2).*(0.1*sin(4*x)+1))-0.06;
    vi = [funx(ti);funy(ti)];
    V_start = Bumpcoefficients;  % Startingcoefficients
    totaltimesteps = M;
    progress = 0;
    qvector = 1./(linspace(1,N,N).^2);
    
    %% Solving SPDEs
    
        dt = T/M;
        V_new = zeros(1,N);
        V_temp = V_start;
        randomnums = importdata('1000x400 Random numbers.mat');
        %randomnums = randn(1000,400);
            for i=1:M
             tic
                for j=1:N
                intarray = fN(linearcombination(V_temp)).*ef(:,:,j);
                c = qvector(j)*(1-exp(-2*dt*ev(j)))/(2*ev(j));
                R = c*randomnums(i,j);
                V_new(j) = exp(-(ev(j))*dt)*V_temp(j) + (1-exp(-ev(j)*dt))/ev(j)*simps(x,simps(y,intarray,2),1) ...  
                         + R;
                end
                h = figure('units','normalized','outerposition',[0 0 1 1]);
                set(h,'Visible', 'off');
                if i==M
                    Z = linearcombination(V_new);
                    for i1=1:reso
                        for i2=1:reso
                            [In,~] = inpolygon(y(i1),x(i2),vi(2,:),vi(1,:));
                            if In==0
                                Z(i1,i2) = NaN;
                            end
                        end
                    end
                end
                V_temp = V_new;
                progress = progress + 1;
                fprintf("Time step %d of %d done. Elapsed time: %d seconds.\n",progress,totaltimesteps,toc)
            end
save(sprintf(strcat('Peanut N=%d, M=%d Peak_%.3f.mat'),N,M,fignum),'Z')
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




    
    