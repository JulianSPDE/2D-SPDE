 %% Generating/Importing eigenfunctions Bump function coefficients 
    mainfolder = fileparts(fileparts(pwd));
    datafolder = strcat(mainfolder,'/Data');
    programsfolder = strcat(mainfolder,'/Programs');
    addpath(genpath(programsfolder)); % Adding the programs folder
    
    ev1error = importdata('Peanut ev 001 6.51554236 facemax=500 reso=81.mat');
    ev1errordeluxe = importdata('Peanut ev 001 6.51554236 facemax=1000 reso=81.mat');
    ev56error = importdata('Peanut ev 056 39.53663872 facemax=500 reso=81.mat');
    ev56errordeluxe = importdata('Peanut ev 056 39.53663871 facemax=1000 reso=81.mat');
    ef1error = importdata('Peanut L2 001 6.51554236 facemax=500 reso=81.mat');
    ef1errordeluxe = importdata('Peanut L2 001 6.51554236 facemax=1000 reso=81.mat');
    ef56error = importdata('Peanut L2 056 39.53663872 facemax=500 reso=81.mat');
    ef56errordeluxe = importdata('Peanut L2 056 39.53663871 facemax=1000 reso=81.mat');
    
    hold on
    plot(ev1error.facevector,ev1error.eigenvaluedifferences)
    plot(ev1errordeluxe.facevector,ev1errordeluxe.eigenvaluedifferences)
    hold off
    set(gca,'Yscale','log')
    saveas(gcf,'Eigenvalue error comparison eigenvalue 1.png')
    clf
    
    hold on
    plot(ev56error.facevector,ev56error.eigenvaluedifferences)
    plot(ev56errordeluxe.facevector,ev56errordeluxe.eigenvaluedifferences)
    hold off
    set(gca,'Yscale','log')
    saveas(gcf,'Eigenvalue error comparison eigenvalue 56.png')
    clf
    
    hold on
    plot(ef1error.facevector,ef1error.normdifferences)
    plot(ef1errordeluxe.facevector,ef1errordeluxe.normdifferences)
    hold off
    set(gca,'Yscale','log')
    saveas(gcf,'Eigenfunction error comparison eigenfunction 1.png')
    clf
    
    hold on
    plot(ef56error.facevector,ef56error.normdifferences)
    plot(ef56errordeluxe.facevector,ef56errordeluxe.normdifferences)
    hold off
    set(gca,'Yscale','log')
    saveas(gcf,'Eigenfunction error comparison eigenfunction 56.png')
    clf
    