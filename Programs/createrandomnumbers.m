


function createrandomnumbers(M,N,samples,randomnumbersfilename)
% INPUTS:
% N number of eigenfunctions
% M number of time steps
% samples number of realizations
% name: string with the name of the random number array.

% Different identifiers:
% 'Disk N', 'Disk M', 'Disk reso', 'Peanut M', 'Peanut N'

    randinfo = dir(randomnumbersfilename);
    isinfolder = length(randinfo);
    if isinfolder == 0
        randomnumbers = randn(M,N,samples);
        save(randomnumbersfilename,'randomnumbers')
    end
end