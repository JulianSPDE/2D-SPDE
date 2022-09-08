% This script imports all vectors containing errors from different samples
% and outputs their average. It also creates an error plot.

function averageerrors(identifier,totalsamples)
% Input: Depending on the parameter and shape, the identifier is 'Disk M',
% 'Disk N', 'Disk reso', 'Peanut M' or 'Peanut N'.
% Number of total samples, the number that the error vector is divided by
% to average it

errorfiles = dir(strcat('*',identifier,' errors*.mat'));
len = length(errorfiles);
xvector = importdata('xvector.mat');
xlen = length(xvector);
errors = zeros(1,xlen);
for i=1:len
    filename = errorfiles(i).name;
    error = importdata(filename);
    errors = errors + error;
end
errors = errors/totalsamples;
save(strcat('Total errors ',identifier,'.mat'),'errors','xvector')
    errors(1)=[];
    xvector(1)=[];
    fig = figure(1);
    set(fig,'Visible','off');
    plot(xvector,errors)
    set(gca,'Yscale','log')
    xlabel('Tested quantity')
    ylabel('L2 error')
    saveas(gcf,sprintf(strcat(identifier,' errors.png')))
end
