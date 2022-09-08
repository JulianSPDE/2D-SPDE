alldata = dir('Peanut*.mat');
num = length(alldata);

maxvalue = 0;
minvalue = 0;
for i=1:num
    filename = alldata(i).name;
    Z = importdata(filename);
    Z = real(Z);
    if i==1
        disp(filename)
    end
    Z(isnan(Z))=0;
    Zmin = min(min(Z));
    Zmax = max(max(Z));
    maxvalue = max(maxvalue,Zmax)
    minvalue = min(minvalue,Zmin)
end
%save('values.mat','minvalue','maxvalue')
for i=1:num
    filename = alldata(i).name;
    savename = num2str(i);
    plot_Peanut_withboundary(filename,savename,minvalue,maxvalue)
end

