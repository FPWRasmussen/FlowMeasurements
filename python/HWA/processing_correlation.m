clc; clear; close all;

%%

fileIn = 'CorrelationTest';

delimiter = ' ';
startRow = 23;
formatSpec = '%s';
try
    fileID = fopen(fileIn,'r');
catch
    fileID = fopen(fileIn{1},'r');
end
tmp = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
tmp2 = strrep(tmp{1},',','.');

tmp2 = str2mat(tmp2);
t = tmp2(:,1:8);
t = str2num(t);
u= tmp2(:,10:17);
u = str2num(u);

meanvoltage = mean(u);
dev = u-meanvoltage;

corr = xcorr(dev);
corr2 = corr(100000:end);
corr3 = corr2/corr2(1);
plot(corr3(1:1000))