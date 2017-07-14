utils = dlmread('results.dat');
f = figure('visible','off');
plot(utils)
print res/utils.jpg;
