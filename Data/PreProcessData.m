%% Pre-Process Data from Lee and Moser (2015)

ReName='0180';
% ReName='2000';
% ReName='5200';

[y,yplus,U]=importmean(['./LM_Channel_' ReName '_mean_prof.dat']);
U_std=importstd(['./LM_Channel_' ReName '_mean_stdev.dat']);

Re = importRe(['./LM_Channel_' ReName '_mean_prof.dat']);

save(['LeeMoser2015Channel' ReName 'Data.mat'],'-v7.3');