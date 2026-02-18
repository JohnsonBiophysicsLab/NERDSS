%%
% fit the pucadyil data using exp_lag functon
% import the pucadyil_data
pucadyil_time_mean_std_sem = importdata(['/Users/sikao/OneDrive - Johns Hopkins/JHU/Pucadyil/ExperimentalDataAndFitP/Pucadyil_time_mean_std_sem.dat']);

% extract the time and mean
pucadyil_time = pucadyil_time_mean_std_sem.data(:,1);
pucadyil_mean = pucadyil_time_mean_std_sem.data(:,2);

% fit the data
ft = fittype('myexp_lag(x, beta1, beta2, beta3, beta4)');
f = fit( pucadyil_time, pucadyil_mean, ft, 'StartPoint', [2000, 100, 10, 300] );
plot( f, pucadyil_time, pucadyil_mean )

%%
% plot exp
offset = 300;
scale = 1.6;
pucadyil_time_mean_std_sem = importdata(['/Users/sikao/OneDrive - Johns Hopkins/JHU/Pucadyil/ExperimentalDataAndFitP/Pucadyil_time_mean_std_sem.dat']);
% extract the time and mean
xexp = pucadyil_time_mean_std_sem.data(:,1);
yexp = (pucadyil_time_mean_std_sem.data(:,2) - offset) .* scale;

plot(xexp,yexp,'b','LineWidth',2)
hold on

%%
% plot the average or each trace (numbers on membrane or solution, need to treat the hist*.dat first),
% the data should be in subdirs from 1 to tracenum
dirname = '/Users/sikao/Downloads/Simulations/0911/pair/data/';
filename1 = '/sim_time.dat';
filename2 = '/sim_clatMem.dat';
filename3 = '/sim_clatSol.dat';
tracenum = 5;
meanoreach = 2; % 2 for mean, 3 for each
minmumlength = 0; % change 0 to 10000 to select the traces longer than 1000
%find minimum time: time_length
time_length = inf;
% for i = 1 : tracenum
for i = [1,3,4,5,6,7,8]
    species_temp = importdata([dirname, num2str(i), filename1]);
    if length(species_temp)-1 < time_length && length(species_temp)-1 > minmumlength
        time_length = length(species_temp)-1;
    end
end
%adjust all traces to time_length
species_all = [];
species_sol = [];
% for i = 1 : tracenum
for i = [1,3,4,5,6,7,8]
    time = importdata([dirname, num2str(i), filename1]);
    clatMem = importdata([dirname, num2str(i), filename2]);
    clatSol = importdata([dirname, num2str(i), filename3]);
    if length(time)-1 >= time_length
        eval(['clatMem_', num2str(i), '=', 'clatMem;']);
        eval(['clatSol_', num2str(i), '=', 'clatSol;']);
        species_all = [species_all, clatMem(1:time_length)];
        species_sol = [species_sol, clatSol(1:time_length)];
        x = time(1:time_length);
    end
end

y_mean = mean(species_all, meanoreach);
y_std = std(species_all, 0, 2);

ysol_mean = mean(species_sol, meanoreach);
ysol_std = std(species_sol, 0, 2);

% % fit the data
% ft = fittype('myexp_lag(x, beta1, beta2, beta3, beta4)');
% f = fit( x, y_mean, ft, 'StartPoint', [2000, 100, 10, 10] );
% plot( f, x, y_mean )

plot(x, y_mean, 'g', 'LineWidth',2)

