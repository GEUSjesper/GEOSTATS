function loopprogress(description,step,step_max)
% Must start tic before loop
percentage = linspace(0,step_max,100);

idx = find(step <= percentage,1);
clc
disp([description,':'])
disp([num2str(idx),' % Done'])
disp(['Estimated looping time remaining: ', num2str(round(toc*(step_max/step-1))), ' seconds'])
