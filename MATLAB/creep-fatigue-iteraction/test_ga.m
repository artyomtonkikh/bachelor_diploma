rng default % For reproducibility
options = optimoptions('ga','PlotFcn','gaplotbestf', 'FitnessLimit', 0)
LB = 0*[1 1 1 1 1 1 1 1 1 1];
UB = 100*[1 1 1 1 1 1 1 1 1 1];
bestvalue=1000000;
save('bestresult','bestvalue');
[solution,objectiveValue] = ga(@Error_3_experiments,10,...
    [],[],[],[],LB, UB,[] ,options)