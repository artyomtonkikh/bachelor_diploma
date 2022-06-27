rng default % For reproducibility
options = optimoptions('ga','PlotFcn','gaplotbestf', 'FitnessLimit', 0)
LB = 0*[1 1];
UB = 4*[1 1];
bestvalue=1000000;
bestparameter=[1 1];
save('bestresult','bestvalue', 'bestparameter');
[solution,objectiveValue] = ga(@my_rastrigin,2,...
    [],[],[],[],LB, UB,[] ,options)