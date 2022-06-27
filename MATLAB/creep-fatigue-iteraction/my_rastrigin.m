function [y] = my_rastrigin(x)
y=rastriginsfcn(x);
load('bestresult','bestvalue', 'bestparameter');
if(y < bestvalue)
    bestvalue=y;
    bestparameter=x;
    save('bestresult','bestvalue', 'bestparameter');
end
end

