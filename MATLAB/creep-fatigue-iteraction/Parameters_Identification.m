er=@(y)Error_3_experiments(y);
bestvalue=10000000;
save('bestresult','bestvalue');
y0=[1 1 1 1 1];
y=fminsearch(er,y0)