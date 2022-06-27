strain_rate=dlmread('C:\Users\artyo\Desktop\data\creep rate 190 MPa.txt');
cyclesArray = log10(strain_rate(:,1));
strainArray = log10(strain_rate(:,2));
plot(cyclesArray,strainArray);
n = 7; % polynomial degree (you can change it as you wish)
p = polyfit(cyclesArray,strainArray,n); % p is coefficient of your polynomial: P(X) = P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1) descending order.
newY = polyval(p,cyclesArray); % function results
plot(cyclesArray,strainArray, 'bo',cyclesArray,newY,'-r');
grid;
legend('Data','Fitted Data');
