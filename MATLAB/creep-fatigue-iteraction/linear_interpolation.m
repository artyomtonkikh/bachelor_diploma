function [funval] = linear_interpolation(x,key_points,key_values)
N=length(key_points)-1;%количество интервалов
for i=1:N
    if(x>=key_points(i) && x<=key_points(i+1))
       funval=key_values(i)+(x-key_points(i))*(key_values(i+1)-key_values(i))/(key_points(i+1)-key_points(i));
       return;
    end
end

