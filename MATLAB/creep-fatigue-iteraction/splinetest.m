clc;
clear;
x=0:0.2:1;
y=sin(x);
for i=1:6
    y_interpolated(i)=linear_interpolation(x(i),x,y);
end
plot(x,y_interpolated,x,y);