function [sigma] = current_stress_table(t)
%input: time
%output: stress_current
key_time=[0 5 85 90 95 175 180 185 265 270 275 355 360];
key_stress=[0 200 200 500 200 200 0 200 200 500 200 200 0];
n=0;
for i=1:12
   if (key_time(i)<=t && t<=key_time(i+1))
      n=i;
   end
end
sigma=key_stress(n)+(t-key_time(n))/(key_time(n+1)-key_time(n))*(key_stress(n+1)-key_stress(n));
end

