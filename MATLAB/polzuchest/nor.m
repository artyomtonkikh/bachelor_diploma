function [n] = nor(tensor)
n=sqrt(trace(tensor*tensor'));
end

