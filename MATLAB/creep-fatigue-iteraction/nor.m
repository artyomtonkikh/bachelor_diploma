function [n] = nor(tensor)
%input: tensor
%output: norm of tensor
n=sqrt(trace(tensor*tensor'));
end

