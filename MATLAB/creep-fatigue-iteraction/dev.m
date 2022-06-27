function [dev_tensor] = dev(tensor)
%input: tensor
%output: deviatoric tensor
dev_tensor=tensor-1/3*trace(tensor)*eye(3);
end

