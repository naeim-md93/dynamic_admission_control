function [V] = LinearFunction(Phi, Theta)
  % size(Phi) = (len_feature_vector, 1) for state=1
  % size(Theta) = (1, parameters) for state=1
  V = transpose(Phi) * Theta;
end