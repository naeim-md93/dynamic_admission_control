function [Theta] = GradientDescent(Theta_prev, Eta, E, lambda, gamma, S)
    n = length(Theta_prev);
    Theta = 0;
    for (k=1:n)
        Theta += (lambda .* gamma)^(n-k) .* LinearFunction(S, Theta_prev)
    end
    
    Theta += Theta_prev + (Eta .* E)
end