function [a, aMinIndex] = ComputeActions(s, AStar, J, V, U, W, costs, D, L, Rho, B, R, phiNext, thetaPrev, gamma)
    cost = CostFunction(s, AStar, J, V, U, W, costs, D, L, Rho, B, R);
    
    [~, aMinIndex] = min(cost+(gamma.*transpose(phiNext)*thetaPrev));
    a = AStar(:,:,:,aMinIndex);
end
