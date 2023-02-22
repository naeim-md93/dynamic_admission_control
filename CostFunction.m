function [cost] = CostFunction(s, A, J, V, U, W, costs, D, L, Rho, B, R)
    
    lenU = size(U, 1);
    UMax = max(max(U));
    lenW = size(W, 1);
    WMax = max(max(W));
    lenJ = size(J, 1);
    lenA = size(A, 4);
    
    cb = costs(1);
    cd = costs(2);
    co = costs(3);
    ce = costs(4);
    
    surgeryCost = 0;
    waitingCost = 0;
    bedExtraUseCost = 0;
    insufficienBedCost = 0;
    
    cost = zeros(lenA, 1);
    
    for (i=1:lenA)
        a = A(:,:,:,i);
        for (j=1:lenJ)
            for (u=1:lenU)
                for (w=1:lenW)
                    surgeryCost = surgeryCost + (V(j).*U(u,j).*w.*a(u,w,j));
                    waitingCost = waitingCost + (V(j).*U(u,j).*w.*(s(u,w,j)-a(u,w,j)));
                    bedExtraUseCost = bedExtraUseCost + max(sum(sum((a(u,w,j).*D(j) - Rho(1).*B(j)), 2), 1), 0);
                    insufficienBedCost = insufficienBedCost + sum(sum(sum(a(u,w,j).*L(j)-Rho(2).*R, 2), 1), 3);
                end
            end
        end
        
        surgeryCost = cb .* sum(sum(surgeryCost, 2), 1);
        waitingCost = cd .* sum(sum(waitingCost, 2), 1);
        bedExtraUseCost = co .* sum(bedExtraUseCost);
        insufficienBedCost = ce .* max(insufficienBedCost, 0);
        
        cost(i)= surgeryCost+waitingCost+bedExtraUseCost+insufficienBedCost;
    end
end
