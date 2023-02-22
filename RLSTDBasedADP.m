function [aTau, zTau, pTau, thetaTau] = RLSTDBasedADP(sTau, N, zPrev, pPrev, thetaPrev, J, V, U, W, costs, D, L, Rho, B, R, gamma, lambda, epsilon)
    
    s1 = sTau;
    z0 = zPrev;
    p0 = pPrev;
    theta0 = thetaPrev;

    lenU = size(U, 1);
    UMax = max(max(U));
    lenW = size(W, 1);
    WMax = max(max(W));
    lenJ = size(J, 2);
    Xi = sum(sum(W));
    
    Phi = poissrnd(ones(Xi, N+1)*5);
    
    zN = z0;
    pN = p0;
    thetaN = theta0;
    
    while (true)
        for (n=1:N)
            phiN = Phi(:,n);
            phiNext = Phi(:,n+1);
            
            AStar = ADP(s1, J, V, U, W, costs, D, L);
            lenI = size(AStar, 4);
            
            nTilda = poissrnd(ones(UMax, lenJ, lenI)*5);
            
            Psi = poissrnd(nTilda);
            
            [a, aMinIndex] = ComputeActions(s1, AStar, J, V, U, W, costs, D, L, Rho, B, R, phiNext, thetaN, gamma);
            
            sNew = zeros(size(s1));
            G = s1 - a;
            for (j=1:lenJ)
                sNew(:,1,j) = Psi(:,j,aMinIndex);
                sNew(:,2:end,j) = G(:,1:end-1,j);
            end
            
            WMask = MaskMaker(sNew, U, W, 0);
            s1 = sNew .* WMask;

            PartA = transpose(phiN - (gamma .* phiNext));

            e = CostFunction(s1, a, J, V, U, W, costs, D, L, Rho, B, R) - (PartA * theta0);
    
            zN = (gamma * lambda * zN) + phiN;
            
            PartB = 1 + (PartA * pN * zN);
            
            pN = pN - ((pN * zN * PartA * pN) ./ (PartB));
            
            thetaN = thetaN + ((pN * zN * e) ./ (PartB));
        end
        norm((thetaN - theta0) ./ theta0)
        if (norm((thetaN - theta0) ./ theta0) < epsilon)
            break
        end
        
        s1 = sTau;
        z0 = zN;
        p0 = pN;
        theta0 = thetaN;
    end
 
    AStar = ADP(sTau, J, V, U, W, costs, D, L);
    zTau = zN;
    pTau = pN;
    thetaTau = thetaN;
    aTau = ComputeActions(s1, AStar, J, V, U, W, costs, D, L, Rho, B, R, phiNext, thetaTau, gamma);
end
