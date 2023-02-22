function [AStar] = ADP(s, J, V, U, W, costs, D, L)
    lenU = size(U, 1);
    UMax = max(max(U));
    lenW = size(W, 1);
    WMax = max(max(W));
    lenJ = size(s, 3);
    
    cb = costs(1);
    cd = costs(2);
    co = costs(3);
    ce = costs(4);
    
    sZeros = zeros(size(s));
    uwZeros = zeros(UMax, WMax);
    
    aPrim = sZeros;
    for (j=1:lenJ)
        for (u=1:lenU)
            for (w=1:W(u,j))
                if (((cd-cb).*V(j).*U(u,j).*w)>(co.*D(j)+ce.*L(j))) || (w==W(u,j))
                    aPrim(U(u,j),W(u,j),j) = s(U(u,j),W(u,j),j);
                end
            end
        end
    end

    G = s - aPrim;
    AStar = aPrim;
    I = 1;
    N = sum(sum(G));
    M = zeros(size(N));
    I = I .* (N + 1);
    lenI = length(I);

    a = zeros(UMax, WMax, lenJ, lenI);
    
    for (i=1:lenI)
        a(:,:,:,i) = aPrim;
        
        for (j=1:lenJ)
            if (M(:,:,j) < N(:,:,j))
                M(:,:,j) = M(:,:,j) + 1;
                for (k=1:j-1)
                    M(:,:,k) = 0;
                end
            end
        end
        
        for (j=1:lenJ)
            GMaxW = ArgMax(G(:,:,j), 1);
            GMaxU = ArgMax(G(:,:,j), 2);
            GMaxUW = GMaxW .* GMaxU;
            x = M(:,:,j);
            
            while (true)
                sMaxValues = GMaxUW .* s(:,:,j);
                sMaxTotal = sum(sum(sMaxValues));
                
                if (sMaxTotal == 0)
                    break
                end
                if (sMaxTotal >= x)
                    a(:,:,j,i) = ((GMaxUW ~= 1) .* a(:,:,j,i)) + (GMaxUW .* x);
                    break
                else
                    a(:,:,j,i) = ((GMaxUW ~= 1) .* a(:,:,j,i)) + sMaxValues;
                    x = x - sMaxTotal;

                    conditionA = s(:,:,j) <= sMaxValues;

                    WMaskz = MaskMaker(s(:,:,j), U, W, 1);
                    conditionB = WMaskz;
                    
                    WMask = MaskMaker(s(:,:,j), U, W, 0);
                    
                    GMaxUBroadcast = uwZeros;
                    GMaxWBroadcast = uwZeros;

                    for (p=1:WMax)
                        for (o=1:UMax)
                            GMaxUBroadcast(:,p) = GMaxUBroadcast(:,p) + GMaxU(:,o);
                            GMaxWBroadcast(o,:) = GMaxWBroadcast(o,:) + GMaxW(p,:);
                        end
                    end
                    conditionC = GMaxUBroadcast == 0;
                    conditionD = GMaxWBroadcast == 0;

                    result = (conditionA.*conditionB.*(conditionC + conditionD)) .* s(:,:,j);
                    GMaxW = ArgMax(result, 1);
                    GMaxU = ArgMax(result, 2);
                    
                    GMaxUW = GMaxW .* GMaxU;
                end
            end
        end
        AStar(:,:,:,i) = a(:,:,:,i);
    end
end
