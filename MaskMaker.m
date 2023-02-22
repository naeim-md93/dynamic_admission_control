function [WMask, UMask] = MaskMaker(sTemp, U, W, z)
    
    lenU = size(U, 1);
    UMax = max(max(U));
    lenW = size(W, 1);
    WMax = max(max(W));
    lenJ = size(sTemp, 3);
    lentau = size(sTemp, 4);

    UMask = zeros(size(sTemp));
    WMask = zeros(size(sTemp));
    
    for (tau=1:lentau)
        for (j=1:lenJ)
            for (u=1:lenU)
                for (w=1:lenW)
                    UMask(U(u,j),:,j,tau) = 1;
                    WMask(U(u,j),1:W(u,j)-z,j,tau) = 1;
                end
            end
        end
    end
end
