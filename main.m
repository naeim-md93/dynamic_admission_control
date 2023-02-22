rng(1);

% Get inputs from user:
%lentau = input('Number of States (S_tau): ');
lentau = 5;

%lenJ = input('Number of J: ');
lenJ = 3;

%U = input('U(lenU, lenJ):');
%U = sort(randi([1, 10], lenU, lenJ), 1, 'ascend');
U = [1,1,1;2,2,2;3,3,3];

%W = input('W(lenU, lenJ):');
%W = sort(randi([1, 10], lenW, lenJ), 1, 'descend');
W = [3,3,3;2,2,2;1,1,1];

%V = input('V = [1,2,...,lenJ]: ');
V = [1, 2, 3];
%V = randi([1, 10], 1, lenJ);

%R = input('R: ');
R = 15;

%N = input('N: ');
N = 5;

%costs = input('[cb;cd;co;ce]: ');
costs = [5;10;50;100];

%sigma = input('[sigmaD, sigmaL]: ');
sigma = [1.72, 2];

%mu = input('[muD, muL]:');
mu = [4, 2];

%gamma = input('gamma: ');
gamma = 0.99;

%lambda = input('lambda: ');
lambda = 0.0001;

%epsilon = input('epsilon: ');
epsilon = 0.001;

%B = input('B([1,2,...,lenJ]): ');
B = [30, 15, 30];
%B = poissrnd(ones(1,lenJ)*10);

%Rho = input('[rho1, rho2]: ');
Rho = [1, 1];
%Rho = poissrnd(ones(1,2)*2);

% Initialization:
J = 1:lenJ;
lenU = size(U, 1);
lenW = lenU;



D = lognrnd(mu(1), sigma(1), [1, lenJ]);

L = lognrnd(mu(2), sigma(2), [1, lenJ]);

% Create S by U and W
WMax = max(max(W));
UMax = max(max(U));

STemp = poissrnd(ones(UMax, WMax, lenJ, lentau)*5);
WMask = MaskMaker(STemp, U, W, 0);
S = WMask .* STemp;

Xi = sum(sum(W));

A = zeros(UMax, WMax, lenJ, lentau);
Theta = ones(Xi, lentau+1);
Z = ones(Xi, lentau+1);
P = ones(Xi, Xi, lentau+1);

for (tau=1:lentau)
    sTau = S(:,:,:,tau);
    thetaPrev = Theta(:,tau);
    zPrev = Z(:,tau);
    pPrev = P(:,:, tau);
    
    [a, z, p, theta] = RLSTDBasedADP(sTau, N, zPrev, pPrev, thetaPrev, J, V, U, W, costs, D, L, Rho, B, R, gamma, lambda, epsilon);
    
    A(:,:,:,tau) = a;
    Theta(:,tau+1) = theta;
    Z(:,tau+1) = z;
    P(:,:,tau+1) = p;
end

for (i=1:size(A,4))
    totalActionPatients = sum(sum(sum(A(:,:,:,i))));
    patients = zeros(totalActionPatients, lenJ);

    startPoint = 1;
    for (j=1:lenJ)
        totalJPatients = sum(sum(A(:,:,j,i)));
        endPoint = startPoint + totalJPatients - 1;
        patients(startPoint:endPoint, j) = 1;
        startPoint = endPoint + 1;
        writematrix(A(:,:,j,i), sprintf('a_%f.xlsx', i), 'sheet', j);
    end
    writematrix(patients, 'patients.xlsx', 'sheet', i);
end