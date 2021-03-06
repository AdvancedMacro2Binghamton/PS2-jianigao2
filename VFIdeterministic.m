close all
<<<<<<< HEAD
delta=0.025              %depreciation rate
alpha=0.35               %capital share of income
beta=0.99                %discount rate
sigma=2                  %utility parameter


nbk=20                 %number of data points int he grid
nba=2                    %number of value for the shocks
crit=1                   %convergence criterion
epsi=1e-1                %convergence parameter


PI=[phh 1-phh;1-pll,pll] %transition matrix
pihh=0.977               %this is pi hh
pill=0.926               %p ll
ah=1.1
al=0.678
A=[ah al]                %stochastic matrix


kmin=0.1
kmax=5
k=linspace(kmin,kmax,nbk)% k grid
c=zeros(nbk,nba)
util=zeros(nbk,nba)
v=zeros(nbk,nba)
Tv=zeros(nbk,nba)


%%%%%%%%%%%%%%%%%%set up%%%%%%%%%%%%%%%%%%%%%

while crit>epsi
    for i=1:nbk
        for j=1:nba
            c=A(j)*k(i)^alpha+(1-delta)*k(i)-k
            neg=find(c<0)
            c(neg)=NaN
            util(:,j)=(c.^(1-sigma))/(1-sigma)
            util(neg,j)=-1e12
        end
        [Tv(i,:),pl(1,:)]=max(util+beta*(v*PI))
    end
    crit=max(abs(Tv-v))
    v=Tv
    i=i+1
=======
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons = k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 

ret = cons .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret(cons < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(1, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat = ret + beta * repmat(v_guess, [num_k 1]);
    
    % find the optimal k' for every k:
    [vfn, pol_indx] = max(value_mat, [], 2);
    vfn = vfn';
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
>>>>>>> parent of 902b1ba... trying to add another loop to calibriate A
end

plot(v(:,1))
hold on 
plot(v(:,2))


