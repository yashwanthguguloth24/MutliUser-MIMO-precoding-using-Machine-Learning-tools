% Date : 16-05-2022
% code written by Yashwanth guguloth

% This code finds the optimal user sets using SUS algorithm
% This code requires data provided in Data folder

% Process for all the data files
for file_num = 1:5
    disp("%%%%%%%%%%%%%%%%%%%%%%%%");
    K = 16; % number of users (based on data)
    M = 64; % number of base station antennas (based on data)
    % path
    str = '/Users/sonuyash/Documents/My docs/My work/IITH/IDP 2022/Codes/Results/Iter';
    path = append(str,num2str(file_num));
    final = append(path,'.mat');
    temp = load(final);
    h = formatInput(temp);
    opt_set = SUS(h);
    Rate = ComputeNewRate(opt_set,K,M,h);
    % Display on screen
    fprintf("File number : %d \n",file_num);
    fprintf('Optimal set = {')
    fprintf('%d ',opt_set);
    fprintf('}\n');
    fprintf('Rate is %d\n',Rate);
    disp("%%%%%%%%%%%%%%%%%%%%%%%%");
end


function So = SUS(h)
% This is the SUS algorithm function which outputs the optimal usersets.
alpha = 0.02999; % alpha can be varied
[K,M] = size(h);
iter = 1;
Tau = zeros(1,K);
for i = 1:K
    Tau(i) = i;
end

% initialization
So = [];
g_k = zeros(M,K);
for user = 1:K
    g_k(:,user) = h(user,:);
end


norms = zeros(1,K);
for user = 1:K
    norms(user) = norm(g_k(:,user));
end
[val,idx] = max(norms);
pi_i = [idx];
So = [So,pi_i];
subspace = [g_k(:,idx)];

nTau = [];
for useri = 1:length(Tau)
    if pi_i(1) ~= Tau(useri)
        innerprod = dot(h(Tau(useri),:),subspace(:,1))/(norm(h(pi_i(1),:))*norm(subspace(:,1)));
        if (abs(innerprod)) < alpha
            nTau = [nTau,Tau(useri)];
        end
    end
end
Tau = nTau;
iter = iter + 1;
while (~isempty(Tau)) && (length(So) < M)
    g_k = zeros(M,length(Tau));
    for k = 1:length(Tau)
        component = zeros(M,1);
        for j = 1:iter-1
            mag = dot(h(Tau(k),:),subspace(:,j))/((norm(subspace(:,j)))^2);  
            component = component + mag*(subspace(:,j));
        end
        g_k(:,k) = h(Tau(k),:)'-component;
    end

    norms = zeros(1,length(Tau));
    for user = 1:length(Tau)
        norms(user) = norm(g_k(:,user));
    end
    [val,idx] = max(norms);
    pi_i = [pi_i,Tau(idx)];
    So = [So,pi_i(length(pi_i))];
    subspace = [subspace,g_k(:,idx)];
    
    nTau = [];
    for useri = 1:length(Tau)
        if pi_i(iter) ~= Tau(useri)
            innerprod = dot(h(Tau(useri),:),subspace(:,iter))/(norm(h(pi_i(iter),:))*norm(subspace(:,iter)));
            if abs(innerprod) < alpha
                nTau = [nTau,Tau(useri)];
            end
        end
    end
    Tau = nTau;
    iter = iter + 1;
end


end


function h = formatInput(a)
% Formatting the data in file for clean input.
h = zeros(16,64);
for i = 1:16
    x = a.ueChannel1{1,i};
    p = reshape(x,[1,64]);
    h(i,:) = p;
end
end


function sz = calculate_size(K,M)
% calcuates number of subsets possible
% kC1 + kC2 + ..... + kCm
sz = 0;
for i = 1:min(K,M)
    sz = sz + nCr(K,i);
end
end

function x = nCr(n,r)
% computes nCr
x = factorial(n)/(factorial(r)*factorial(n-r));
end 

function rate = ComputeNewRate(S,K,M,h)
% This function computes the rate for a given h,optimal set using the formula.
set = h(S,:);
pseudo_inv = pinv(set);
sec_term = set*pseudo_inv;
[m,n] = size(sec_term);
% final_mat = eye(m) + sec_term;
rate = log2(det(abs(eye(m)+sec_term*sec_term')));
end
