% Date : 16-05-2022
% code written by Yashwanth guguloth

% This code finds the optimal user sets using bruteforce search
% This code requires data provided in Data folder

% processing for all the files
for file_num = 1:5
    disp("%%%%%%%%%%%%%%%%%%%%%%%%");
    K = 16;
    M = 64;
    str = '/Users/sonuyash/Documents/My docs/My work/IITH/IDP 2022/Codes/Results/Iter';
    % path of data file
    path = append(str,num2str(file_num));
    final = append(path,'.mat');
    temp = load(final);
    h = formatInput(temp);
    [user_comb,indices,rates,lengths,opt_set] = orthogonal_bruteforce(h,K,M);
    Rate = ComputeNewRate(opt_set,K,M,h);
    fprintf("File number : %d \n",file_num);
    fprintf('Optimal set = {')
    fprintf('%d ',opt_set);
    fprintf('}\n')
    fprintf('Rate = %d \n',Rate);
    % Displaying all the orthogonal sets for reference
    fprintf("Orthogonal-sets  |   Cardinality |    rate\n")
    for j = 1:length(indices)
        fprintf('{')
        fprintf('%d ',filter_set(user_comb(:,indices(j))));
        fprintf('}  |  %d |  %d\n',lengths(j),rates(j))
    end
    disp("%%%%%%%%%%%%%%%%%%%%%%%%");
end

fprintf('\n');
disp('*************************************')
fprintf('\n');

function [user_comb,indices,rates,lengths,optimal_set] = orthogonal_bruteforce(h,K,M)
% This function outputs optimal sets based on bruteforce search.
sz = calculate_size(K,M); % size
is_orthogonal = zeros(1,sz);
user_comb = zeros(M,sz);
cntr = 1;
for set_size = 1:min(M,K)
    combs = nchoosek(1:K,set_size);
    [rows,cols] = size(combs);
    for comb = 1:rows
        mat = h(combs(comb,:),:);
        isorthogonal = CheckOrthogonality(mat);
        is_orthogonal(1,cntr) = isorthogonal;
        user_comb(1:length(combs(comb,:)),cntr) = combs(comb,:);
        cntr = cntr+1;
    end
end

L = sum(is_orthogonal);
rates = zeros(1,L);
indices = zeros(1,L);
lengths = zeros(1,L);

cnt = 1;
for i = 1:sz
    if is_orthogonal(1,i) == 1
        user_set = user_comb(:,i);
        fuser_set = filter_set(user_set);
        rates(1,cnt) = ComputeNewRate(fuser_set,K,M,h);
        indices(1,cnt) = i;
        lengths(1,cnt) = length(fuser_set);
        cnt = cnt+1;
    end
end

[max_val,max_idx] = max(rates);
optimal_set = filter_set(user_comb(:,indices(max_idx)));

end


function isorthogonal = CheckOrthogonality(h)
% This function checks whether H(S)*H(S).T ~= I or not
% We allow some tolerance epsilon because H(S)*H(S).T will not be equal to
% I due to computational constraints.

epsilon = 0.05;
[R,C] = size(h);
for row = 1:R
    h(row,:) = h(row,:)/(norm(h(row,:)));
end
res = h*h';
% 
flag = 1;
for i = 1:R
    for j = 1:R
        if i == j
            if (abs(1-res(i,j))) > epsilon
                flag = 0;
            end
        else 
            if (abs(res(i,j))) > epsilon
                flag = 0;
            end
        end
    end
end

if flag == 0
    isorthogonal = 0;
else
    isorthogonal = 1;
end
end




function rate = ComputeNewRate(S,K,M,h)
% Computes rate of a  given set
set = h(S,:);
pseudo_inv = pinv(set);
sec_term = set*pseudo_inv;
[m,n] = size(sec_term);
% final_mat = eye(m) + sec_term;
rate = log2(det(abs(eye(m)+sec_term*sec_term')));
end




function arr = filter_set(A)
% extracts non zero values from a set.
arr = [];
for i = 1:length(A)
    if ~isequal(A(i),0)
        arr = [arr A(i)];
    end
end
end


function h = formatInput(a)

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