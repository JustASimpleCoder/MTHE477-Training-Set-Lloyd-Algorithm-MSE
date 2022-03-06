% Jacob Cohen (20069127)


% generate 1 by 5000 array of random samples from the standard 
% normal distirbution to simulate a source 
rng('default')
trainingSet = randn(1,5000);
% Apply Lloyd algorithm for training data and MSE:

%Step 1: (inputs)
epsilon = 0.001;
m = 1;
% MSE with quantizer rate R = 1, recalling R = log_{2}(N), then  N = 2
% given the training set we approximate the distribution of X by
% a vector of the probabilities
sourceDistApprx = [];
for i = 1:5000
    sourceDistApprx(i) = (1/i);
end
% For MSE optimal, initialize codebook using Lloyd-Max conidtions
% first construct 2-level quantizer with c being the expectation of X
% given X > 0 where X ~ sourceDistApprx, since it is discreet we use
% discrete expectation to approximate the continous expectation:
c = 0;
for i = 1:5000 
   c = c+ 2*trainingSet(i)*sourceDistApprx(i);
end
% thus the inital codebook 1 is
codeBookRate1_1 = [-c ,c];
% with corresponding bins, R_1 = (-inf, 0] and R_2 = (0,inf)
% with the 2-level quantizer Q(x) = c if x >  0 and Q(x) = -c if x<=0
% with distortion

% step 2: 
[y1_m, y2_m] = partitionCodebook(c,m);

% Step3:
while ((meanDistortion(c,m)- meanDistortion(c,m+1))/(meanDistortion(c,m))) < epsilon
    [y1_m, y2_m] = partitionCodebook(c,m);
    m = m + 1;
end
disp('The codebook yielding an approximation to the optimal code is:')
disp([y1_m, y2_m])

% we can find how close there are to the optimal codebook, which in class
% was shown to be 
c_optimal = sqrt(2/pi);

%hence we can compare the percentatges to the optimal
percentageToOptimal_y1 = 100 - 100*abs(c_optimal - y1_m)/c_optimal;
percentageToOptimal_y2 = 100 - 100*abs(c_optimal - y1_m)/c_optimal;

fprintf('y1_m is %.2f  perecnt optimal \n', percentageToOptimal_y1)
fprintf('y2_m is is %.2f  percent  optimal \n', percentageToOptimal_y2)
fprintf(" Hence the code  is %.2f percent optimal \n ", max(percentageToOptimal_y1,percentageToOptimal_y2) )



% fucnction that partitions the codeboopk for the m-th codeword, given
% the initial codebook [-c, c] (follows from step 2 of LLoyds algorithm)
function [y1_m, y2_m] = partitionCodebook(c,m)
    trainingSet = randn(1,5000);
    %initialize bins
    bin1 = [];
    bin2 = [];
    j = 1;
    s = 1;
    % Generate partitions using the NNC
    for i = 1:5000
        if abs(trainingSet(i) - c^m) <= abs(trainingSet(i) + c^m)
            bin1(j) = trainingSet(i);
            j = j + 1;
        elseif abs(trainingSet(i) - c^m) > abs(trainingSet(i) + c^m)
            bin2(s) = trainingSet(i);
            s = s + 1;
        end
    end

    %calculate the m+1-th codebook C_m
    sz_bin1 = size(bin1);
    sum_bin1 = 0;
    for i = 1:sz_bin1(2)    
        sum_bin1 = sum_bin1 + bin1(i);
    end
    y1_m = (1/(sz_bin1(2)))*(sum_bin1);
    
    sum_bin2 = 0;
    sz_bin2 = size(bin2);
    for i = 1:sz_bin2(2)    
        sum_bin2 = sum_bin2 + bin2(i);
    end
    y2_m = (1/(sz_bin1(2)))*(sum_bin2);
end
  
% computes the mean distortion after the m-th iteration of the code
% (follows from step 3 of LLoyds algorithm)
function D_m = meanDistortion(c, m)
    trainingSet = randn(1,5000);
    sum_distortion = 0;
    for i = 1:5000
        if trainingSet(i) > 0
            sum_distortion = sum_distortion + (trainingSet(i) - c^m);
        else
            sum_distortion = sum_distortion + (trainingSet(i) + c^m );
        end
    end
    D_m = (1/5000)*sum_distortion;
end








