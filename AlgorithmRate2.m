% Jacob Cohen (20069127)


% generate 1 by 5000 array of random samples from the standard 
% normal distirbution to simulate a source 
rng('default')
trainingSet = randn(1,5000);
% Apply Lloyd algorithm for training data and MSE:

%Step 1: (inputs)
epsilon = 0.001;
m = 1;
% MSE with quantizer rate R = 2, recalling R = log_{2}(N), then  N = 4
% given the training set we approximate the distribution of X by
% a vector of the probabilities
sourceDistApprx = [];
for i = 1:5000
    sourceDistApprx(i) = (1/i);
end
% For MSE optimal, initialize codebook using Lloyd-Max conidtions
% first construct 4-level quantizer with c being the expectation of X
% given X > 0 where X ~ sourceDistApprx, since it is discreet we use
% discrete expectation to approximate the continous expectation:
c = 0;
for i = 1:5000 
   c = c+ 2*trainingSet(i)*sourceDistApprx(i);
end
% thus the inital codebook 1 is
codeBookRate1_1 = [-c, -c/2 ,c/2, c/2];
% with corresponding bins, R_1 = (-inf, (-3/2)c] and R_2 = ((-3/2)c, 0]
% R_3 = (0, (3/2)c] and r_4 = ((3/2)c, inf)
% with the 2-level quantizer Q(x) = c if x >  0 and Q(x) = -c if x<=0
% with distortion

% step 2: 
[y1, y2, y3, y4] = partitionCodebook(c,1);
% Step3:
while ((meanDistortion(c,m) - meanDistortion(c,m+1))/(meanDistortion(c,m))) < epsilon
    [y1, y2, y3, y4] = partitionCodebook(c,m);
    m = m + 1;
end
disp('The codebook yielding an approximation to the optimal code is: ')
disp([y1, y2, y3, y4])



% fucnction that partitions the codeboopk for the m-th codeword, given
% the initial codebook [-c, c] (follows from step 2 of LLoyds algorithm)
function [y1_m, y2_m, y3_m, y4_m]  = partitionCodebook(c,m)
    trainingSet = randn(1,5000);
    %initialize bins
    bin1 = [];
    bin2 = [];
    bin3 = [];
    bin4 = [];
    j = 1;
    s = 1;
    t = 1;
    l = 1;
    % Generate partitions using the NNC
    for i = 1:5000
        if (abs(trainingSet(i) - c^m) <= abs(trainingSet(i) + c^m)) && (abs(trainingSet(i) - c^m) <= abs(trainingSet(i) + (c/2)^m)) && (abs(trainingSet(i) - c^m) <= abs(trainingSet(i) - (c/2)^m))
            bin1(j) = trainingSet(i);
            j = j + 1;
        elseif (abs(trainingSet(i) - (c/2)^m) <= abs(trainingSet(i) + c^m)) && (abs(trainingSet(i) - (c/2)^m) <= abs(trainingSet(i) + (c/2)^m)) && (abs(trainingSet(i) - (c/2)^m) <= abs(trainingSet(i) - (c)^m))
            bin2(s) = trainingSet(i);
            s = s + 1;
        elseif (abs(trainingSet(i) + (c/2)^m) <= abs(trainingSet(i) + c^m)) && (abs(trainingSet(i) + (c/2)^m) <= abs(trainingSet(i) - (c/2)^m)) && (abs(trainingSet(i) + (c/2)^m) <= abs(trainingSet(i) - (c)^m))
            bin3(t) = trainingSet(i);
            t = t + 1;     
        elseif (abs(trainingSet(i) + (c)^m) <= abs(trainingSet(i) + (c/2)^m)) && (abs(trainingSet(i) + (c)^m) <= abs(trainingSet(i) - (c/2)^m)) && (abs(trainingSet(i) + (c)^m) <= abs(trainingSet(i) - (c)^m))
            bin4(l) = trainingSet(i);
            l = l + 1;
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
    y2_m = (1/(sz_bin2(2)))*(sum_bin2);

    sum_bin3 = 0;
    sz_bin3 = size(bin2);
    for i = 1:sz_bin2(2)    
        sum_bin3 = sum_bin3 + bin3(i);
    end
    y3_m = (1/(sz_bin3(2)))*(sum_bin3);

    sum_bin4 = 0;
    sz_bin4 = size(bin4);
    for i = 1:sz_bin4(2)    
        sum_bin4 = sum_bin4 + bin4(i);
    end
    y4_m = (1/(sz_bin4(2)))*(sum_bin4);

end
  
% computes the mean distortion after the m-th iteration of the code
% (follows from step 3 of LLoyds algorithm)
function D_m = meanDistortion(c, m)
    trainingSet = randn(1,5000);
    sum_distortion = 0;
    for i = 1:5000
        if trainingSet(i) > c
            sum_distortion = sum_distortion + (trainingSet(i) - c^m);
        elseif trainingSet(i) <= c && trainingSet(i) > c/2
            sum_distortion = sum_distortion + (trainingSet(i) - (c/2)^m);
        elseif trainingSet(i) >= -c && trainingSet(i) <= -c/2
            sum_distortion = sum_distortion + (trainingSet(i) + (c/2)^m);
        elseif trainingSet(i) <= -c 
            sum_distortion = sum_distortion + (trainingSet(i) + (c)^m);
        end
    end
    D_m = (1/5000)*sum_distortion;
end








