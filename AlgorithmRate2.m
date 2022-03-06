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
% given the training set we approximate the distribution of X 
for i = 1:5000 
   c = c + 2*trainingSet(i)*(1/5000);
end

% thus the inital codebook 1 is
initialCodeBook = [-c, -c/2 ,c/2, c];
% with corresponding bins, R_1 = (-inf, (-3/2)c] and R_2 = ((-3/2)c, 0]
% R_3 = (0, (3/2)c] and r_4 = ((3/2)c, inf)
% with the 2-level quantizer Q(x) = c if x >  0 and Q(x) = -c if x<=0
% with distortion


[c1_m, c2_m, c3_m, c4_m] =  partitionCodebook(initialCodeBook, trainingSet, m);
while ((meanDistortion([c1_m, c2_m, c3_m, c4_m], trainingSet,m) - meanDistortion([c1_m, c2_m, c3_m, c4_m], trainingSet, m+1))/(meanDistortion([c1_m, c2_m, c3_m, c4_m], trainingSet, m))) < epsilon
    [c1_m, c2_m, c3_m, c4_m] = partitionCodebook([-c, -c/2 ,c/2, c/2],m);
    m = m + 1;
end
disp('The codebook yielding an approximation to the optimal code is: ')
disp([c1_m, c2_m, c3_m, c4_m])
disp(meanDistortion([c1_m, c2_m, c3_m, c4_m], trainingSet,m))



% fucnction that partitions the codeboopk for the m-th codeword, given
% the initial codebook [-c, c] (follows from step 2 of LLoyds algorithm)
function [y1_m, y2_m, y3_m, y4_m]  = partitionCodebook(codeBook, tset, m)
    trainingSet = tset;
    c1 = codeBook(1);
    c2 = codeBook(2);
    c3 = codeBook(3);
    c4 = codeBook(4);
    
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
        if (abs(trainingSet(i) - c1^m) <= abs(trainingSet(i) - c2^m)) && abs(trainingSet(i) - c1^m) <= abs(trainingSet(i) - c3^m) && abs(trainingSet(i) - c1^m) <= abs(trainingSet(i) - c4^m)
            bin1(j) = trainingSet(i);
            j = j + 1;
        elseif abs(trainingSet(i) - c2^m) <= abs(trainingSet(i) - c1^m) && abs(trainingSet(i) - c2^m) <= abs(trainingSet(i) - c3^m) && abs(trainingSet(i) - c2^m) <= abs(trainingSet(i) - c4^m)
            bin2(s) = trainingSet(i);
            s = s + 1;
        elseif abs(trainingSet(i) - c3^m) <= abs(trainingSet(i) - c1^m) && abs(trainingSet(i) - c3^m) <= abs(trainingSet(i) - c2^m) && abs(trainingSet(i) - c3^m) <= abs(trainingSet(i) - c4^m)
            bin3(t) = trainingSet(i);
            t = t + 1;     
        elseif abs(trainingSet(i) - c4^m) <= abs(trainingSet(i) - c1^m) && abs(trainingSet(i) - c4^m) <= abs(trainingSet(i) - c2^m) && abs(trainingSet(i) - c4^m) <= abs(trainingSet(i) - c3^m)
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
function D_m = meanDistortion(codeBook, tset, m)
    trainingSet = tset;
    c1 = codeBook(1);
    c2 = codeBook(2);
    c3 = codeBook(3);
    c4 = codeBook(4);


    sum_distortion = 0;
    for i = 1:5000
        if trainingSet(i) > c4
            sum_distortion = sum_distortion + (trainingSet(i) - c1^m)^2;
        elseif trainingSet(i) <= c4 && trainingSet(i) > c3
            sum_distortion = sum_distortion + (trainingSet(i) - (c2)^m)^2;
        elseif trainingSet(i) <= c3 && trainingSet(i) <= c2
            sum_distortion = sum_distortion + (trainingSet(i) - (c3)^m)^2;
        elseif trainingSet(i) <= c1 
            sum_distortion = sum_distortion + (trainingSet(i) - (c4)^m)^2;
        end
    end
    D_m = (1/5000)*sum_distortion;
end








