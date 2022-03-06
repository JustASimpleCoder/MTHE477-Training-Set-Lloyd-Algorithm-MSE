% Jacob Cohen (20069127)
% Apply Lloyd algorithm for training data and MSE:
%Step 1: (inputs)
% generate 1 by 5000 array of random samples from the standard normal
% distirbution to simulate a training set with a gaussian distribution
rng('shuffle')
trainingSet = randn(1,5000);

% Quantizer with rate R = 1, recalling R = log_{2}(N), then  N = 2, thus we want
% a 2-level MSE quantizer. In class we have seen, the 2-level quantizer 
% unique optimal for a gausian source is Q(x) = c if x >  0 and Q(x) = -c if x<=0
% with bins  R_1 = (-inf, 0] and R_2 = (0,inf).
% Since X~1/5000,  then the optimal codeBook symbols can be approximated by, 

for i = 1:5000 
   c = c + 2*trainingSet(i)*(1/5000);
end

% thus the inital codebook for m =1 1 is
initialCodeBook = [-c,c];

%initialize algorithm parameters
epsilon = 0.001;
m = 1;

% step 2 and step 3:
[c1_m, c2_m] =  partitionCodebook(initialCodeBook, trainingSet, m);
while ( ( meanDistortion([c1_m, c2_m], trainingSet, m) - meanDistortion([c1_m, c2_m], trainingSet, m+1) )/( meanDistortion([c1_m, c2_m],trainingSet, m)) ) < epsilon
    m = m + 1;
    [c1_m, c2_m] = partitionCodebook([c1_m, c2_m], trainingSet, m);
end
disp('The codebook yielding an approximation to the optimal code is:')
disp(initialCodeBook)

% we can find how close there are to the optimal codebook, which in class
% was shown to be 
c_optimal = sqrt(2/pi);

%hence we can compare the percentatges to the optimal
percentageToOptimal_y1 = 100 - 100*abs(c_optimal + y1_m)/c_optimal;
percentageToOptimal_y2 = 100 - 100*abs(c_optimal - y2_m)/c_optimal;

fprintf('y1_m is %.2f  perecnt optimal \n', percentageToOptimal_y1)
fprintf('y2_m is is %.2f  percent  optimal \n', percentageToOptimal_y2)



% fucnction that takes in codebook for the m-th codeword, the simulated training set
% and the m-th iteration of the Lloyd algorithm, and partitions the training set into N=2 bins using the Nearest Neighbour Condition 
% it returns the empiracal average inside each bin as the codeword for the (m+1)-th iteration [y1_(m+1), y2_(m+1)]
% (follows from step 2 of LLoyds algorithm)
function [y1, y2] = partitionCodebook(codeBook,tSet, m)
    trainingSet = tSet;
    c1 = codeBook(1);
    c2 = codeBook(2);
    %initialize bins
    bin1 = [];
    bin2 = [];
    j = 1;
    s = 1;
    % Generate partitions using the NNC
    for i = 1:5000
        if abs(trainingSet(i) - c1^m) <= abs(trainingSet(i) - c2^m)
            bin1(j) = trainingSet(i);
            j = j + 1;
        else %assign lower index with equality
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
    y1 = (1/(sz_bin1(2)))*(sum_bin1);
    
    sum_bin2 = 0;
    sz_bin2 = size(bin2);
    for i = 1:sz_bin2(2)    
        sum_bin2 = sum_bin2 + bin2(i);
    end
    y2 = (1/(sz_bin2(2)))*(sum_bin2);
end
  
% computes the mean distortion after the m-th iteration of the code
% (follows from step 3 of LLoyds algorithm)
function D_m = meanDistortion(codeBook,tSet, m)
    trainingSet = tSet;
    c1 = codeBook(1);
    c2 = codeBook(2);
    
    sum_distortion = 0;
    for i = 1:5000
        if trainingSet(i) > 0
            sum_distortion = sum_distortion + (trainingSet(i) - c1^m)^2;
        else
            sum_distortion = sum_distortion + (trainingSet(i) - c2^m)^2;
        end
    end
    D_m = (1/5000)*sum_distortion;
end








