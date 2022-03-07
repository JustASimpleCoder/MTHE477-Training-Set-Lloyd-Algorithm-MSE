% Jacob Cohen (20069127)
clear all
clc


% Apply Lloyd algorithm for training data and MSE:
%Step 1: (inputs)

%training set of iid normally (gaussian) distributed source
rng('default')
trainingSet = randn(1,5000);

%estimate average value of X being postitive
c = 0;
for i = 1:5000 
   c = c + 2*trainingSet(i)*(1/5000);
end
% set up initial codebook using paritions (-inf, -c], (-c,0] (0,c], (c,inf)
% so that c1 = E[X|X <= -c] c2 = E[X|-c< X <= 0] 
% c3 = E[X|0 < X <= c] c4 = E[X| X > c ]

bin_1 = [];
bin_2 = [];
bin_3 = [];
bin_4 = [];
j = 1;
k = 1;
l = 1;
n = 1;
for i = 1:5000
    if trainingSet(i) <= -c
        bin_1(j) = trainingSet(i);
        j = j + 1;
    elseif trainingSet(i) > -c && trainingSet(i) <= 0
        bin_2(k) = trainingSet(i);
        k = k + 1;
    elseif trainingSet(i) > 0 && trainingSet(i) <= c
        bin_3(l) = trainingSet(i);
        l = l + 1;
    else 
        bin_4(n) = trainingSet(i);
        n = n + 1;
    end
end

sz_bin1 = size(bin_1);
sz_bin2 = size(bin_2);
sz_bin3 = size(bin_3);
sz_bin4 = size(bin_4);
c1 = (1/sz_bin1(2))*sum(bin_1);
c2 = (1/sz_bin2(2))*sum(bin_2);
c3 = (1/sz_bin3(2))*sum(bin_3);
c4 = (1/sz_bin4(2))*sum(bin_4);

initialCodeBook = [c1, c2, c3, c4];

epsilon = 0.001;
m = 1;

%Step 2 & 3:
[y1, y2, y3, y4]  =  partitionCodebook(initialCodeBook, trainingSet, m);
[y1_next, y2_next, y3_next, y4_next] = partitionCodebook(initialCodeBook, trainingSet, m+1);
%repeat step 2 & 3 until within the mean distortions being within a threshold of epsilon
while ((meanDistortion([y1, y2, y3, y4], trainingSet) - meanDistortion([y1_next, y2_next, y3_next, y4_next], trainingSet))/(meanDistortion([y1, y2, y3, y4], trainingSet))) >= epsilon
    m = m + 1;
    [y1, y2, y3, y4] = partitionCodebook(initialCodeBook, trainingSet, m);
    [y1_next, y2_next, y3_next, y4_next] = partitionCodebook(initialCodeBook, trainingSet, m+1);
    disp(m)
end
disp('The codebook yielding an approximation to the optimal code is: ')
disp([y1_next, y2_next, y3_next, y4_next])
disp('With corresponding Distoriton in bits/sourcesymbol: ')
disp(meanDistortion(initialCodeBook, trainingSet))

% fucnction that partitions the codeboopk for the m-th codeword, given
% the initial codebook (follows from step 2 of LLoyds algorithm)
function [y1, y2, y3, y4]  = partitionCodebook(codeBook, tset, m)
    trainingSet = tset;
    c1 = - abs(codeBook(1)^m);
    c2 = - abs(codeBook(2)^m);
    c3 = abs(codeBook(3)^m);
    c4 = abs(codeBook(4)^m);
   
    %initialize bins
    bin1 = [];
    bin2 = [];
    bin3 = [];
    bin4 = [];
    p = 1;
    q = 1;
    r = 1;
    s = 1;
    % Generate partitions using the NNC and sum values in training set 
    for i = 1:5000
        if (abs(trainingSet(i) -  c1) <= abs(trainingSet(i) - c2)) && abs(trainingSet(i) - c1) <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c1) <= abs(trainingSet(i) - c4)
            bin1(p) = trainingSet(i);
            p = p + 1;
        elseif abs(trainingSet(i) - c2) < abs(trainingSet(i) - c1) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c4)
            bin2(q) = trainingSet(i);
            q = q + 1;
        elseif abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c3) < abs(trainingSet(i) - c2) && abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c4)
            bin3(r) = trainingSet(i);
            r = r + 1;
        elseif abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c2) && abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c3)
            bin4(s) = trainingSet(i);
            s = s + 1;
        end
    end

    sz_bin1 = size(bin1);
    sz_bin2 = size(bin2);
    sz_bin3 = size(bin3);
    sz_bin4 = size(bin4);
    %calculate the m+1-th codebook C_m

    y1 = 1/(sz_bin1(2))*sum(bin1);
    y2 = 1/(sz_bin2(2))*sum(bin2);
    y3 = 1/(sz_bin3(2))*sum(bin3);
    y4 = 1/(sz_bin4(2))*sum(bin4);

end
  
% computes the mean distortion after the m-th iteration of the code
% (follows from step 3 of LLoyds algorithm)
function D_m = meanDistortion(codeBook, tset)
    trainingSet = tset;
    c1 = codeBook(1);
    c2 = codeBook(2);
    c3 = codeBook(3);
    c4 = codeBook(4);
    sum_distortion = 0;
    for i = 1:5000
        if trainingSet(i) <= (c1+c2)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c1)^2;
        elseif trainingSet(i) > (c1+c2)/2 && trainingSet(i) <= 0
            sum_distortion = sum_distortion + (trainingSet(i) - c2)^2;
        elseif trainingSet(i) > 0 && trainingSet(i) <= (c3+c4)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c3)^2;
        elseif trainingSet(i) > (c3+c4)/2 
            sum_distortion = sum_distortion + (trainingSet(i) - c4)^2;
        end
    end
    D_m = (1/5000)*sum_distortion;
end









