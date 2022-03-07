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
bin_5 = [];
bin_6 = [];
bin_7 = [];
bin_8 = [];
j = 1;
k = 1;
l = 1;
n = 1;
o = 1;
p = 1;
q = 1;
r = 1;
for i = 1:5000
    if trainingSet(i) <= -c
        bin_1(j) = trainingSet(i);
        j = j + 1;
    elseif trainingSet(i) > -c && trainingSet(i) <= -c/2
        bin_2(k) = trainingSet(i);
        k = k + 1;
    elseif trainingSet(i) > -c/2  && trainingSet(i) <= -c/4
        bin_3(l) = trainingSet(i);
        l = l + 1;
    elseif trainingSet(i) > -c/4  && trainingSet(i) <= 0 
        bin_4(n) = trainingSet(i);
        n = n + 1;
    elseif trainingSet(i) > 0  && trainingSet(i) <= c/4
        bin_5(n) = trainingSet(i);
        o = o + 1;
    elseif trainingSet(i) > c/4  && trainingSet(i) <= c/2
        bin_6(n) = trainingSet(i);
        p = p + 1;
    elseif trainingSet(i) > c/2  && trainingSet(i) <= c
        bin_7(n) = trainingSet(i);
        q = q + 1;
    else 
        bin_8(n) = trainingSet(i);
        r = r + 1;
    end
end

sz_bin1 = size(bin_1);
sz_bin2 = size(bin_2);
sz_bin3 = size(bin_3);
sz_bin4 = size(bin_4);
sz_bin5 = size(bin_5);
sz_bin6 = size(bin_6);
sz_bin7 = size(bin_7);
sz_bin8 = size(bin_8);
c1 = (1/sz_bin1(2))*sum(bin_1);
c2 = (1/sz_bin2(2))*sum(bin_2);
c3 = (1/sz_bin3(2))*sum(bin_3);
c4 = (1/sz_bin4(2))*sum(bin_4);
c5 = (1/sz_bin5(2))*sum(bin_5);
c6 = (1/sz_bin6(2))*sum(bin_6);
c7 = (1/sz_bin7(2))*sum(bin_7);
c8 = (1/sz_bin8(2))*sum(bin_8);

initialCodeBook = [c1, c2, c3, c4, c5, c6, c7, c8];

epsilon = 0.001;
m = 1;

%Step 2 & 3:
[y1, y2, y3, y4, y5, y6, y7, y8]  =  partitionCodebook(initialCodeBook, trainingSet, m);
[y1_next, y2_next, y3_next, y4_next, y5_next, y6_next, y7_next, y8_next] = partitionCodebook(initialCodeBook, trainingSet, m+1);
%repeat step 2 & 3 until within the mean distortions being within a threshold of epsilon
while ((meanDistortion([y1, y2, y3, y4, y5, y6, y7, y8] , trainingSet) - meanDistortion([y1_next, y2_next, y3_next, y4_next, y5_next, y6_next, y7_next, y8_next], trainingSet))/(meanDistortion([y1, y2, y3, y4, y5, y6, y7, y8] , trainingSet))) >= epsilon
    m = m + 1;
    [y1, y2, y3, y4, y5, y6, y7, y8]  = partitionCodebook(initialCodeBook, trainingSet, m);
    [y1_next, y2_next, y3_next, y4_next, y5_next, y6_next, y7_next, y8_next] = partitionCodebook(initialCodeBook, trainingSet, m+1);
end

disp('The codebook yielding an approximation to the optimal code is: ')
disp([y1_next, y2_next, y3_next, y4_next, y5_next, y6_next, y7_next, y8_next])
disp('With corresponding Distoriton in bits/sourcesymbol: ')
disp(meanDistortion([y1_next, y2_next, y3_next, y4_next, y5_next, y6_next, y7_next, y8_next], trainingSet))

% fucnction that partitions the codeboopk for the m-th codeword, given
% the initial codebook (follows from step 2 of LLoyds algorithm)
function [y1, y2, y3, y4, y5, y6, y7, y8] = partitionCodebook(codeBook, tset, m)
    trainingSet = tset;
    c1 = - abs(codeBook(1)^m);
    c2 = - abs(codeBook(2)^m);
    c3 = - abs(codeBook(3)^m);
    c4 = - abs(codeBook(4)^m);
    c5 = abs(codeBook(5)^m);
    c6 = abs(codeBook(6)^m);
    c7 = abs(codeBook(7)^m);
    c8 = abs(codeBook(8)^m);
   
    %initialize bins
    bin1 = [];
    bin2 = [];
    bin3 = [];
    bin4 = [];
    bin5 = [];
    bin6 = [];
    bin7 = [];
    bin8 = [];
    p = 1;
    q = 1;
    r = 1;
    s = 1;
    t = 1;
    u = 1;
    v = 1;
    w = 1;
    % Generate partitions using the NNC and sum values in training set 
    for i = 1:5000
        if (abs(trainingSet(i) -  c1) <= abs(trainingSet(i) - c2)) && abs(trainingSet(i) - c1) <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c1) <= abs(trainingSet(i) - c4) && (abs(trainingSet(i) -  c1) <= abs(trainingSet(i) - c5)) && abs(trainingSet(i) - c1) <= abs(trainingSet(i) - c6) && abs(trainingSet(i) - c1) <= abs(trainingSet(i) - c7) && abs(trainingSet(i) - c1) <= abs(trainingSet(i) - c8)
            bin1(p) = trainingSet(i);
            p = p + 1;
        elseif abs(trainingSet(i) - c2) < abs(trainingSet(i) - c1) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c4) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c5) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c6) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c7) && abs(trainingSet(i) - c2) <= abs(trainingSet(i) - c8)
            bin2(q) = trainingSet(i);
            q = q + 1;
        elseif abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c3) < abs(trainingSet(i) - c2) && abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c4) && abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c5) && abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c6) && abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c7) && abs(trainingSet(i) - c3) <= abs(trainingSet(i) - c8)
            bin3(r) = trainingSet(i);
            r = r + 1;
        elseif abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c2) && abs(trainingSet(i) - c4) < abs(trainingSet(i) - c3) && abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c5) && abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c6) && abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c7) && abs(trainingSet(i) - c4) <= abs(trainingSet(i) - c8)
            bin4(s) = trainingSet(i);
            s = s + 1;
        elseif abs(trainingSet(i) - c5) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c5) <= abs(trainingSet(i) - c2) && abs(trainingSet(i) - c5) <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c5) < abs(trainingSet(i) - c4) && abs(trainingSet(i) - c5) <= abs(trainingSet(i) - c6) && abs(trainingSet(i) - c5) <= abs(trainingSet(i) - c7) && abs(trainingSet(i) - c5) <= abs(trainingSet(i) - c8)
            bin5(t) = trainingSet(i);
            t = t + 1;
        elseif abs(trainingSet(i) - c6) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c6) <= abs(trainingSet(i) - c2) && abs(trainingSet(i) - c6) <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c6) <= abs(trainingSet(i) - c4) && abs(trainingSet(i) - c6) < abs(trainingSet(i) - c5) && abs(trainingSet(i) - c6) <= abs(trainingSet(i) - c7) && abs(trainingSet(i) - c6) <= abs(trainingSet(i) - c8)
            bin6(u) = trainingSet(i);
            u = u + 1;
        elseif abs(trainingSet(i) - c7) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c7)  <= abs(trainingSet(i) - c2) && abs(trainingSet(i) - c7)  <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c7) <= abs(trainingSet(i) - c4) && abs(trainingSet(i) - c7)  <= abs(trainingSet(i) - c5) && abs(trainingSet(i) - c7)  < abs(trainingSet(i) - c6) && abs(trainingSet(i) - c7) <= abs(trainingSet(i) - c8)
            bin7(v) = trainingSet(i);
            v = v + 1;
        elseif abs(trainingSet(i) - c8) <= abs(trainingSet(i) - c1) && abs(trainingSet(i) - c8)  <= abs(trainingSet(i) - c2) && abs(trainingSet(i) - c8)  <= abs(trainingSet(i) - c3) && abs(trainingSet(i) - c8) <= abs(trainingSet(i) - c4) && abs(trainingSet(i) - c8)  <= abs(trainingSet(i) - c5) && abs(trainingSet(i) - c8)  <= abs(trainingSet(i) - c6) && abs(trainingSet(i) - c8)  < abs(trainingSet(i) - c7)
            bin8(w) = trainingSet(i);
            w = w + 1;
        end
    end

    sz_bin1 = size(bin1);
    sz_bin2 = size(bin2);
    sz_bin3 = size(bin3);
    sz_bin4 = size(bin4);
    sz_bin5 = size(bin5);
    sz_bin6 = size(bin6);
    sz_bin7 = size(bin7);
    sz_bin8 = size(bin8);
    %calculate the m+1-th codebook C_m

    y1 = 1/(sz_bin1(2))*sum(bin1);
    y2 = 1/(sz_bin2(2))*sum(bin2);
    y3 = 1/(sz_bin3(2))*sum(bin3);
    y4 = 1/(sz_bin4(2))*sum(bin4);
    y5 = 1/(sz_bin5(2))*sum(bin5);
    y6 = 1/(sz_bin6(2))*sum(bin6);
    y7 = 1/(sz_bin7(2))*sum(bin7);
    y8 = 1/(sz_bin8(2))*sum(bin8);

end
  
% computes the mean distortion after the m-th iteration of the code
% (follows from step 3 of LLoyds algorithm)
function D_m = meanDistortion(codeBook, tset)
    trainingSet = tset;
    c1 = codeBook(1);
    c2 = codeBook(2);
    c3 = codeBook(3);
    c4 = codeBook(4);
    c5 = codeBook(5);
    c6 = codeBook(6);
    c7 = codeBook(7);
    c8 = codeBook(8);
    sum_distortion = 0;
    for i = 1:5000
        if trainingSet(i) <= (c1+c2)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c1)^2;

        elseif trainingSet(i) > (c1+c2)/2 && trainingSet(i) <= (c2+c3)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c2)^2;

        elseif trainingSet(i) > (c2+c3)/2 && trainingSet(i) <= (c3+c4)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c3)^2;

        elseif trainingSet(i) > (c3+c4)/2 && trainingSet(i) <= (c4+c5)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c4)^2;

        elseif trainingSet(i) > (c4+c5)/2 && trainingSet(i) <= (c5+c6)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c5)^2;

        elseif trainingSet(i) > (c5+c6)/2 && trainingSet(i) <= (c6+c7)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c6)^2;

        elseif trainingSet(i) > (c6+c7)/2 && trainingSet(i) <= (c7+c8)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c7)^2;

        elseif trainingSet(i) > (c7+c8)/2
            sum_distortion = sum_distortion + (trainingSet(i) - c8)^2;

        end
    end
    D_m = (1/5000)*sum_distortion;
end



