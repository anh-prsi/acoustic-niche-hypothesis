function [] = null_new(f1, freq_range, s3)

%null model, without mass constraint and no noise/rounding

close all;

load freq_mass_all.mat;

f = f_all;

binrange = 10:30;

Einv = NaN(1,1000);
S3 = NaN(1,1000);

f(isnan(f)) = [];

f(f>freq_range(2)) = [];
f(f<freq_range(1)) = [];

f1(isnan(f1)) = [];

f1(f1>freq_range(2)) = [];
f1(f1<freq_range(1)) = [];

samplefraction = 1;

     %subsampling
   %  ix = randperm(length(f1));
   %  ix = ix(1:ceil(samplefraction*length(f1)));
   %  f1 = f1(ix);
%     m1 = m1(ix);
   %  [f1,ix] =  sort(f1);
 %    m1= m1(ix);

nsample = length(f1);

ftemp = f1;
ftemp = (ftemp-min(ftemp))/(max(ftemp) - min(ftemp));
p = diff(sort(ftemp));
D = sum(p.^2);
S = nsample;
Einv0 = (1./D)/(S-1);

switch s3
    case 'number_av'
        temp =[];
        for bins = binrange
            n1 = hist(ftemp, bins);
            p1 = n1/(S-1);
            temp=[temp, 2/(1 +bins*sum(p1.^2))];
        end
        S30 = mean(temp);
    case 'spatial_av'
        bins = 10;
        delta = 10;
        temp =[];
        for q = 0:1/(delta*bins):1/bins
            f0 = mod(ftemp-q,1);
            n1 = hist(f0, bins);
            p1 = n1/(S-1);
            temp=[temp, 2/(1+bins*sum(p1.^2))];
        end
        S30 = mean(temp);
end

fsample = NaN(nsample, 1000);

for j = 1:10000
        ixrand = randperm(length(f));
        ixrand = ixrand(1:nsample);
        ftemp = sort(f(ixrand));
        
        ftemp = (ftemp-min(ftemp))/(max(ftemp) - min(ftemp));
        p = diff(sort(ftemp));
        D = sum(p.^2);
        S = nsample;
        Einv(j) = (1./D)/(S-1);
        
    switch s3
        case 'number_av'
            temp =[];
            for bins = binrange
                n1 = hist(ftemp, bins);
                p1 = n1/(S-1);
                temp=[temp, 2/(1 +bins*sum(p1.^2))];
            end
            S3(j) = mean(temp);
        case 'spatial_av'
            temp =[];
            for q = 0:1/(delta*bins):1/bins
                f0 = mod(ftemp-q,1);
                n1 = hist(f0, bins);
                p1 = n1/(S-1);
                temp=[temp, 2/(1+bins*sum(p1.^2))];
            end
            S3(j) = mean(temp);
    end
         
end

Einv_score = length(Einv(Einv < Einv0))/ length(Einv);
S3_score = length(S3(S3 < S30))/ length(S3);

[mean(S3_score), mean(Einv_score)]

% ixmax=find(S3 == max(S3),1);
% hold on;scatter(1:length(fsample(:,ixmax)), fsample(:,ixmax), [], 'r.')
% scatter(1:length(fsample(:,ixmax)), fsample(:,ixmax), [], ixsample(:, ixmax))
% scatter(1:length(f1), f1, 'gs')
