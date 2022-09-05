function [Einv_score, Eprime_score, S3_score] = fwd_test(sigmarange, samplingfraction)

binrange = 10:30;

Einv_score = NaN(100, length(sigmarange));
Eprime_score = NaN(100, length(sigmarange));
S3_score = NaN(100, length(sigmarange));

Einv = NaN(1, 100);
Eprime = NaN(1, 100);
S3 = NaN(1, 100);

for i = 1:length(sigmarange)
    for k = 1:1000
        f  = 10:210;
        f = f + sigmarange(i)*rand(size(f));
        f = f(30:170);
        
        ix = randperm(length(f));
        ix = ix(1:ceil(samplingfraction*length(f)));
        f = f(ix);
        f = (f-min(f))/(max(f)-min(f));
        p = diff(sort(f));
        D = sum(p.^2);
        S = length(f);
        pair = nchoosek(p, 2);
        
        Einv0 = (1./D)/(S-1);
        Eprime0 = 1-sum(abs(pair(:,1) - pair(:,2))/(S-1));
        
        temp =[];
        for bins = binrange
            n1 = hist(f, bins);
            p1 = n1/(S-1);
            temp=[temp, 2/(1+bins*sum(p1.^2))];
        end
        S30 = mean(temp);        
        
        for j = 1:100
            
            f = rand(1,141);
            ix = randperm(length(f));
            ix = ix(1:ceil(samplingfraction*length(f)));
            f = f(ix);
            f = (f-min(f))/(max(f)-min(f));
            p = diff(sort(f));
            D = sum(p.^2);
            S = length(f);
            pair = nchoosek(p, 2);
            
            Einv(j) = (1./D)/(S-1);
            Eprime(j) = 1-sum(abs(pair(:,1) - pair(:,2))/(S-1));
            
            temp =[];
            for bins = binrange
                n1 = hist(f, bins);
                p1 = n1/(S-1);
                temp=[temp, 2/(1+bins*sum(p1.^2))];
            end
            S3(j) = mean(temp); 
            
        end
        
        Einv_score(k, i) = length(Einv(Einv < Einv0))/ length(Einv);
        Eprime_score(k, i) = length(Eprime(Eprime < Eprime0))/ length(Eprime);
        S3_score(k, i) = length(S3(S3 < S30))/ length(S3);
        
    end
end
