function [Jprime, Eheip, Ecomp, Einv, Eln, F21, E, Eprime, Cav, Csd, Csd2, Cnnsd, S1, S2, S3, S4, S5, S6] = fwd_test_simple(sigmarange, bins)

samplingfraction = 1;

Jprime = NaN(length(sigmarange), 1000);
Eheip = NaN(length(sigmarange), 1000);
Ecomp = NaN(length(sigmarange), 1000);
Einv = NaN(length(sigmarange), 1000);
Eln = NaN(length(sigmarange), 1000);
F21 = NaN(length(sigmarange), 1000);
E = NaN(length(sigmarange), 1000);
Eprime = NaN(length(sigmarange), 1000);
Cav = NaN(length(sigmarange), 1000);
Csd = NaN(length(sigmarange), 1000);
Csd2 = NaN(length(sigmarange), 1000);
Cnnsd = NaN(length(sigmarange), 1000);
S1 = NaN(length(sigmarange), 1000);
S2 = NaN(length(sigmarange), 1000);
S3 = NaN(length(sigmarange), 1000);
S4 = NaN(length(sigmarange), 1000);
S5 = NaN(length(sigmarange), 1000);
S6 = NaN(length(sigmarange), 1000);

for i = 1:length(sigmarange)
    for k = 1:1000
        if (i < length(sigmarange))
            f  = 10:210;
            f = f + sigmarange(i)*randn(size(f));
            f = f(30:170);
        %  f = f(85:115);
        else
            f = rand(1,141); 
        %  f = rand(1,31);
        end
            
        ix = randperm(length(f));
         ix = ix(1:ceil(samplingfraction*length(f))); 
         f = f(ix); 
        f = (f-min(f))/(max(f)-min(f));
        p = diff(sort(f));
        H = - sum(p .*log(p));
        D = sum(p.^2); 
         S = length(f); 
        pair = nchoosek(p, 2);
        inter = nchoosek(f,2); 
        Jprime(i,k) = H/log(S-1); 
        Eheip(i,k) = (exp(H)-1)/(S-2);
        Ecomp(i,k) = (1-D)/(1-1/(S-1));
        Einv(i,k) = (1./D)/(S-1); 
        Eln(i,k) = -log(D)/log(S-1);
        F21(i,k) = (1./D - 1)./(exp(H)-1);
        E(i,k) = (sum(min(p, 1/(S-1))) - 1/(S-1))/(1-1/(S-1));
        Eprime(i,k) = 1-sum(abs(pair(:,1) - pair(:,2))/(S-1));
        Cav(i,k) = mean(abs(inter(:,1) - inter(:,2))); 
        Csd(i,k) = std(abs(inter(:,1) - inter(:,2)));
        Csd2(i,k) = std(abs(pair(:,1) - pair(:,2)));
        Cnnsd(i,k) = std(p); 
        
        n1 = hist(f, bins);
        p1=n1/(S-1);
        S1(i,k)=sum(sqrt(p1/bins));
        S2(i,k)=2*sum(p1./(1+bins*p1));            
        r =-p.*log(p)-1/bins*log(1/bins) + (p+1/bins) .*log(p+1/bins);
        r(isnan(r)) = 0;
        S4(i,k) = sum(r)/(2*log(2));
  %      S4(i,k) = 2*(S-1)/((S-bins) + bins*sum(p1.*max(n1-1, 0)));
        S5(i,k)= 1/(sqrt(bins) * sqrt(sum(p1.^2)));
        [temp, ix] = sort(n1, 'descend');
        S6(i,k) = 1-sum(-diff(p1(ix)) ./ (1:(length(p1)-1)));
        
        temp =[]; 
        for bins = 10:30
            n1 = hist(f, bins); 
            p1 = n1/(S-1);
            temp=[temp, 2/(1+bins*sum(p1.^2))];
        end
        S3(i,k) = mean(temp); 
    end
end

% errorbar(x, mean(Jprime(1:end-1,:)'), std(Jprime(1:end-1,:)'), std(Jprime(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Jprime(end,:)'), mean(Jprime(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Jprime(end,:)')+std(Jprime(end,:)'), mean(Jprime(end,:)')+std(Jprime(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Jprime(end,:)')-std(Jprime(end,:)'), mean(Jprime(end,:)')-std(Jprime(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Eheip(1:end-1,:)'), std(Eheip(1:end-1,:)'), std(Eheip(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Eheip(end,:)'), mean(Eheip(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Eheip(end,:)')+std(Eheip(end,:)'), mean(Eheip(end,:)')+std(Eheip(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Eheip(end,:)')-std(Eheip(end,:)'), mean(Eheip(end,:)')-std(Eheip(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Ecomp(1:end-1,:)'), std(Ecomp(1:end-1,:)'), std(Ecomp(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Ecomp(end,:)'), mean(Ecomp(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Ecomp(end,:)')+std(Ecomp(end,:)'), mean(Ecomp(end,:)')+std(Ecomp(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Ecomp(end,:)')-std(Ecomp(end,:)'), mean(Ecomp(end,:)')-std(Ecomp(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Einv(1:end-1,:)'), std(Einv(1:end-1,:)'), std(Einv(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Einv(end,:)'), mean(Einv(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Einv(end,:)')+std(Einv(end,:)'), mean(Einv(end,:)')+std(Einv(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Einv(end,:)')-std(Einv(end,:)'), mean(Einv(end,:)')-std(Einv(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Eln(1:end-1,:)'), std(Eln(1:end-1,:)'), std(Eln(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Eln(end,:)'), mean(Eln(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Eln(end,:)')+std(Eln(end,:)'), mean(Eln(end,:)')+std(Eln(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Eln(end,:)')-std(Eln(end,:)'), mean(Eln(end,:)')-std(Eln(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(F21(1:end-1,:)'), std(F21(1:end-1,:)'), std(F21(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(F21(end,:)'), mean(F21(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(F21(end,:)')+std(F21(end,:)'), mean(F21(end,:)')+std(F21(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(F21(end,:)')-std(F21(end,:)'), mean(F21(end,:)')-std(F21(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(E(1:end-1,:)'), std(E(1:end-1,:)'), std(E(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(E(end,:)'), mean(E(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(E(end,:)')+std(E(end,:)'), mean(E(end,:)')+std(E(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(E(end,:)')-std(E(end,:)'), mean(E(end,:)')-std(E(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Eprime(1:end-1,:)'), std(Eprime(1:end-1,:)'), std(Eprime(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Eprime(end,:)'), mean(Eprime(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Eprime(end,:)')+std(Eprime(end,:)'), mean(Eprime(end,:)')+std(Eprime(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Eprime(end,:)')-std(Eprime(end,:)'), mean(Eprime(end,:)')-std(Eprime(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Cav(1:end-1,:)'), std(Cav(1:end-1,:)'), std(Cav(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Cav(end,:)'), mean(Cav(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Cav(end,:)')+std(Cav(end,:)'), mean(Cav(end,:)')+std(Cav(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Cav(end,:)')-std(Cav(end,:)'), mean(Cav(end,:)')-std(Cav(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Csd(1:end-1,:)'), std(Csd(1:end-1,:)'), std(Csd(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Csd(end,:)'), mean(Csd(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Csd(end,:)')+std(Csd(end,:)'), mean(Csd(end,:)')+std(Csd(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Csd(end,:)')-std(Csd(end,:)'), mean(Csd(end,:)')-std(Csd(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Csd2(1:end-1,:)'), std(Csd2(1:end-1,:)'), std(Csd2(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Csd2(end,:)'), mean(Csd2(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Csd2(end,:)')+std(Csd2(end,:)'), mean(Csd2(end,:)')+std(Csd2(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Csd2(end,:)')-std(Csd2(end,:)'), mean(Csd2(end,:)')-std(Csd2(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(Cnnsd(1:end-1,:)'), std(Cnnsd(1:end-1,:)'), std(Cnnsd(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(Cnnsd(end,:)'), mean(Cnnsd(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(Cnnsd(end,:)')+std(Cnnsd(end,:)'), mean(Cnnsd(end,:)')+std(Cnnsd(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(Cnnsd(end,:)')-std(Cnnsd(end,:)'), mean(Cnnsd(end,:)')-std(Cnnsd(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(S1(1:end-1,:)'), std(S1(1:end-1,:)'), std(S1(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(S1(end,:)'), mean(S1(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(S1(end,:)')+std(S1(end,:)'), mean(S1(end,:)')+std(S1(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(S1(end,:)')-std(S1(end,:)'), mean(S1(end,:)')-std(S1(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(S2(1:end-1,:)'), std(S2(1:end-1,:)'), std(S2(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(S2(end,:)'), mean(S2(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(S2(end,:)')+std(S2(end,:)'), mean(S2(end,:)')+std(S2(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(S2(end,:)')-std(S2(end,:)'), mean(S2(end,:)')-std(S2(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(S3(1:end-1,:)'), std(S3(1:end-1,:)'), std(S3(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(S3(end,:)'), mean(S3(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(S3(end,:)')+std(S3(end,:)'), mean(S3(end,:)')+std(S3(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(S3(end,:)')-std(S3(end,:)'), mean(S3(end,:)')-std(S3(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(S4(1:end-1,:)'), std(S4(1:end-1,:)'), std(S4(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(S4(end,:)'), mean(S4(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(S4(end,:)')+std(S4(end,:)'), mean(S4(end,:)')+std(S4(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(S4(end,:)')-std(S4(end,:)'), mean(S4(end,:)')-std(S4(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(S5(1:end-1,:)'), std(S5(1:end-1,:)'), std(S5(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(S5(end,:)'), mean(S5(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(S5(end,:)')+std(S5(end,:)'), mean(S5(end,:)')+std(S5(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(S5(end,:)')-std(S5(end,:)'), mean(S5(end,:)')-std(S5(end,:)')], 'r--', 'linewidth', 1);
% errorbar(x, mean(S6(1:end-1,:)'), std(S6(1:end-1,:)'), std(S6(1:end-1,:)')); hold on; plot([min(x), max(x)], [mean(S6(end,:)'), mean(S6(end,:)')], 'r', 'linewidth', 2); plot([min(x), max(x)], [mean(S6(end,:)')+std(S6(end,:)'), mean(S6(end,:)')+std(S6(end,:)')], 'r--', 'linewidth', 1); plot([min(x), max(x)], [mean(S6(end,:)')-std(S6(end,:)'), mean(S6(end,:)')-std(S6(end,:)')], 'r--', 'linewidth', 1);