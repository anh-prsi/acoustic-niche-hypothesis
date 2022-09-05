function [h0, h, Sim30, Sim3] = similarity_subLA_swap(tot_it, loc1, loc2)

load('sub_LA.mat')

binrange = 3:10;

ftemp1 = f(find(h(:,loc1) == 1));
ftemp2 = f(find(h(:,loc2) == 1));

%ix1 = find(h(:,loc1) == 1);
%ix1 = ix1(randperm(length(ix1))); 
%ix1 = ix1(1:length(ix1)/2); 
%ix2 = find(h(:,loc2) == 1);
%ix2 = ix2(randperm(length(ix2))); 
%ix2 = ix2(1:length(ix2)/2); 

%ftemp2 = f([ix1; ix2]); 

f1=ftemp1;
f2=ftemp2;
S1 = length(f1);
S2 = length(f2);

temp =[];
temp4=[];
for bins = binrange
    [n0, b] = hist([f1, f2], bins);
    n1 = hist(f1, b);
    p1 = n1/S1;
    n2 = hist(f2, b);
    p2 = n2/S2;
    temp=[temp, 2*sum(p1.*p2)/(sum(p1.^2) + sum(p2.^2))];
    temp4a = ((p1+p2).*log(p1+p2) - p1.*log(p1) - p2.*log(p2))/(2*log(2));
    temp4a(isnan(temp4a)) = 0;
    temp4 = [temp4, sum(temp4a)];
end
Sim30 = mean(temp);
Sim40 = mean(temp4);


Sim3 = NaN(1, tot_it);
Sim4 = NaN(1, tot_it);

it=1;

h0 = h;

while (it <= tot_it)
    h=h0;
    nswaps = 1;
    
    while (nswaps < 50000)
        hclass= h;
        if (size(hclass,1) >1)
            s_i = randi(size(hclass,1));
            occ_i = hclass(s_i,:);
            h_i = find(occ_i==1);
            if (isempty(h_i)) 
                break;
            end
            h_i = h_i(randi(length(h_i)));
            s_j = find(hclass(:,h_i) ==0);
            if (~isempty(s_j))
                s_j = s_j(randi(length(s_j)));
                h_j = find(hclass(s_j,:) ==1);
                if (isempty(h_j)) 
                    break;
                end
                h_j = h_j(randi(length(h_j)));
                if (hclass(s_i,h_j)==0)
                    hclass(s_i, h_i) = 0;
                    hclass(s_i, h_j) = 1;
                    hclass(s_j, h_j) = 0;
                    hclass(s_j, h_i) = 1;
                    nswaps = nswaps+1;
                    h = hclass;
                end
            end
        end
    end
    
    ftemp1 = f(find(h(:,loc1) == 1));
    ftemp2 = f(find(h(:,loc2) == 1));
    
    f1 = ftemp1;
    S1 = length(f1);
    f2 = ftemp2;
    S2 = length(f2);
    
    temp =[];
    temp4 = [];
    for bins = binrange
        [n0, b] = hist([f1, f2], bins);
        n1 = hist(f1, b);
        p1 = n1/S1;
        n2 = hist(f2, b);
        p2 = n2/S2;
        temp=[temp, 2*sum(p1.*p2)/(sum(p1.^2) + sum(p2.^2))];
        temp4a = ((p1+p2).*log(p1+p2) - p1.*log(p1) - p2.*log(p2))/(2*log(2));
        temp4a(isnan(temp4a)) = 0;
        temp4 = [temp4, sum(temp4a)];
    end
    Sim3(it) = mean(temp);
    Sim4(it) = mean(temp4);
    
    it = it+1;
end

[length(Sim3(Sim3<Sim30))/length(Sim3), length(Sim4(Sim4<Sim40))/length(Sim4)]