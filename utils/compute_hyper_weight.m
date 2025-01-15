function ECI = compute_hyper_weight(baseClsSegs,n)
para_theta = 0.4;
%M = size(mClsLabels,2);
ETs = getAllClsEntropy(baseClsSegs);
ECI = exp(-ETs./para_theta./n);
end

function Es = getAllClsEntropy(baseClsSegs)



[~, nCls] = size(baseClsSegs);

Es = zeros(nCls,1);
for i = 1:nCls
    temp = baseClsSegs(:,i);
    idx = temp~=0;
    set = baseClsSegs(idx,:);
    set(:, all(set==0)) = [];
    Es(i) = getOneClsEntropy(set);
    clear set
end
end
function E = getOneClsEntropy(set)

E = 0;
count = sum(set); 


count = count./size(set,1);
E = E-sum(count.*log2(count));
end

