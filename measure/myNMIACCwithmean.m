function [resmax]= myNMIACCwithmean(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,size(U,2));
%maxIter = 1;

for iter = 1:1
    indx = kmeans(U_normalized,numclass,'MaxIter',100, 'Replicates',3); %3
    %indx = litekmeans(U_normalized,numclass,'MaxIter',100, 'Replicates',3);
    indx = indx(:);
    result(iter,:) = Clustering8Measure(Y,indx); % result = [ACC nmi Purity Fscore Precision Recall AR Entropy];
end
resmax = result;
%resmax = result;
%resmax = mean(result,1);
%resstd = std(result);
%resmax = mean(result,1);

end