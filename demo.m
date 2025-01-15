clear
data = 'tr12_base_clustering.mat';
load(data);
n = 20;
k = length(unique(gt));
%Two main parameters k and alpha
anchor_all = [k-1:1:k+3];  %[k-1:k+3]
alpha_all = [0.9:0.01:1]; 
for anchor_ind = 1:length(anchor_all)   
    anchor = anchor_all(anchor_ind);  % parameter k: the number of anchors
for alpha_ind = 1:length(alpha_all)
    alpha = alpha_all(alpha_ind);  % parameter alpha: similarity threshold
    tic
parfor i= 1:10 
    % 10 ensemble experiments were 
    % conducted using different combinations of base clustering results.
    [a,b] = size(members);
    zz = RandStream('mt19937ar','Seed',i);
    RandStream.setGlobalStream(zz);
    indx = randperm(b);
    EC_end = members(:,indx(1:n));% Base clustering result
    H = Gbe(EC_end);  % Initial hypergraph
    [~, mClsLabels] = computeMicroclusters(H);
    Y = generate_local_consensus(mClsLabels,k);  %Generate local consensus Matrix Y

    simOfCluster = full(simxjac(H'));
    RWofCluster = RandomWalkofCluster(simOfCluster); % Hyperedge random walks
    H_enhance = hyper_enhance(H,RWofCluster,alpha);  %Hypergraph enhancement
    edge_weight = compute_hyper_weight(H_enhance,n);  % Hyperedge weights
    [U,F] = solve(H_enhance,anchor,edge_weight',Y); %optimization
    res = myNMIACCwithmean(U,gt,k);
    ACC(i) = res(1);
    NMI(i) = res(2);
    ARI(i) = res(7);
    F1(i) = res(4);
end

mean_acc(anchor_ind,alpha_ind)= mean(ACC);
mean_nmi(anchor_ind,alpha_ind) = mean(NMI);
mean_ari(anchor_ind,alpha_ind) = mean(ARI);
mean_f1(anchor_ind,alpha_ind) = mean(F1);
std_acc(anchor_ind,alpha_ind)= std(ACC);
std_nmi(anchor_ind,alpha_ind)= std(NMI);
std_ari(anchor_ind,alpha_ind) = std(ARI);
std_f1(anchor_ind,alpha_ind) = std(F1);

toc
end

end
