function C_out = relabel_cluster_by_cardinality(C_in)
% Relabel Clusters by Cardinality
% So that clusters with smaller index are larger in cardinality

labels      = unique(C_in);
cardinality = zeros(size(labels));
C_out       = zeros(size(C_in));

for i = 1:length(labels)
    cardinality(i) = sum(C_in == labels(i));
end

[~, idxSort] = sort(cardinality, 'descend');
labels = labels(idxSort);
for i = 1:length(labels)
    loc = C_in == labels(i);
    C_out(loc) = i;
end