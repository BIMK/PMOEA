function community =label2community(label)
    community = {};
    communityLength = length(unique(label));
    uniqueCommunity = unique(label);
    for i = 1:communityLength
       community{i} = find(label==uniqueCommunity(1,i));
    end
    
end