function [ nmi ] = normalized_mutual_information( clu_assignment1,clu_assignment2 )
%NMI Compute the Normalized Mutual Information of the two input partitions.
%clu_assignment2=decode(clu_assignment2);
%% Compute the confusion matrix
num_nodes = length(clu_assignment1);
clu_num1 = max(clu_assignment1);
clu_num2 = max(clu_assignment2);
cmat = zeros(clu_num1,clu_num2); % confusion matrix.
for i = 1:clu_num1
    for j = 1:clu_num2
        index1 = find(clu_assignment1 == i);
        index2 = find(clu_assignment2 == j);
        num1 = length(index1);
        num2 = length(index2);
        numtemp = 0;
        for m = 1:num1
            for n = 1:num2
                if index1(m) == index2(n)
                    numtemp = numtemp + 1;
                end
            end
        end       
        cmat(i,j) = numtemp;
    end
end


%% =============Compute the Normalized Mutual Information==============%

%% =============the numerator part=================%
nmi_numerator = 0;
for i = 1:clu_num1
    for j = 1:clu_num2
        if cmat(i,j) ~= 0
            nmi_numerator = nmi_numerator + ...
            cmat(i,j) * ...
            log10( (cmat(i,j)*num_nodes) / (sum(cmat(i,:))*sum(cmat(:,j))) );
        end       
    end
end
nmi_numerator = -2 * nmi_numerator;


%% ============the denominator part================%
nmi_denominator1 = 0;
nmi_denominator2 = 0;
for i = 1:clu_num1
    nmi_denominator1 = nmi_denominator1 + sum(cmat(i,:)) * log10( sum(cmat(i,:)) / num_nodes );
end
for j = 1:clu_num2
    nmi_denominator2 = nmi_denominator2 + sum(cmat(:,j)) * log10( sum(cmat(:,j)) / num_nodes );
end
nmi_denominator = nmi_denominator1 + nmi_denominator2;

%% final NMI
nmi = nmi_numerator/nmi_denominator;

%====================================================================%
end

