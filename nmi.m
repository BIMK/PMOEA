
function MIhat = nmi( A, B )
%NMI Normalized mutual information
% http://en.wikipedia.org/wiki/Mutual_information
% http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
% Author: http://www.cnblogs.com/ziqiao/   [2011/12/15] 
if length( A ) ~= length( B)
    error('length( A ) must == length( B)');
end
total = length(A);
A_ids = unique(A);
A_class = length(A_ids);
B_ids = unique(B);
B_class = length(B_ids);
% Mutual information
idAOccur = double (repmat( A, A_class, 1) == repmat( A_ids', 1, total ));
idBOccur = double (repmat( B, B_class, 1) == repmat( B_ids', 1, total ));
idABOccur = idAOccur * idBOccur';
Px = sum(idAOccur') / total;
Py = sum(idBOccur') / total;
Pxy = idABOccur / total;
MImatrix = Pxy .* log2(Pxy ./(Px' * Py)+eps);
MI = sum(MImatrix(:));
% Entropies
Hx = -sum(Px .* log2(Px + eps),2);
Hy = -sum(Py .* log2(Py + eps),2);
%Normalized Mutual information
%MIhat = 2 * MI / (Hx+Hy);


  MIhat = MI / sqrt(Hx*Hy); %another version of NMI

end

% Example :  
% (http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html)
% A = [1 1 1 1 1 1   2 2 2 2 2 2    3 3 3 3 3];
% B = [1 2 1 1 1 1   1 2 2 2 2 3    1 1 3 3 3];
% nmi(A,B) 

% ans =  0.3646

