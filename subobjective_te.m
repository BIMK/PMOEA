function [ obj ] = subobjective_te( weight, objectives, idealpoint )
%SUBOBJECTIVE 
s = size(weight, 1);
objsize = size(objectives,1);
    
weight((weight == 0))=0.00001;
    
if objsize==s 
    part2 = abs(objectives-idealpoint(ones(objsize,1),:));
%     part2 = objectives;
    obj = max(weight.*part2,[],2);
elseif objsize ==1
    part2 = abs(objectives-idealpoint);
%  part2 = objectives;
    obj = max(weight.*part2(ones(s, 1),:),[],2);   
else
    error('individual size must be same as weight size, or equals 1');
end


end

