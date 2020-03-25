function code = formatCode(code)
    uniqueLable = unique(code);
    if size(uniqueLable,2) ~= size(code,2)
        for i = 1:size(code,2)
            for j = 1:size(uniqueLable,2)
                if code(1,i)==uniqueLable(1,j)
                    code(1,i) = j;
                end
            end
        end
    end
end