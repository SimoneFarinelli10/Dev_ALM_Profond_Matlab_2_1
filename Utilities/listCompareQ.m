function result = listCompareQ(list1, list2)
%listCompareQ compares two lists given as cell array of strings
 result = logical(prod([ismember(list1, list2);ismember(list2, list1)]));
end

