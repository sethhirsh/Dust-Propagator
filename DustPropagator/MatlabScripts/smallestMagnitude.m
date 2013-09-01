function [ result ] = smallestMagnitude(list1, list2)
%inputs two vectors of equal length
%compares pair of elements of the two vectors and stores the element of
%smallest magnitude in vector

%returns a vector of the same length as the two inputted vectors

finalList = zeros(size(list1));


for i = 1:size(list1)
    
    magnitude1 = abs(list1(i));
    magnitude2 = abs(list2(i));
    
    if magnitude1 >= magnitude2
        finalList(i) = list2(i);
    else
        finalList(i) = list1(i);
    end
    
    
    
end

    result = finalList;
end

