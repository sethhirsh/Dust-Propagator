function [ result ] = smallestMagnitude(list1, list2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

    result = finalList
end

