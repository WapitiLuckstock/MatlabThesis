function A = turnHalfHourtoHourBySumming(B)
for i = 1:2:48
    for j = 1:365
    A(j,i) = B(j,i)+B(j,i+1);
    end
end

for i = 2:24
    A(:,i) = [];
end

end