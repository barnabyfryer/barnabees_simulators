function [F,k_global] = Fix_x(F,k_global,nodes)

%Take all rows corresponding to the x updates
Rows = nodes*2-1;
%Set these rows to be zero
k_global(Rows,:) = 0;
%Set corresponding force to be zero too
F(Rows,1) = 0;

%Add a one in the diagonal so that the update is forced to be du = 0
for i = 1:length(nodes)
    k_global(Rows(i),Rows(i)) = 1;
end

end

