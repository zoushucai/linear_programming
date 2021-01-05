function tab = nma_simplex(A,b,f)

[A,b] = make_phase_one(A,b);
tab   = simplex(A,b,f)
print_solution_vector(tab);
end
%==========================
function [A,b] = make_phase_one(A,b)
[m,n]              = size(A);
tab                = zeros(m+1,n+m+1);
tab(1:m,1:n)       = A;
tab(end,n+1:end-1) = 1;
tab(1:m,end)       = b(:);
tab(1:m,n+1:n+m)   = eye(m);
fprintf('>>>> Current Tableau \n');
disp(tab);

    for i = 1:m     %now make all entries in bottom row zero
        tab(end,:) = tab(end,:)-tab(i,:);
    end
tab = simplex(tab(1:m,1:n+m),tab(1:m,end),tab(end,1:n+m));
A = tab(1:m,1:n);
b = tab(1:m,end);
end
%=================================
function tab = simplex(A,b,f)
[m,n]        = size(A);
tab          = zeros(m+1,n+1);
tab(1:m,1:n) = A;
tab(m+1,1:n) = f(:);
tab(1:m,end) = b(:);

keep_running = true;
while keep_running
    
    fprintf('***********************\n');
    fprintf('Current Tableau \n');
    disp(tab);
    
    
    if any(tab(end,1:n)<0)%check if there is negative cost coeff.
        [~,J] = min(tab(end,1:n)); %yes, find the most negative
        % now check if corresponding column is unbounded
        if all(tab(1:m,J)<=0)
            error('problem unbounded. All entries <= 0 in column %d',J);
            %do row operations to make all entries in the column 0
            %except pivot
        else
            pivot_row = 0;
            min_found = inf;
            for i = 1:m
                if tab(i,J)>0
                    tmp = tab(i,end)/tab(i,J);
                    if tmp < min_found
                        min_found = tmp;
                        pivot_row = i;
                    end
                end
            end
            
            fprintf('Pivot row is %d\n',pivot_row);   
            %normalize
            tab(pivot_row,:) = tab(pivot_row,:)/tab(pivot_row,J);
            %now make all entries in J column zero.
            for i=1:m+1
                if i ~= pivot_row
                    tab(i,:)=tab(i,:)-sign(tab(i,J))*...
                        abs(tab(i,J))*tab(pivot_row,:);
                end
            end
        end
        %print current basic feasible solution
        fprintf('Current basic feasible solution is\n');
        disp(get_current_x());
        
    else
        keep_running=false;
    end
end

%internal function, finds current basis vector
    function current_x = get_current_x()
        current_x = zeros(n,1);
        for j=1:n
            if length(find(tab(:,j)==0))==m
                idx= tab(:,j)==1;
                current_x(j)=tab(idx,end);
            end
        end
    end
end

%=================================
function print_solution_vector(tab)
%tab(1:m,1:n) = A;

[nRow,nCol] = size(tab);

A = tab(1:nRow-1,1:nCol-1);
b = tab(1:nRow-1,nCol);

q = A ~= 0;
q = find(sum(q,1)==1); %find all columns with one non-zero entry;

solution_vector = zeros(nCol-1,1);

for n=1:length(q)
    j = find(A(1:nRow-1,q(n))==1);
    if isempty(j)
        solution_vector(q(n)) = 0;
    else
        solution_vector(q(n)) = b(j);
    end
end

solution_vector

end