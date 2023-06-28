function [wnew, index] = resamplingIndex(weights,N)

Ns = length(weights); 

if nargin ==1
  N = length(weights);
end

index = zeros(1,N);
w = zeros(1,N);

c = cumsum(weights);

i = 1;
u = zeros(1,N);
u(1) = rand/N;


for j = 1:N
    u(j) = u(1) + (j-1)/N;
    
    while (u(j)>c(i))
        i = i + 1;
    end
    
    index(j) = i;
    
end

wnew = ones(1,N)/N;

end

