%% binomial random draw function optimized for speed.

function [ res ] = mybinornd( N, p )
[row_cnt, col_cnt] = size(N);
res = zeros(row_cnt, col_cnt);
for ii=1:row_cnt
   for jj=1:col_cnt
      res(ii, jj) = sum(rand(1,N(ii,jj))<p(ii,jj));
   end
end
end