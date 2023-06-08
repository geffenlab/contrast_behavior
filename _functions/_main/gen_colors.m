function colors = gen_colors(sd,noisev,minv,maxv,n)

if nargin < 5
    n = 7;
end

colors = zeros(n,3);
if strcmp(sd,'lohi')
    % red
    colors(1,:) = repmat(noisev,1,3);
    colors(2:end,:) = [[ones(n-2,1); .8] repmat(linspace(minv,maxv,n-1)',1,2)];
    
else
    % blue
    colors(1,:) = repmat(noisev,1,3);
    colors(2:end,:) = [repmat(linspace(minv,maxv,n-1)',1,2) [ones(n-2,1); .8]];
end