function [ obj ] = getdiff( zAtrue, nsam, d )

diff = nan( d, d, nsam*(nsam-1)/2 );
count = 1;
for i = 1:(nsam-1)
    for j = (i+1):nsam
        diff(:,:,count) = abs( zAtrue(:,:,i) - zAtrue(:,:,j) );
        count = count + 1;
    end
end

obj = diff;

