function [original_result allresults mean_result blo bhi count] = bootstrap_for_vector(B, alpha, bootfun, varargin)

% This function deal with possible NaN value occassionally appear across data set
% The way it gets around this problem is that if the RESAMPLED dataset appear to have all NaN, 
% i.e. The whole column is NaN after resampling, that particular run is discarded 

% Since it will always discard the resampled value if they are all NaN, if the
% original data has a COLUMN that has all NaN. This program will not work.
% (That means you need to collect at least one valid point for all the samples)

% The first argument of varargin MUST be the ONLY data for resampling

% This function can deal with a 2D matrix. It views each row as an observation
% and use the function defined by the user to output a statistics that can
% also be in the form of a vector

% original_result is the mean of the result evaluated by the function
% allresults is just all the data evaluated by the function 
% mean_result is the mean using bootstrap
% blo is the lower bound adjusted by this algo
% bhi is the upper bound adjusted by this algo
% count gives the number of valid permutation used in this method after NaN is
% discarded.
%
% B is the number of permutations, e.g. 1000
% alpha is the confidence interval, e.g. 0.1
% bootfun is the function defined by the user. The result of this
% function is permutated.
% The first input of varargin MUST be the data used in the function defined
% by the user

% Example:
% [original_result allresults mean_result blo bhi count] = bootstrap_for_vector(1000,0.10, @impulse_response_for_bootci, all_in_complex(1:50,:,1)');

% written by Kai Sum Li and other members in Dr. Jonathan Simon's lab at
% University of Maryland, College Park
% contact: jzsimon@umd.edu

data = varargin{1};
if length(varargin) > 1
    for i = 2: length(varargin)
        input{i-1} = varargin{i};
    end
else
    input = {};
end
[row col] = size(data);
count = 1;


try  
    original_result = feval(bootfun,data,input{:});
catch 
    error('The following error occurred while trying to evaluate bootfun ''%s'':\n\n%s', func2str(bootfun),lasterr');
end

if (sum(sum(isnan(data)) ~= row) ~= col) 
    error('Observations provided has at least 1 column that is all NaN');
elseif ((sum(sum(isnan(data)) == (row-1))) > 0)
    warning('Observations provided has at least 1 column that has only 1 value. Resamping will throw away results and bca may not be able to calculate correctly');
end





%original_result = feval(bootfun,data,input{:});

for i = 1:B
    index = ceil(row.*rand(row,1));
    permutated = data(index,:);
    % Checking if resample observation have all NaN. If so, resample it
    % again
    if (sum(sum(isnan(permutated)) ~= row) == col)  
        allresults(count,:) = feval(bootfun,permutated,input{:});
        count = count+1;
    end
end
count = count-1;

mean_result = nanmean(allresults);
  
 
 
 
% Computing Bootstrap Confidence Interval using bca
% Martinez W.L. and Martinez A.R.
% Computational Statistics Handbook with Matlab 2002 pp247-249
 
[Arow Acol] = size(allresults);
numberofresultlessthanoriginal = zeros(1,Acol);
cc = 1;
 
for i = 1:Arow
   numberofresultlessthanoriginal = numberofresultlessthanoriginal + (allresults(i,:) < original_result);
end
z0 = norminv(numberofresultlessthanoriginal./(count),0,1); 
 
 
% ahat = zeros(1,102);
zlo = norminv(alpha/2,0,1);  % this is the z^(a/2)
zup = norminv(1-alpha/2,0,1);  % this is the z^(1-a/2)
