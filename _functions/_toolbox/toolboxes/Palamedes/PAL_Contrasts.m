%
%PAL_Contrasts  Linear, orthogonal set of contrasts
%
%   syntax: contrasts = PAL_Contrasts(N, {optional arguments} )
%
%   contrasts = PAL_Contrasts(N) where N is a scalar, returns an N x N 
%   model matrix containing an intercept term and, if they can be found,
%   N - 1 orthogonal polynomial contrasts.
%
%   Other sets of contrasts may be generated by passing the type of the
%   contrasts as an optional argument. Supported: 'Helmert', 'Periodic',
%   'Polynomial' and 'Identity'. The latter returns an identity matrix.
%
%   In the case of polynomial contrasts equally spaced values of the IV are
%   assumed. This may be overridden by supplying a third argument, which 
%   should be a vector specifying IV values (see example below). When a
%   vector is supplied that contains non-integers, orthogonality of the
%   matrix is not verified and not guaranteed (and a warning is issued 
%   stating this). The warning will include the largest absolute value of
%   the dot-products of all possible pairings of contrasts. If this number 
%   is real small (e.g., 1e-08) the detected non-orthogonality may simply
%   be rounding error. Note that any linear transformation on IV
%   levels will not affect the contrast coefficients. In other words, in 
%   case IV levels are non-integers, apply a linear transformation in order 
%   to create values that are integers if possible (compare examples 2 and 
%   3 below).
%
%   In the case of polynomial contrasts and larger values of N, a full set 
%   of orthogonal contrasts may not necessarily be found. In these cases, 
%   PAL_Contrasts returns only those contrasts that correspond to the 
%   highest degree polynomial that could be found. It also issues a warning 
%   that the contrast matrix is incomplete.
%
%   Example 1: PAL_Contrasts(3) returns:
%
%    1    1    1
%   -1    0    1
%    1   -2    1
%
%   Example 2: PAL_Contrasts(3,'Poly',[.2 .3 .5]) issues a warning and:
%
%    1.0000    1.0000    1.0000
%   -0.1333   -0.0333    0.1667
%    0.0086   -0.0129    0.0043
%
%   Example 3: PAL_Contrasts(3,'poly',10*[.2 .3 .5]) returns:
%
%      1     1     1
%     -4    -1     5
%      2    -3     1
%
%   Example 4: PAL_Contrasts(3,'Helmert') returns:
%
%    1    1    1
%    2   -1   -1
%    0    1   -1
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.6.2, 1.6.3 (see History.m)

function contrasts = PAL_Contrasts(N, varargin)

ws = warning('off','MATLAB:gcd:largestFlint'); %Results in non-orthogonality:
                                               %Palamedes issues its own
                                               %warning.

option = 'polynomial';

if ~isempty(varargin)
    valid = 0;
    if strncmpi(varargin{1}, 'polynomial',4)
        option = 'polynomial';
        valid = 1;
    end
    if strncmpi(varargin{1}, 'periodic',4)
        option = 'periodic';
        valid = 1;
    end
    if strncmpi(varargin{1}, 'helmert',4)
        option = 'helmert';
        valid = 1;
    end
    if strncmpi(varargin{1}, 'identity',4)
        option = 'identity';
        valid = 1;
    end
    if valid == 0
        warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{1});    
    end        
end            

switch lower(option(1:4))
    case{'iden'}
        contrasts = eye(N);
    case{'helm'}
        contrasts = zeros(N,N);
        contrasts(1,:) = 1;
        for index = 2:N
            x = zeros(1,N);
            x(index-1) = N-(index-1);
            x(index:N) = -1;
            contrasts(index,:) = x;
        end
    case{'peri'}
        contrasts = zeros(N,N);
        contrasts(1,:) = 1;
        for index = 2:N
            f = (index-1)/2;
            x = [0:2*pi/N:2*pi-2*pi/N] + pi/N;
            x = cos(f.*x);
            contrasts(index,:) = x;
        end
    case{'poly'}    
        levels = 1:N;
        Integer = 1;
        if length(varargin) == 2
            if isnumeric(varargin{2}) && size(varargin{2},1) == 1 && size(varargin{2},2) == N
                levels = varargin{2};
                if sum(levels == round(levels)) ~= length(levels);
                    Integer = 0;
                end
            else
                warning('PALAMEDES:invalidOption','''%s'' is not a numerical vector of length %s. Ignored.',num2str(varargin{2}),num2str(N));
            end
        end
        contrasts = ones(1,N);
        orthogonal = 1;
        index = 2;
        while (orthogonal == 1 || ~Integer) && index <= N
            x = (levels).^(index-1);
            temp = x;
            for w = 1:index-1   %Gram-Schmidt orthogonalization
                bN = sum(temp.*contrasts(w,:)); 
                bD = sum(contrasts(w,:).^2);
                temp = temp - (bN./bD).*contrasts(w,:);
                if Integer
                    temp = round(temp.*bD./gcd(bN,bD));
                    Divider = temp(1);
                    for coef = 2:N
                        Divider = gcd(Divider,temp(coef));
                    end
                    temp = temp./Divider;
                else
                    temp = temp./max(abs(min(temp)), max(temp));
                end
            end
            orthogonal = PAL_isOrthogonal([contrasts; temp]);
            if orthogonal || ~Integer
                contrasts = [contrasts; temp];
                index = index + 1;
            end
        end                
        if index - 1 < N
            warning('PALAMEDES:PAL_Contrasts:incompletePolynomials','complete set of orthogonal contrasts could not be found');
        end
        if ~Integer
            [orthogonal, nonOrtho] = PAL_isOrthogonal(contrasts);
            message = ['''levels'' contains numbers that are not integers. Returned contrasts may not be orthogonal. Largest absolute value of dot-product: ' num2str(max(abs(nonOrtho(:,3)))) '. Type help PAL_Contrasts for more information and a possible solution.'];            
            warning('PALAMEDES:PAL_Contrasts:mayNotBeOrthogonal',message);
        end
end

warning(ws);