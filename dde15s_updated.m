function sol = dde15s_updated(ddefun,lags,history,tspan,options,varargin) 
% Modification of dde23 to solve stiff DDEs by the method of steps
% with ode15s as the integrator, using its solution structure output.
%
% At present, some features of dde23 are not available:
%
% ** There is no provision for history in the form of a structure,
%    hence continuation is not available.
% ** Event location is not present.
% ** Stats are not aggregated.
%
% This code was written by Dr. Lawrence F. Shampine (Southern Methodist 
% University, Dallas, TX) for the manuscript: 
% Agrawal, Vikas, et al. "A dynamic mathematical model to clarify signaling 
% circuitry underlying programmed cell death control in Arabidopsis disease 
% resistance." Biotechnology progress 20.2 (2004): 426-442.
% $$$ % The original code works with MATLAB 6.5 (R13). In particular, it 
% $$$ % uses ODE15S Revision 1.83.  Some modifications to SOLEXTEND may be 
% $$$ % required to work with other versions of ODE15S.
%
% The code was updated to run in MATLAB R2018b by Jacek Kierzenka
% (Mathworks) in December 2020.

% Check inputs.
if nargin < 4
  error('Not enough input arguments. Must specify DDEs, lags, history, and tspan.');
elseif nargin == 4
  options = [];
end

t0 = tspan(1);
tfinal = tspan(end);   % Ignore all entries of tspan except first and last.
if tfinal <= t0
  error('Must have tspan(1) < tspan(end).')
end

if isnumeric(history)
  temp = history;
else
  temp = feval(history,t0,varargin{:});
end 
y0 = temp(:);
maxlevel = 4;
initialy = ddeget(options,'InitialY',[],'fast');
if ~isempty(initialy)
  y0 = initialy(:);
  maxlevel = 5;
end

t = t0;
y = y0;
neq = length(y);

% If solving a DDE, locate potential discontinuities. We need to step to each of
% the points of potential lack of smoothness. Because we start at t0, we can
% remove it from discont.  The solver always steps to tfinal, so it is
% convenient to add it to discont.
if isempty(lags)
  discont = tfinal;
  minlag = Inf;
else
  lags = lags(:)';
  minlag = min(lags);
  if minlag <= 0
    error('The lags must all be positive.')
  end
  vl = t0;
  maxlag = max(lags);
  jumps = ddeget(options,'Jumps',[],'fast');
  if ~isempty(jumps)
    indices = find( ((t0 - maxlag) <= jumps) & (jumps <= tfinal) );
    if ~isempty(indices)
      jumps = jumps(indices);
      vl = sort([vl jumps(:)']);
      maxlevel = 5;
    end
  end
  jumps = ddeget(options,'Jumps',[],'fast');
  if ~isempty(jumps)
    indices = find( ((t0 - maxlag) <= jumps) & (jumps <= tfinal) );
    if ~isempty(indices)
      jumps = jumps(indices);
      vl = sort([vl jumps(:)']);
      maxlevel = 5;
    end
  end
  discont = vl;
  for level = 2:maxlevel
    vlp1 = vl(1) + lags;
    for i = 2:length(vl)
      vlp1 = [vlp1 (vl(i)+lags)];
    end % Restrict to tspan.
    indices = find(vlp1 <= tfinal);
    vl = vlp1(indices);
    if isempty(vl)
      break;
    end
    nvl = length(vl);
    if nvl > 1 % Purge duplicates in vl.
      vl = sort(vl);
      indices = find(abs(diff(vl)) <= 10*eps*abs(vl(1:nvl-1))) + 1;
      vl(indices) = [];
    end
    discont = [discont vl];
  end
  if length(discont) > 1
    discont = sort(discont); % Purge duplicates.
    indices = find(abs(diff(discont)) <= 10*eps*abs(discont(1:end-1))) + 1;
    discont(indices) = [];
  end
end

% Add tfinal to the list of discontinuities if it is not already included.
if abs(tfinal - discont(end)) <= 10*eps*abs(tfinal)
  discont(end) = tfinal;
else
  discont = [discont tfinal];
end

% Discard t0 and discontinuities in the history.
indices = find(discont <= t0);
discont(indices) = [];

% Add t0 to the list of discontinuities:
discont = [t0 discont];

% Special cases to get going.
sol.x = t0;                  
sol.y = y0;
for i = 2:length(discont)
  distance = discont(i) - discont(i-1);
  nsteps = ceil(distance/minlag);
  stepsize = distance/nsteps;
  e = discont(i-1);
  for j = 1:nsteps
    b = e;
    e = b + stepsize;
    % Hit the discontinuity exactly.
    if j == nsteps
      e = discont(i);
    end
    solex = ode15s(@ddefcn,[b,e],sol.y(:,end),options,...
                   ddefun,lags,history,sol,varargin{:});
    if (i == 2) & (j == 1)    % Initialize solution structure.
      sol = solex;
    else
      sol = solextend(sol,solex);
    end
  end
end
sol.discont = discont;

%---------------------------------------------------------------------------

function solout = solextend(sol,solex)
% Concatenate two solution structures produced by ODE15S.
% $$$ % This function works with ODE15S from MATLAB 6.5, but 
% $$$ % it may require modifications to work with other versions.
% The code was updated to run in MATLAB R2018b

solout.solver = 'ode15s';
solout.x = [sol.x, solex.x(2:end)];
solout.y = [sol.y, solex.y(:,2:end)];

% This code was updated to run in MATLAB R2018b
% $$$ solout.idata.kvec = [sol.idata.kvec; solex.idata.kvec(2:end)];  
% $$$ dim3   = size(  sol.idata.dif3d,3);
% $$$ dim3ex = size(solex.idata.dif3d,3);
% $$$ if dim3ex > dim3
% $$$   sol.idata.dif3d(:,:,dim3+1:dim3ex) = 0;
% $$$   dim3 = dim3ex;
% $$$ elseif dim3ex < dim3
% $$$   solex.idata.dif3d(:,:,dim3ex+1:dim3) = 0;       
% $$$ end
% $$$ solout.idata.dif3d = [sol.idata.dif3d; solex.idata.dif3d(2:end,:,:)];
solout.idata.kvec = cat(2,sol.idata.kvec, solex.idata.kvec(2:end));  
[s1,s2,s3] = size(sol.idata.dif3d);
[~,s2n,s3n] = size(solex.idata.dif3d);
dataType = superiorfloat(sol.x,solex.x);
solout.idata.dif3d = zeros(s1,max(s2,s2n),s3+s3n+1-2,dataType);
solout.idata.dif3d(:,1:s2,1:s3) = sol.idata.dif3d;
solout.idata.dif3d(:,1:s2n,s3+1:end) = solex.idata.dif3d(:,:,2:end);
solout.extdata = [];

%---------------------------------------------------------------------------

function v = ddefcn(t,y,ddefun,lags,history,sol,varargin)
if isempty(lags)
  Z = [];
else
  nlags = length(lags);
  Z = zeros(length(y),nlags);
  for i = 1:nlags
    tmlag = min(t - lags(i),sol.x(end));  
    if tmlag <= sol.x(1)   
      if isnumeric(history)
        Z(:,i) = history(:);
      else
        Z(:,i) = feval(history,t,varargin{:});
      end
    else
      Z(:,i) = deval(sol,tmlag);
    end
  end
end
v = feval(ddefun,t,y,Z,varargin{:});

            
