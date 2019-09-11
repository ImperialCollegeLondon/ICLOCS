function b=callback(t,f,x)

% t is the current iteration of the algorithm. 
% f is the current value of the objective 
%  x contains the extra information  passed  to ipopt.
% The default @callback is empty.

% Example of a possible callaback 
  global sol
  sol.iterates=t;
  fprintf('%3d  %0.3g \n',t,f);
  b = true;