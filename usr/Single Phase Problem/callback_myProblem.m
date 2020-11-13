function b=callback_myProblem(t,f,data)

% t is the current iteration of the algorithm. 
% f is the current value of the objective 
% data contains the extra information  passed  to ipopt.
% global variable sol contain solution information regarding current iteration

% b returns true for solver to continue, retures false to terminate

% The default @callback is empty.

% Example of a possible callaback 
  global sol
  sol.iterates=t;
  fprintf('%3d  %0.3g \n',t,f);
  b = true;