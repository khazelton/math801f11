function setdemorandstream(n)
%SETDEMORANDOMSTREAM Set default stream for reliable demo results.
%
%  SETDEMORANDOMSTREAM(N) is less distracting in demo code, but
%  equivalent to:
%
%    rs = RandStream('mcg16807','Seed',n);
%    RandStream.setGlobalStream(rs);

% Copyright 2011 The MathWorks, Inc.

rs = RandStream('mcg16807','Seed',n);
RandStream.setGlobalStream(rs);
