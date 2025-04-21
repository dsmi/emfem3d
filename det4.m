function d = det4(a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3,d4)
% d = det4(a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4)
%
% Calculate determinants of
%  | a1 a2 a3 a4 |
%  | b1 b2 b3 b4 |
%  | c1 c2 c3 c4 |
%  | d1 d2 d3 d4 |
% where the a1 a2 and the other arguments can be vectors.
%

% Determinant of
%  | b2 b3 b4 |
%  | c2 c3 c4 |
%  | d2 d3 d4 |
dt1 = b2.*(c3.*d4-d3.*c4) - b3.*(c2.*d4-d2.*c4) + b4.*(c2.*d3-d2.*c3);

% Determinant of
%  | b1 b3 b4 |
%  | c1 c3 c4 |
%  | d1 d3 d4 |
dt2 = b1.*(c3.*d4-d3.*c4) - b3.*(c1.*d4-d1.*c4) + b4.*(c1.*d3-d1.*c3);

% Determinant of
%  | b1 b2 b4 |
%  | c1 c2 c4 |
%  | d1 d2 d4 |
dt3 = b1.*(c2.*d4-d2.*c4) - b2.*(c1.*d4-d1.*c4) + b4.*(c1.*d2-d1.*c2);

% Determinant of
%  | b1 b2 b3 |
%  | c1 c2 c3 |
%  | d1 d2 d3 |
dt4 = b1.*(c2.*d3-d2.*c3) - b2.*(c1.*d3-d1.*c3) + b3.*(c1.*d2-d1.*c2);

% The final result
d = ( a1.*dt1 - a2.*dt2 + a3.*dt3 - a4.*dt4 );
