function s = fredMDtcodeNote
% returns text explaining FRED-MD tcodes

s = 'The column tcode denotes the following data transformation for a series $x$: (1) no transformation; (2) $\Delta x_t$; (3) $\Delta^2 x_t$; (4) $\log(x_t)$; (5) $\Delta\log(x_t)$; (6) $\Delta^2\log(x_t)$; (7) $\Delta(x_t/x_{t-1} -1.0)$.';