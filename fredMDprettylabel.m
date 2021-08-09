function pcode = fredMDprettylabel(ncode, doShort)
% FREDMDPRETTYLABEL ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 09-Dec-2019 18:55:22 $
% $Revision : 1.00 $
% DEVELOPED : 9.7.0.1247435 (R2019b) Update 2
% FILENAME  : fredmdPrettylabel.m


if nargin <2 || isempty(doShort)
    doShort = true;
end

pcode = cell(size(ncode));

if doShort
    for n = 1 : length(ncode)
        switch ncode{n}
            case 'RPI'
                pcode{n} = 'Real Income';
            case 'DPCERA3M086SBEA'
                pcode{n} = 'Real Consumption';
            case 'CMRMTSPLx'
                pcode{n} = 'Real M\&T Sales';
            case 'INDPRO'
                pcode{n} = 'IP';
            case 'CUMFNS'
                pcode{n} = 'Capacity Utilization';
            case 'UNRATE'
                pcode{n} = 'Unemployment Rate';
            case 'PAYEMS'
                pcode{n} = 'Nonfarm Payrolls';
            case 'CES0600000007'
                pcode{n} = 'Hours';
            case 'CES0600000008'
                pcode{n} = 'Hourly Earnings';
            case 'WPSFD49207'
                pcode{n} = 'PPI (Fin. Goods)';
            case 'PPICMM'
                pcode{n} = 'PPI (Metals)';
            case 'PCEPI'
                pcode{n} = 'PCE Prices';
            case 'CPIAUCSL'
                pcode{n} = 'CPI';
            case 'FEDFUNDS'
                pcode{n} = 'Federal Funds Rate';
            case 'HOUST'
                pcode{n} = 'Housing Starts';
            case {'S_P500', 'SP500'}
                pcode{n} = 'S\&P 500';
            case 'EXUSUKx'
                pcode{n} = 'USD / GBP FX Rate';
            case 'GS1'
                pcode{n} = '1-Year Yield';
            case 'GS5'
                pcode{n} = '5-Year Yield';
            case 'GS10'
                pcode{n} = '10-Year Yield';
            case 'BAAFFM'
                pcode{n} = 'Baa Spread';
            case 'CYCLE1S_ASSETS'
                pcode{n} = 'Cyc. Assets';
            otherwise
                pcode{n}= ncode{n};
        end % switch
    end % for n
else
    for n = 1 : length(ncode)
        switch ncode{n}
            case 'RPI'
                pcode{n} = 'Real Personal Income';
            case 'DPCERA3M086SBEA'
                pcode{n} = 'Real Personal Consumption Expenditures';
            case 'CMRMTSPLx'
                pcode{n} = 'Real Manu. and Trade Industries Sales';
            case 'INDPRO'
                pcode{n} = 'IP Index';
            case 'CUMFNS'
                pcode{n} = 'Capacity Utilization: Manufacturing';
            case 'UNRATE'
                pcode{n} = 'Civilian Unemployment Rate';
            case 'PAYEMS'
                pcode{n} = 'All Employees: Total nonfarm';
            case 'CES0600000007'
                pcode{n} = 'Avg Weekly Hours : Goods-Producing';
            case 'CES0600000008'
                pcode{n} = 'Avg Hourly Earnings : Goods-Producing';
            case 'WPSFD49207'
                pcode{n} = 'PPI: Finished Goods';
            case 'PPICMM'
                pcode{n} = 'PPI: Metals and metal products';
            case 'PCEPI'
                pcode{n} = 'Personal Cons. Expend.: Chain Index';
            case 'CPIAUCSL'
                pcode{n} = 'CPI';
            case 'FEDFUNDS'
                pcode{n} = 'Effective Federal Funds Rate';
            case 'HOUST'
                pcode{n} = 'Housing Starts: Total New Privately Owned';
            case {'S_P500', 'SP500'}
                pcode{n} = 'S\&P 500';
            case 'EXUSUKx'
                pcode{n} = 'USD / GBP Foreign Exchange Rate';
            case 'GS1'
                pcode{n} = '1-Year Treasury Rate';
            case 'GS5'
                pcode{n} = '5-Year Treasury Rate';
            case 'GS10'
                pcode{n} = '10-Year Treasury Rate';
            case 'BAAFFM'
                pcode{n} = 'Moody''s Baa Corporate Bond Minus FEDFUNDS';
            case 'CYCLE1S_ASSETS'
                pcode{n} = 'Cycle1s Assets';
            otherwise
                pcode{n}= ncode{n};
        end % switch
    end % for n
end % if doShort