function wrapthisfigure(this, figurename, wrap, captionname, figurecomment, landscape, doJPG, noWrap)
%--------------------------------------------------------------
% Prints the current figure to file 'figurename' as fig, eps, jpg and pdf.
% inserts into wrap
%--------------------------------------------------------------
% function wrapthisfigure(this, figurename, wrap, captionname, figurecomment, landscape)

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 4
   captionname = [];
end
if nargin < 5
   figurecomment = [];
end
if nargin < 6
   landscape = false;
end
if nargin < 7
   doJPG = false;
end
if nargin < 8
   noWrap = false;
end

if ~isempty(figurename)
    % replace any "." by "-" in figurename, otherwise, print will not append
    % the appropriate file name extension
    figurename = strrep(figurename, '.', '-');
end

if ~isempty(captionname)
   set(this, 'name', captionname)
elseif ~isempty(figurename)
   set(this, 'name', figurename)
end

if isempty(wrap) % do nothing, except changing the name of the figure as above
    return
end

% set(gcf, 'Renderer', 'painters') % to fix date axis bug

if nargin > 1 && ~isempty(wrap)
    if ~noWrap
        if (wrap.id ~= 0)
            if landscape
                latexwrapper(wrap, 'add', 'figure', figurename, captionname, figurecomment);
            else
                latexwrapper(wrap, 'add', 'sidewaysfigure', figurename, captionname, figurecomment);
            end
        else
            warning('em:msg', 'Cannot append figure to latexwrapper since wrap file is closed (fid=0)')
        end
    end
   if isfield(wrap, 'dir')
      figurename = fullfile(wrap.dir, figurename);
   end
end


if landscape
   orient(this, 'landscape')
end

drawnow

if doJPG
    print(this, '-djpeg', '-r500', figurename);
else
    print(this, '-deps', '-r300', '-loose', figurename);
end
