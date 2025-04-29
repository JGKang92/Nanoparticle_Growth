function f = figtile(varargin)
    TileSpacing = "tight";
    Padding = "tight";
    if nargin == 0
        f = figure();
        clf
        tiledlayout("flow", TileSpacing = TileSpacing, Padding = Padding)
    elseif nargin == 1
        if isscalar(varargin{1})
            if isnan(varargin{1})
                f = figure(NumberTitle = 'off');
            else
                f = figure(varargin{1});
            end
            clf
            tiledlayout("flow", TileSpacing = TileSpacing, Padding = Padding)
        else
            f = figure();
            clf
            row = varargin{1}(1);
            col = varargin{1}(2);
            tiledlayout(row, col, ...
                TileSpacing = TileSpacing, Padding = Padding)
        end
    elseif nargin ==2
        if isscalar(varargin{1})
            if isscalar(varargin{2})
                error("wrong argument format")
            end
            rowcol = varargin{2};
            row = rowcol(1);
            col = rowcol(2);
            figNum = varargin{1};
        elseif isscalar(varargin{2})
            if isscalar(varargin{1})
                error("wrong argument format")
            end
            rowcol = varargin{1};
            row = rowcol(1);
            col = rowcol(2);
            figNum = varargin{2};
        end
        if isnan(figNum)
            f = figure(NumberTitle = "off");
        else
            f = figure(figNum);
        end
        clf
        tiledlayout(row, col, TileSpacing = TileSpacing, Padding = Padding)
    else
        error("improper number of arguments\n")
    end
end