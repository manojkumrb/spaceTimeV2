%% ML_FigureSize - Standard Figure Size good for pdf and eps
%
% ML_FigureSize(Dimensione)
% Dimensione == 2 --> good for landscape orientation
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function ML_FigureSize(Dimensione)

if nargin==0; Dimensione=1; end;

set(gcf, 'Units', 'centimeters');
set(gcf, 'PaperPositionMode', 'auto');
if Dimensione==1;
    set(gcf,'PaperOrientation','portrait');
    set(gcf, 'Position',  [5   5   20.9920   13.6157]);
elseif Dimensione==2;
    set(gcf,'PaperOrientation','landscape');
    set(gcf, 'Position',  [5   2.5   25   16.2159]);    
end
set(gca, 'Units','normalized','Position',[0.15 0.2 0.75 0.7]);

