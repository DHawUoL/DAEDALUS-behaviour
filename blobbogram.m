function blobbogram( filename,varargin )
%This function takes a text file and produces a Forest Plot or Blobbogram.
%The input should be a tab separated text file with columns name, OR, L95
%and U95, where name will be used as the tick labels on the axis, OR is the
%point estimate, and L95 and U95 are the 95% confidence intervals.
%The output is stored in a publication-ready png format and an editable
%.fig format.

%A row can be defined as a title in the original input text file by a row
% that contains string 1 NaN NaN
%  optional parameters:

% 'scale' takes inputs 'linear' and 'log' for the x
%  axis. Default is linear

%  'title' allows for custom title on figure. Defaults to filename
%
%   Usage: blobbogram('GORD') will take the association
%   analysis in GORD.txt, generate a Forest Plot, and store it
%   in GORD_ForestPlot.png, which should be publication-ready, and
%   GORD_ForestPlot.png for minor readjustments in MATLAB's GUI.

%   Harry Green (2018) Genetics of Complex Traits

%% Parsing optionals

p = inputParser;
defaultscale = 'linear';
expectedscale = {'linear','log'};
defaulttitle = filename;
defaultpagesize = [25 20];

addRequired(p,'filename');
addParameter(p,'scale',defaultscale,...
                 @(x) any(validatestring(x,expectedscale)));
addParameter(p,'title',defaulttitle);
             
parse(p,filename,varargin{:});

filename=p.Results.filename;
scale=p.Results.scale;
graphtitle=p.Results.title;


%%

%First read in data and store it in a table T
opts = detectImportOptions( strcat(filename),'NumHeaderLines',0, 'delimiter', '\t');
T = readtable(strcat(filename),opts);

hold on % Lots of plots needed, need to add TWO lines for each title, except when one of the titles is on line 1 of the table
h=height(T)+sum(isnan(T.U95))+1;
if isnan(T.U95(1))
    h=h-1;
end

plot([1,1],[0,h],':k') %first plot a vertical line for x=1
ylim([0 h]) %this allows for space at top and bottom

variables=T.name; %this is used as y tick labels later

j=h-1; %j is the y axis level at which the lines are plotted, which starts at j=h (because that's the top line and it works down)
for i=1:height(T)
    if isnan(T.U95(i))
        match=find(strcmp(variables,T.name(i)));
        variables={variables{1:match-1},' ',variables{match:end}}'; %if U95 is NAN (line i of table is a header) move j down one to leave space and add a space to the variable vector so there's a gap in the tick labels above the header
        if i~=1
            j=j-1;
        end
    else
        plot([-10^5,eps,10^5],[j,j,j],'Color',[.8 .8 .8],'LineWidth',.25) %plot guide line
        plot([T.L95(i),T.U95(i)],[j,j],'b','LineWidth',1) %plot 95% CI
        %plot(T.OR(i),i,'b.','MarkerSize',20) %plot point estimate
        plot(T.OR(i),j,'bd','MarkerFaceColor','b') %plot point estimate
    end
    j=j-1;
end

variables=flip(variables); %whoops

box on

set(gca,'TickLabelInterpreter', 'latex'); %LaTeX font used for figure
yticklabels(variables); %this reads the name column from the input file for the tick labels
yticks(1:h-1);

title(graphtitle,'Interpreter','latex')
%set(gca, 'YGrid', 'on', 'XGrid','off') %the YGrid is helpful for matching axis labels with plot objects

xlabel('Odds Ratio','Interpreter','latex') %more LaTeX

if strcmp(scale,'linear')
    xlim([min(T.L95)-0.1*abs(min(T.L95)) max(T.U95)+0.1*abs(min(T.U95))])
end

if strcmp(scale,'log')
    set(gca,'XScale','log')
    xlim([0.9*min(T.L95) 1/0.9*max(T.U95)]) %this allows for space at top and bottom
%     n=-5:5; %allows for axis labels from 10^-5:10^5
%     s=sort([10.^n, 1.2*10.^n, 1.5*10.^n, 2*10.^n, 3*10.^n, 6*10.^n]);
%     s(s<0.75*min(T.L95))=[];
%     s(s>1.33*max(T.U95))=[];
    a=-10:0.2:10;
    ea=exp(a);
    ea=sort([0.25,0.5,2,3,ea]);
    xticks(ea);
    xticklabels(string(round(ea,2)))
end

%Save the MATLAB fig
figname=[filename,'_ForestPlot.fig'];
savefig(figname);

%Save the png
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 20 12])
print([filename,'_ForestPlot.png'],'-dpng')

