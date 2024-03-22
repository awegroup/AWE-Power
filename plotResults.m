function fig = plotResults(processedOutputs, yAxesOptions, subplotLayout, options)
    
    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'markerList', {'o','^','s','v','*','+','x','pentagram','diamond'});
    addParameter(p, 'LineWidth', 1);
    addParameter(p, 'MarkerSize', 6);
    addParameter(p, 'LineStyle', ':');

    if nargin<4
        parse(p,struct())
    else
        parse(p,options)
    end
    lineOptions = p.Results;

    if nargin < 3
        subplotLayout = [];
    end

    if isempty(subplotLayout)
        numSubplots = 1;
    else
        numSubplots = prod(subplotLayout);
    end

    fig=figure; hold on; grid on; box on

    if isfield(options,'colors')
        colororder(options.colors);
    end
    
    for i = 1:numSubplots
        if numSubplots==1
            subplot(1,1,1);
        else
            subplot(subplotLayout(1), subplotLayout(2), i);
        end
        plotSinglePlot(processedOutputs, yAxesOptions{i}, lineOptions)
    end

    if numSubplots>1
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='off';
        if isfield(options,'yLabel')
            ylabel(options.yLabel);
        end
        if isfield(options,'xLabel')
            xlabel(options.xLabel);
        end
        if isfield(options,'title')
            title(options.title);
        end
    end
end