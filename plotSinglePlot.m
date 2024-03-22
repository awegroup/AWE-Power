function plotSinglePlot(processedOutputs, subplotOptions, lineOptions)
    % Input:
    % processedOutputs: Structure containing processed data
    % subplotOptions: Structure specifying options for a single subplot
    
    hold on; grid on; box on
    % Plot each data series
    for i = 1:size(subplotOptions,1)
        % Determine yyaxis position
        if i==1 && size(subplotOptions,1)==2
            yyaxis left
        elseif i==2
            yyaxis right
        end
        ylabel(subplotOptions(i).yLabel);
        % Plot data
        markerInd = 1;
        if isfield(subplotOptions(i),'xData')
            for j = 1:numel(subplotOptions(i).yData)
                if iscell(subplotOptions(i).xData)
                    if isnumeric(subplotOptions(i).xData{j})
                        x = subplotOptions(i).xData{j};
                    else
                        x = processedOutputs.(subplotOptions(i).xData{j});
                    end
                else
                    if isnumeric(subplotOptions(i).xData)
                        x = subplotOptions(i).xData;
                    else
                        x = processedOutputs.(subplotOptions(i).xData);
                    end
                end
                if isnumeric(subplotOptions(i).yData{j})
                    y = subplotOptions(i).yData{j};
                else
                    if isfield(subplotOptions(i),'yFcn')
                        if size(processedOutputs.(subplotOptions(i).yData{j}),1) == 1
                            y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j})');
                        else
                            y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j}));
                        end
                    else
                        y = processedOutputs.(subplotOptions(i).yData{j});
                    end
                end
                plot(x, y, 'LineStyle', lineOptions.LineStyle, 'Marker', lineOptions.markerList{markerInd}, ...
                    'LineWidth', lineOptions.LineWidth, 'MarkerSize', lineOptions.MarkerSize);
                markerInd = markerInd+1;
            end
        else
            for j = 1:numel(subplotOptions(i).yData)
                if isfield(subplotOptions(i),'windIndex')
                    ws = subplotOptions(i).windIndex;
                    if isnumeric(subplotOptions(i).yData{j})
                        y = subplotOptions(i).yData{j};
                    else
                        if isfield(subplotOptions(i),'yFcn')
                            y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j})(ws,:));
                        else
                            y = processedOutputs.(subplotOptions(i).yData{j})(ws,:);
                        end
                    end
                else
                    if isnumeric(subplotOptions(i).yData{j})
                        y = subplotOptions(i).yData{j};
                    else
                        if isfield(subplotOptions(i),'yFcn')
                            y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j}));
                        else
                            y = processedOutputs.(subplotOptions(i).yData{j});
                        end
                    end
                end
                plot(y, 'LineStyle', lineOptions.LineStyle, 'Marker', lineOptions.markerList{markerInd}, ...
                    'LineWidth', lineOptions.LineWidth, 'MarkerSize', lineOptions.MarkerSize);
                markerInd = markerInd+1;
            end
        end  
    end
    
    % Set x-axis label and legend
    if isfield(subplotOptions(i),'xLabel')
        xlabel(subplotOptions(i).xLabel);
    end
    if isfield(subplotOptions(i),'Legend')
        legendLabels = [subplotOptions.Legend];
        legend(legendLabels, 'location', 'northwest');
    end

    if isfield(subplotOptions(i),'xLim')
        xlim(subplotOptions(i).xLim);
    end
    if isfield(subplotOptions(i),'yLim')
        ylim(subplotOptions(i).yLim);
    end
    if isfield(subplotOptions(i),'title')
        title(subplotOptions(i).title);
    end
    hold off;
end