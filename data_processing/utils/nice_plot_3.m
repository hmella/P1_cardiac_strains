function nice_plot_3(varargin)
  %NICE_PLOT_MAG Summary of this function goes here
%   Detailed explanation goes here

    %% PARSE INPUT
    % Default arguments
    defapi = struct(...
            'AxesFontSize',  20,...
            'AxesLineWidth', 2,...
            'LegendFontSize', 13,...
            'LegendLocation', 'northwest',...
            'Axis', [],...
            'XLabel', true,...
            'YLabel', true,...
            'XLabelStr', 'displacement (in wavelengths)',...
            'YLabelStr', [],...
            'YAxisTickValues', 0:10:100,...
            'YAxisExponent', 0,...
            'LegendHandle', []);

    % Check input
    api = parseinputs(defapi, [], varargin{:});

    % Axes options
    a = (1:6)-0.45;
    b = (1:6)+0.45;
    d = 1e-02;

    % Get current axes
    ax = gca;

    % Set XAxis labels
    ax.XAxis.TickValues = [1 10 20];
    if api.XLabel
        xlabel(api.XLabelStr, 'interpreter', 'LaTeX');
        ax.XAxis.TickLabels = [1 10 20];
    else
        ax.XAxis.TickLabels = [];
    end

    % Set YAxis labels
    ax.YAxis.TickValues = api.YAxisTickValues;
    if api.YLabel
        ylabel(api.YLabelStr, 'interpreter', 'LaTeX')
    else
        ax.YAxis.TickLabels = [];
    end
    ax.YAxis.Exponent = api.YAxisExponent;
    ax.XGrid = 'off';
    ax.YGrid = 'off';
%         ax.XMinorGrid = 'on';
%         ax.YMinorGrid = 'off';
%         ax.MinorGridLineStyle = '-';
%         ax.MinorGridAlpha = 0.1;
%         ax.GridAlpha = ax.MinorGridAlpha;
    ax.Box = 'on';
    ax.FontWeight = 'bold';
    ax.FontSmoothing = 'on';
    ax.FontSize = api.AxesFontSize;
    ax.LineWidth = api.AxesLineWidth;
    ax.XAxis.MinorTick = 'off';
    ax.TickLength = [0.025, 0.25];

    ax.XAxis.TickLength = [0.025, 0.25];

    ax.TickDir = 'in';
    ax.TickLabelInterpreter = 'latex';
    axis(api.Axis);

%     l.Location = 'northeast';
%     l.LineWidth = 1.5;
%     l.FontSize = api.LegendFontSize;
%     l.Interpreter = 'latex';

end

