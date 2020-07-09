function nice_plot(varargin)
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
            'Legend1', [],...
            'Legend2', []);

    % Check input
    api = parseinputs(defapi, [], varargin{:});

    % Axes options
    a = (1:6)-0.45;
    b = (1:6)+0.45;
    d = 1e-02;

    % Get current axes
    ax = gca;

    % Set XAxis labels
    ax.XAxis.TickValues = [1 1.5 2 2.5 3];
    if api.XLabel
        xlabel(api.XLabelStr, 'interpreter', 'LaTeX');
        ax.XAxis.TickLabels = [1 1.5 2 2.5 3];
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
        ax.XMinorGrid = 'on';
        ax.YMinorGrid = 'off';
        ax.MinorGridLineStyle = '-';
        ax.MinorGridAlpha = 0.1;
        ax.GridAlpha = ax.MinorGridAlpha;
    ax.Box = 'on';
    ax.FontWeight = 'bold';
    ax.FontSmoothing = 'on';
    ax.FontSize = api.AxesFontSize;
    ax.LineWidth = api.AxesLineWidth;
    ax.XAxis.MinorTick = 'off';
    ax.TickLength = [0.025, 0.25];

        ax.XAxis.TickLength = [0.0, 0.25];
        a = (1:0.5:3)-0.225;
        b = (1:0.5:3)+0.225;
        d = 1e-02;
        ax.XAxis.MinorTickValues = [a(1):d:b(1), a(2):d:b(2),...
                                    a(3):d:b(3), a(4):d:b(4),...
                                    a(5):d:b(5)];
    ax.TickDir = 'in';
    ax.TickLabelInterpreter = 'latex';
    axis(api.Axis);

    % Legends
    h1 = api.Legend1;
    h2 = api.Legend2;
%     l1 = legend(h1,sprintf('s=%.0f mm',8),...
%                    sprintf('s=%.0f mm',10),...
%                    sprintf('s=%.0f mm',12),...
%                    sprintf('s=%.0f mm',14),...
%                    sprintf('s=%.0f mm',16));
%     l1.NumColumns = 5;
%     l1.Location = 'northwest';
%     l1.FontSize = 0.5*api.LegendFontSize;
%     l1.Interpreter = 'LaTeX';

    l2 = legend(ax,h2,'HARP','SinMod','DENSEanalysis');
    l2.Location = 'northwest';
    l2.LineWidth = 1.5;
    l2.FontSize = api.LegendFontSize;
    l2.Interpreter = 'latex';

end

