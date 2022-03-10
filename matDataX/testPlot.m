h2 = xMATERIAL('name', 'H2');
h2.meanFreePath;

colorMap = parula(5);
pressure = [0.1 0.5 1 5 10 50];

for i = 1:5
    MFP(i,:) = h2.meanFreePath('T',h2.TPlot, 'p', pressure(i));
    hold on
    plot(h2.TPlot, MFP(i,:),'Color',colorMap(i,:),'Linewidth',2)
end
hold off
xlabel('Temperature T / K')
ylabel('Mean free path / m')
title(strcat('Mean free path of ',h2.material))
set(gca,'linewidth',2)
lgd = legend(string(pressure));
lgd.Title.String = 'Pressure / bar';
legend('Location','best')