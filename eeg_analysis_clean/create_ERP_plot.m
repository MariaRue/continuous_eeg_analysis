function create_ERP_plot(time,Avg,Se,colour,LineWidth)


h = shadedErrorBar(time,Avg,Se, 'lineprops', '-k');
h.patch.FaceColor = colour;
h.mainLine.Color = colour;
h.mainLine.LineWidth = LineWidth; 





end 