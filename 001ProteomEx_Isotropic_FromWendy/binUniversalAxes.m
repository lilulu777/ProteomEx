    measResults = imgPair.measResults;
    xmin = 4;%in pixels
    xmax = size(measResults,1);
    xvals = [xmin:xmax];
    notNaNindices = xvals(~isnan(measResults([xmin:xmax],3)));
    xvals = notNaNindices * imgPair.info.pixelWidthExpanded/imgPair.info.expFactor;
    edges = 0:.5:80;
    [bincounts] = histc(xvals,edges);
    prevcounts = 1;
    yvals_binned = []
    for i = 1:length(bincounts)-1;
        i
        notNaNindices(prevcounts:(prevcounts + bincounts(i) -1))
        yvals_binned(i) = sum(measResults(notNaNindices(prevcounts: (prevcounts + bincounts(i) -1)),3))/bincounts(i);
        prevcounts = bincounts(i) + prevcounts
    end
    yvals_binned = yvals_binned*imgPair.info.pixelWidthExpanded/imgPair.info.expFactor;