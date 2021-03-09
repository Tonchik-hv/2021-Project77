function html = demo(filename)
txtTsSize = [];
txtLoadData = [];
%%
if nargin == 0
    filename = 'walk_6.xml';    
end


ext = filename(end-2:end);
if strcmp(ext, 'xml')
    ts = getXmlData(filename(1:end-4));    
elseif strcmp(ext, 'mat')
    ts = importdata(filename);
else
    txtLoadData = 'The file extention is not ".mat", presenting demo version.';
end
%}
[n, m] = size(ts);
if n < m
    ts = ts';
end
if size(ts,2) > 1
    txtTsSize = 'The time series should be one-dimensional. Processing the first dimension';
    ts = ts(:,1);
end
%%
ts = center(ts);
figPCs = 'pc.png';
figLam = 'Lam.png';
[idxSegm, structT, txtSegmMsg] = calcSegms(ts, 'figname', {figPCs, figLam});


figSegm = 'segmTS.png';

plotPeriods(ts(1:min(length(idxSegm), 500)), idxSegm(1:min(length(idxSegm), 500)), [], figSegm);

%
html='<ul>';
if ~isempty(txtLoadData)
    html = [html, '<p>', txtLoadData, '</p>'];
end
if ~isempty(txtTsSize) && ~isempty(strfind(txtTsSize, 'short'))
    html = [html, '<p>', txtTsSize, '</p>'];
end
html = [html, '<table width = "400"> <col width = "200"><col width = "200">'];
html = [html, '<tr> <td ><p>', txtTsSize, 'The square roots of eigenvalues $\sqrt{\lambda_j}$. </p>', ...
    '<img src="' figLam '"  WIDTH="400" alt="Uploaded image"/> </td>'];
if ~isempty(txtSegmMsg)
    html = [html, '<td>', txtSegmMsg, '</td></tr></table> </ul>'];
else
html = [html, '<td> The red cicles mark the eigenvalues of the selected components </td></tr> </ul>'];
html = [html, '<tr><td><p> The trajetory of seleted principal components. </p>', ...
    '<img src="' figPCs '"  WIDTH="400" alt="Uploaded image"/></td>',...
    '<td> The blue line represents the trajectory of selected principal components,' ...
    'plotted againt each other. Red cicles denote the ending points of extracted periods. </td></tr>'];
html = [html, '<tr><td> <p> The figure shows the results of time series segmentation. </p>', ...
     '<img src="' figSegm '"  WIDTH="400" alt="Uploaded image"/></td>',... 
     '<td> The blue line is the historical time series, the red circles mark end points of the periods. </td></tr></table></ul>'];
end

fid = fopen('out.html', 'w');
fprintf(fid, '%s', html);
fclose(fid);
%}

end

function ts = center(ts)

centr = mean(ts);
ts = ts - repmat(centr, length(ts),1);

end
