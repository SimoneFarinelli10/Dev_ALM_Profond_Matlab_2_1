function xlsrange = GetExcelRange(startrow, endrow, startcol, endcol)
%GetExcelRange creates an excel range for Matlab import in the format
%required by xlsread

%helper function
function colstring = GetExcelColumn(colindex)
   colstring = dec2base(colindex-1, 26);
   digit = colstring < 'A';
   colstring(digit) = colstring(digit) + 'A' - '1';
   colstring(~digit) = colstring(~digit) + 9;
   colstring(end) = colstring(end) + 1;
end

xlsrange = sprintf('%s%d:%s%d', GetExcelColumn(startcol), startrow, GetExcelColumn(endcol), endrow);
      
end



