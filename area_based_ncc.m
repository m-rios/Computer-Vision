function [disp_map] = area_based_ncc(image1, image2, kernel, max_disp)  
    image1 = rgb2gray(image1);
    image2 = rgb2gray(image2);
    image1 = cast(image1,'double');
    image2 = cast(image2,'double');
    
    disp_map = zeros(size(image1,1), size(image1,2));
    trim_cols = (size(kernel,2)-1)/2;
    trim_rows = (size(kernel,1)-1)/2;
    for col1 = trim_cols+1 +134 :size(image1,2)-trim_cols-1
        for row1 = trim_rows+1 +174:size(image1,1)-trim_rows-1
            win1 = image1(row1-trim_rows:row1+trim_rows,col1-trim_cols:col1+trim_cols);
            %win2 = image2(row1-trim_rows:row1+trim_rows, 1+trim_cols:size(image2,2)-trim_cols);
            c1 = max(trim_cols+1,col1-max_disp) ;
            c2 = min(size(image1,2)-(trim_cols+1),col1+max_disp);
            win2 = image2(row1-trim_rows:row1+trim_rows,c1:c2);
            c = normxcorr2(win1,win2);
            figure, surf(c), shading flat 
            [ypeak, xpeak] = find(c==max(c(:)));
            if c1 == (trim_cols+1)
                xdisp = abs(xpeak-trim_cols + trim_cols +1 - col1);
            else
                xdisp = abs(xpeak-trim_cols + col1-max_disp - col1);
            end
            disp_map(row1,col1) = min(xdisp);
        end
    end
    disp_map = cast(disp_map, 'double');
    disp_map = (disp_map - min(disp_map(:)))/(max(disp_map(:) - min(disp_map(:))));
end