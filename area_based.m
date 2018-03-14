function [disp_map] = area_based(image1, image2, kernel, method, max_disp)
    %kernel should be odd size
%     assert(sum(size(kernel)/2 ~= 0)==2);
    image1 = cast(image1,'double');
    image2 = cast(image2,'double');
    assert(strcmp(method,'sd') |strcmp(method,'ad') |strcmp(method,'ncc'))
    disp_map = zeros(size(image1,1), size(image,2));
    trim_cols = (size(kernel,2)-1)/2;
    trim_rows = (size(kernel,1)-1)/2;
    for col1 = trim_cols+1:size(image1,2)-trim_cols-1
        for row1 = trim_rows+1:size(image1,1)-trim_rows-1
            pos_x=trim_cols+1;pos_y=row1;
            diff = inf;
            if strcmp(method,'ncc')
                diff = -inf;
            end
            win1 = image1(row1-trim_rows:row1+trim_rows,col1-trim_cols:col1+trim_cols).*kernel;
            for col2 = max(trim_cols + 1, col1-max_disp): min(size(image1,2) - (trim_cols + 1), col1+max_disp)
                win2 = image2(row1-trim_rows:row1+trim_rows,col2-trim_cols:col2+trim_cols).*kernel;
                switch(method)
                    case 'sd'
                        diff_curr = sd(win1,win2);
                        if diff_curr < diff
                            diff = diff_curr;
                            pos_x = col2;
                        end
                    case 'ad'
                        diff_curr = ad(win1,win2);
                        if diff_curr < diff
                            diff = diff_curr;
                            pos_x = col2;
                        end
                    case 'ncc'
                        diff_curr = ncc(win1,win2);
                        if diff_curr > diff
                            diff = diff_curr;
                            pos_x = col2;
                        end
                end
            end
            disp_map(row1,col1) = abs(col1 - pos_x) ;
        end
    end
    disp_map = cast(disp_map, 'uint8');
end

