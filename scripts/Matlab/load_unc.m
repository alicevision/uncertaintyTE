function [ unc_cam, unc_pts ] = load_unc( path, ncams, npts )
    unc_cam = cell(ncams,1);
    unc_pts = cell(npts,1);

    fid = fopen(path,'r');
    frewind(fid);
    
    counter = 1;
    while ( true ) 
        % load line
        line_in = fgets(fid);
        if (line_in == -1)
            break;
        end
        % skip #
        if (line_in(1) == '#')
           continue; 
        end
        
        % parse line to a matrix
        X = str2num(line_in);
        n = ceil((sqrt(8*size(X,2))-1)/2);
        A = triu(ones(n));
        A = A';
        A(A~=0) = 1:n*(n+1)/2;
        C = reshape( X(A + A' - diag(diag(A))) , n, n );
            
        % save it to the output array
        if counter <= ncams     % cameras
            unc_cam{counter} = C;
            
        else                    % points
            unc_pts{counter - ncams} = C;
        end
        
        counter = counter + 1;
    end
    fclose(fid);
end