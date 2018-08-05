function [out] = btc_func (infile,bx,by,outfile)
% Image Compression Using Block Truncation Coding.
%  infile is input file name
%  bx is a positive integer (x block size)
%  by is a positive integer (y block size)
%  outfile is output file name which will be created
%*****************************************************************

if (exist(infile)==2)
    a = imread(infile);
    figure('Name','Input image');
    imshow(a);
else
    warndlg('The file does not exist.',' Warning: ');
    out=[];
    return
end

    double_a=double(a);
    ax=size(a,1)-mod(size(a,1),bx);
    ay=size(a,2)-mod(size(a,2),by);
    out_rgb=zeros(ax,ay,3);

    % -----------------------  RED component  ------------------   
    dvalue=double_a(:,:,1);
    max_r = max(max(dvalue));
    min_r = min(min(dvalue));
    thre_r = (max_r+min_r)/2;
    dx=size(dvalue,1);
    dy=size(dvalue,2);
    % if input image size is not a multiple of block size image is resized
    modx=mod(dx,bx);
    mody=mod(dy,by);
    dvalue=dvalue(1:dx-modx,1:dy-mody);
    % the new input image dimensions (pixels)
    dx=dx-modx;
    dy=dy-mody;
    % number of non overlapping blocks required to cover the entire input image
    nbx=size(dvalue,1)/bx;
    nby=size(dvalue,2)/by;
    
    % the output compressed image
    matrice=zeros(bx,by);
    % the compressed data
    m_u=zeros(nbx,nby);
    m_l=zeros(nbx,nby);
    mat_log=logical(zeros(bx,by));
    
    posbx=1;
    for ii=1:bx:dx
        posby=1;
        for jj=1:by:dy
            % the current block
            blocco=dvalue(ii:ii+bx-1,jj:jj+by-1);
            
            % min and max values
            max1 = max(max(blocco));
            min1 = min(min(blocco));
            
            if (max1-min1)>=thre_r
                Ymax = max(blocco(:));
                Ymin = min(blocco(:));
                Th = (3*Ymax+Ymin)/4;
                Tl = (Ymax+3*Ymin)/4;
                p=0;
                q=0;
                for i=1:bx
                    for j=1:by
                       if blocco(i,j)>=Th
                           p=p+1;
                       elseif blocco(i,j)<=Tl
                           q=q+1;
                       end
                    end    
                end    

                sum1 = 0;
                sum2 = 0;
                for i = 1:bx
                    for j = 1:by
                        if blocco(i,j)>=Th
                            sum1 = sum1+(blocco(i,j))
                        elseif blocco(i,j)<=Tl
                            sum2 = sum2+(blocco(i,j))
                        end
                    end
                end  
                 M=0.5*((1/p)*sum1)+((1/q)*sum2); 
                 R = ((1/p)*sum1)-((1/q)*sum2);
                 T1=M+(-1/3)*R;
                 T2=M;
                 T3=M+(1/3)*R;
                % the logical matrix correspoending to the current block
                blocco_binario=zeros(bx,by);
                for i = 1:bx
                    for j = 1:by
                        if blocco(i,j)>=T3
                            blocco_binario(i,j)=3;
                        elseif blocco(i,j)>=T2 && blocco(i,j)<T3
                            blocco_binario(i,j)=2;
                        elseif blocco(i,j)>=T1 && blocco(i,j)<=T2
                            blocco_binario(i,j)=1;
                        else
                            blocco_binario(i,j)=0;
                        end
                    end
                end
                % the average gray level of pixels whose level is GREATER than the block average gray level
                K=sum(sum(double(blocco_binario)));
                mu=sum(sum(double(blocco_binario).*blocco))/K;
                % the average gray level of pixels whose level is SMALLER than the block average gray level
                matrice(ii:ii+bx-1,jj:jj+by-1)=double(M+(double(2*blocco_binario)-3)/6*R);            
                posby=posby+1;
            
            else
                %disp(min1);
                % the average gray level of the current block
                m=mean(mean(blocco));
                % the logical matrix correspoending to the current block
                blocco_binario=(blocco>=m);
                % the number of pixel (of the current block) whose gray level is greater than the average gray level of the current block
                K=sum(sum(double(blocco_binario)));
                % the average gray level of pixels whose level is GREATER than the block average gray level
                mu=sum(sum(double(blocco_binario).*blocco))/K;
                % the average gray level of pixels whose level is SMALLER than the block average gray level
                if K==bx*by
                    ml=0;
                else
                    ml=sum(sum(double(~blocco_binario).*blocco))/(bx*by-K);
                end
                % the COMPRESSED DATA which correspond to the input image
                m_u(posbx,posby)=mu;                            %---> the m_u matrix (see the cited reference)
                m_l(posbx,posby)=ml;                            %---> the m_l matrix 
                mat_log(ii:ii+bx-1,jj:jj+by-1)=blocco_binario;  %---> the logical matrix
                % the compressed image

                matrice(ii:ii+bx-1,jj:jj+by-1)=(double(blocco_binario).*mu)+(double(~blocco_binario).*ml);            
                posby=posby+1;
            end
            end
            posbx=posbx+1;
        end
    out_rgb(:,:,1)=matrice;

    % -----------------------  GREEN component  ------------------   
    dvalue=double_a(:,:,2);
    
    max_g = max(max(dvalue));
    min_g = min(min(dvalue));
    thre_g = (max_g+min_g)/2;
    dx=size(dvalue,1);
    dy=size(dvalue,2);
    % if input image size is not a multiple of block size image is resized
    modx=mod(dx,bx);
    mody=mod(dy,by);
    dvalue=dvalue(1:dx-modx,1:dy-mody);
    % the new input image dimensions (pixels)
    dx=dx-modx;
    dy=dy-mody;
    % number of non overlapping blocks required to cover the entire input image
    nbx=size(dvalue,1)/bx;
    nby=size(dvalue,2)/by;
    
    % the output compressed image
    matrice=zeros(bx,by);
    % the compressed data
    m_u=zeros(nbx,nby);
    m_l=zeros(nbx,nby);
    mat_log=logical(zeros(bx,by));
    
    posbx=1;
    for ii=1:bx:dx
        posby=1;
        for jj=1:by:dy
            % the current block
            blocco=dvalue(ii:ii+bx-1,jj:jj+by-1);
            % the average gray level of the current block
            max1 = max(max(blocco));
            min1 = min(min(blocco));
            
            if (max1-min1)>=thre_g
                Ymax = max(blocco(:));
                Ymin = min(blocco(:));
                Th = (3*Ymax+Ymin)/4;
                Tl = (Ymax+3*Ymin)/4;
                p=0;
                q=0;
                for i=1:bx
                    for j=1:by
                       if blocco(i,j)>=Th
                           p=p+1;
                       elseif blocco(i,j)<=Tl
                           q=q+1;
                       end
                    end    
                end    

                sum1 = 0;
                sum2 = 0;
                for i = 1:bx
                    for j = 1:by
                        if blocco(i,j)>=Th
                            sum1 = sum1+(blocco(i,j))
                        elseif blocco(i,j)<=Tl
                            sum2 = sum2+(blocco(i,j))
                        end
                    end
                end  
                 M=0.5*((1/p)*sum1)+((1/q)*sum2); 
                 R = ((1/p)*sum1)-((1/q)*sum2);
                 T1=M+(-1/3)*R;
                 T2=M;
                 T3=M+(1/3)*R;
                % the logical matrix correspoending to the current block
                blocco_binario=zeros(bx,by);
                for i = 1:bx
                    for j = 1:by
                        if blocco(i,j)>=T3
                            blocco_binario(i,j)=3;
                        elseif blocco(i,j)>=T2 && blocco(i,j)<T3
                            blocco_binario(i,j)=2;
                        elseif blocco(i,j)>=T1 && blocco(i,j)<=T2
                            blocco_binario(i,j)=1;
                        else
                            blocco_binario(i,j)=0;
                        end
                    end
                end
                % the average gray level of pixels whose level is GREATER than the block average gray level
                K=sum(sum(double(blocco_binario)));
                mu=sum(sum(double(blocco_binario).*blocco))/K;
                % the average gray level of pixels whose level is SMALLER than the block average gray level
                matrice(ii:ii+bx-1,jj:jj+by-1)=double(M+(double(2*blocco_binario)-3)/6*R);            
                posby=posby+1;
            
            else
                %disp(min1);
                % the average gray level of the current block
                m=mean(mean(blocco));
                % the logical matrix correspoending to the current block
                blocco_binario=(blocco>=m);
                % the number of pixel (of the current block) whose gray level is greater than the average gray level of the current block
                K=sum(sum(double(blocco_binario)));
                % the average gray level of pixels whose level is GREATER than the block average gray level
                mu=sum(sum(double(blocco_binario).*blocco))/K;
                % the average gray level of pixels whose level is SMALLER than the block average gray level
                if K==bx*by
                    ml=0;
                else
                    ml=sum(sum(double(~blocco_binario).*blocco))/(bx*by-K);
                end
                % the COMPRESSED DATA which correspond to the input image
                m_u(posbx,posby)=mu;                            %---> the m_u matrix
                m_l(posbx,posby)=ml;                            %---> the m_l matrix 
                mat_log(ii:ii+bx-1,jj:jj+by-1)=blocco_binario;  %---> the logical matrix
                % the compressed image

                matrice(ii:ii+bx-1,jj:jj+by-1)=(double(blocco_binario).*mu)+(double(~blocco_binario).*ml);            
                posby=posby+1;
            end
        end
        posbx=posbx+1;
    end
    out_rgb(:,:,2)=matrice;

    % -----------------------  BLUE component  ------------------   
    dvalue=double_a(:,:,3);
    max_b = max(max(dvalue));
    min_b = min(min(dvalue));
    thre_b = (max_b+min_b)/2;
   
    dx=size(dvalue,1);
    dy=size(dvalue,2);
    % if input image size is not a multiple of block size image is resized
    modx=mod(dx,bx);
    mody=mod(dy,by);
    dvalue=dvalue(1:dx-modx,1:dy-mody);
    % the new input image dimensions (pixels)
    dx=dx-modx;
    dy=dy-mody;
    % number of non overlapping blocks required to cover the entire input image
    nbx=size(dvalue,1)/bx;
    nby=size(dvalue,2)/by;
    
    % the output compressed image
    matrice=zeros(bx,by);
    % the compressed data
    m_u=zeros(nbx,nby);
    m_l=zeros(nbx,nby);
    mat_log=logical(zeros(bx,by));
    
    posbx=1;
    for ii=1:bx:dx
        posby=1;
        for jj=1:by:dy
            % the current block
            blocco=dvalue(ii:ii+bx-1,jj:jj+by-1);
            % the average gray level of the current block
            max1 = max(max(blocco));
            min1 = min(min(blocco));
            
            if (max1-min1)>=thre_b
                Ymax = max(blocco(:));
                Ymin = min(blocco(:));
                Th = (3*Ymax+Ymin)/4;
                Tl = (Ymax+3*Ymin)/4;
                p=0;
                q=0;
                for i=1:bx
                    for j=1:by
                       if blocco(i,j)>=Th
                           p=p+1;
                       elseif blocco(i,j)<=Tl
                           q=q+1;
                       end
                    end    
                end    

                sum1 = 0;
                sum2 = 0;
                for i = 1:bx
                    for j = 1:by
                        if blocco(i,j)>=Th
                            sum1 = sum1+(blocco(i,j))
                        elseif blocco(i,j)<=Tl
                            sum2 = sum2+(blocco(i,j))
                        end
                    end
                end  
                 M=0.5*((1/p)*sum1)+((1/q)*sum2); 
                 R = ((1/p)*sum1)-((1/q)*sum2);
                 T1=M+(-1/3)*R;
                 T2=M;
                 T3=M+(1/3)*R;
                % the logical matrix correspoending to the current block
                blocco_binario=zeros(bx,by);
                for i = 1:bx
                    for j = 1:by
                        if blocco(i,j)>=T3
                            blocco_binario(i,j)=3;
                        elseif blocco(i,j)>=T2 && blocco(i,j)<T3
                            blocco_binario(i,j)=2;
                        elseif blocco(i,j)>=T1 && blocco(i,j)<=T2
                            blocco_binario(i,j)=1;
                        else
                            blocco_binario(i,j)=0;
                        end
                    end
                end
                % the average gray level of pixels whose level is GREATER than the block average gray level
                K=sum(sum(double(blocco_binario)));
                mu=sum(sum(double(blocco_binario).*blocco))/K;
                matrice(ii:ii+bx-1,jj:jj+by-1)=double(M+(double(2*blocco_binario)-3)/6*R);            
                posby=posby+1;
            
            else
                %disp(min1);
                % the average gray level of the current block
                m=mean(mean(blocco));
                % the logical matrix correspoending to the current block
                blocco_binario=(blocco>=m);
                % the number of pixel (of the current block) whose gray level is greater than the average gray level of the current block
                K=sum(sum(double(blocco_binario)));
                % the average gray level of pixels whose level is GREATER than the block average gray level
                mu=sum(sum(double(blocco_binario).*blocco))/K;
                % the average gray level of pixels whose level is SMALLER than the block average gray level
                if K==bx*by
                    ml=0;
                else
                    ml=sum(sum(double(~blocco_binario).*blocco))/(bx*by-K);
                end
                % the COMPRESSED DATA which correspond to the input image
                m_u(posbx,posby)=mu;                            %---> the m_u matrix (see the cited reference)
                m_l(posbx,posby)=ml;                            %---> the m_l matrix 
                mat_log(ii:ii+bx-1,jj:jj+by-1)=blocco_binario;  %---> the logical matrix
                % the compressed image

                matrice(ii:ii+bx-1,jj:jj+by-1)=(double(blocco_binario).*mu)+(double(~blocco_binario).*ml);            
                posby=posby+1;
            end
        end
        posbx=posbx+1;
    end
    out_rgb(:,:,3)=matrice;

    %-----------------------------------------------------------
    if isa(a,'uint8')
        out=uint8(out_rgb);
        figure('Name','Compressed image');
        imshow(out);
        imwrite(out, outfile);
        return
    end
    
    if isa(a,'uint16')
        out=uint16(out_rgb);
        figure('Name','Compressed image');
        imshow(out);
        imwrite(out, outfile);
        return
    end
    
    if isa(a,'double')
        out=(out_rgb);
        figure('Name','Compressed image');
        imshow(out);
        imwrite(out, outfile);
        return
    end 
end
