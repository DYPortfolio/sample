function cmirrorsCxR2
clear

niterations = 1000;
radius1=1;

radius2 = 1.9;
center1 = 0;
center2 = [1.1 ,0];
% xys = [2 2];
xys = [3 .01];
histres = 400;
histarea = 9;
sc = (histres/2)/histarea;

rthetahist = zeros(histres,histres);
xyhist =  zeros(histres,histres);
xyhist2 =  zeros(histres,histres);
xyhisttmp = zeros(histres,histres);
xyhisttmp2 = zeros(histres,histres);

fig2 = figure;
ax2 = axes(fig2);
imagescptr = imagesc(ax2,xyhist);
set(ax2,'YDir', 'normal');
set(ax2,'YTick', []);
set(ax2,'YTickLabel', []);
set(ax2,'XTick', []);
set(ax2,'XTickLabel',[]);
set(ax2,'SortMethod','depth');
imagescptr = imagesc(ax2,'CData',log(xyhist));

C2 = [];
hitcoodrstmp = [];
circoverlayplot = [];
circoverlayplot0 = [];
xysplot=[];
rthetahistarrayindex=1;


ctmp = cos(0:pi/20:pi*2);
stmp = sin(0:pi/20:pi*2);

Rctrarray = 1./([0:.01:4]+.01);
Cctrarray = 1./([0:.01:4+.01]+.1);
 xctrarray = [-3 :1.5: 3]+.2;
 yctrarray = [-3 :1: 3]+.1;

gradientmap1 = xctrarray-min(xctrarray(:));
gradientmap1 = gradientmap1./max(gradientmap1(:));
gradientmap2 = yctrarray-min(yctrarray(:));
gradientmap2 = gradientmap2./max(gradientmap2(:));

histscellsize = [length(Rctrarray) length(Cctrarray)];
loghistscell = cell(histscellsize(1),histscellsize(2));
Chistscell =  cell(histscellsize(1),histscellsize(2));
Rhistscell = cell(histscellsize(1),histscellsize(2));

rmap = zeros(histres,histres);
bmap = zeros(histres,histres);
rmap2 = zeros(histres,histres);
bmap2 = zeros(histres,histres);

trans=@(c,r,x) (x-c)/norm((x-c)) .* (2 * r -sum((x-c).^2).^.5) ;

%trans=@(c,r,x) (c  +  ((   (r  )   ^2)/sum((   (x   - c ).^2 )))  *  ( x - c));

%trans=@(c,r,x) (c  +  ((   (r )   ^2)/sum((   (x  -c).^2 )))  *  ( x -c));

tic
for Rctrind = 1:length(Rctrarray)
    Rctr = Rctrarray(Rctrind);
    for Cctrind = 1:length(Cctrarray)
        Cctr = Cctrarray(Cctrind);
        %
        center2 =[Cctr, 0];
        radius2 = Rctr;
        %         xys = [xctr yctr];
        
        
        rmap = 0;
        bmap = 0;
        %
        
        if exist('circoverlayplot0')
            delete(circoverlayplot0)
        end
        
        if exist('circoverlayplot')
            delete(circoverlayplot)
        end
        hold on
        circoverlayplot0 = plot(ax2,ctmp*(histres/2)/histarea + histres/2 ,stmp*(histres/2)/histarea + histres/2,'w');
        circoverlayplot = plot(ax2,ctmp*radius2*(histres/2)/histarea + histres/2 + sc*center2(1),stmp*radius2*(histres/2)/histarea + histres/2 + sc*center2(2),'w');
        for yctrind = 1:length(yctrarray)
            yctr = yctrarray(yctrind);
            for xctrind = 1:length(xctrarray)
                xctr = xctrarray(xctrind);
                
                %
                xys = [xctr yctr];
                xys2=[];
                %
                for ictr = 1:niterations
                    
                    % radius1=1 center1 = 1;
                    c = center1;
                    R=radius1;
                    
                    xytmp = trans(c,R, xys(ictr,:) );
                    

%trans=@(c,r,x) (c  +  ((r^2)/sum((   x  -c).^2))  *  ( x -c));


%                     xytmp = center1 + ((radius1^2)/sum((xys(ictr,:)-center1).^2))*(xys(ictr,:)-center1);
                    
                    xytmp(isnan(xytmp))=0;
                    xys2(end+1,:) = xytmp;
                    
                    % %
                    %               plot(ax2,xytmp(1,1)*(histres/2)/histarea + histres/2,...
                    %  xytmp(1,2)*(histres/2)/histarea + histres/2,'.r');
                    % %
                    
                    xyhistxy = fXYtohistT(xytmp);
                    xyhistxy(xyhistxy<=0)=1;%add to function above
                    try
                       % xyhisttmp2(xyhistxy(1),xyhistxy(2)) =  xyhisttmp2(xyhistxy(1),xyhistxy(2)) +1;
                    catch me
                        disp('1!')
                        
                      %  return
                    end
                    
                    c = center2;
                    R=radius2;


                 %   xytmp= c  +  ((R^2)/sum((   xys(ictr,:)   -c).^2))  *  ( xys(ictr,:) -c);
                 %   xytmp2= c  +  ((R^2)/sum((   xytmp        -c).^2))  *  ( xytmp       -c);

                    xytmp2 = trans(c,R,xytmp);

%                     xytmp2 = center2 + ((radius2^2)/sum((xytmp-center2).^2))*(xytmp-center2);
                    
                    
                    xytmp2(isnan(xytmp2))=0;
                    xys(end+1,:) = xytmp2;
                    %
                    %                plot(ax2,xytmp2(1,1)*(histres/2)/histarea + histres/2,...
                    %  xytmp2(1,2)*(histres/2)/histarea + histres/2,'.g');
                    %
                    xyhistxy = fXYtohistT(xytmp2);
                    xyhistxy(xyhistxy<=0)=1;%add to function above
                    
                    try
                        xyhisttmp(xyhistxy(1),xyhistxy(2)) =  xyhisttmp(xyhistxy(1),xyhistxy(2)) +1;
                    catch me
                        disp('2!');
                        % return
                    end
                    
                    
                end
                %
                %
                rmap = rmap + (xyhisttmp+xyhisttmp2) * gradientmap1(xctrind);
                bmap = bmap + (xyhisttmp+xyhisttmp2) * (gradientmap2(yctrind));
                
                rmap2 = rmap2 + xyhisttmp2 * gradientmap1(xctrind);
                bmap2 = bmap2 + xyhisttmp2 * (gradientmap2(yctrind));
                
                xyhist = xyhist + xyhisttmp;
                xyhist2 = xyhist2 + xyhisttmp2;
                xyhisttmp = zeros(histres,histres);
                xyhisttmp2 = zeros(histres,histres);
                %
                
                % set(imagescptr,'CData',xyhist'~=0)
            end%xstarting point in (2) color(s)
            
        end
        
        rbnorm = ((rmap.^2 + bmap.^2).^.5);
        rmapnormd=(rmap)./rbnorm;
        bmapnormd=(bmap)./rbnorm;
        
        rgbmap = ones(histres,histres,3);
        
        
        rgbmap(:,:,1)=rmapnormd';
        rgbmap(:,:,3)=bmapnormd';
        
        % set(imagescptr,'CData',log(xyhist+xyhist2)');
        %add circles to hist images
        %xyhist
        
        %could be generated given index/counter
        Rhistscell{Rctrind,Cctrind} = radius2;
        Chistscell{Rctrind,Cctrind} = center2;
        loghistscell{Rctrind,Cctrind} = log(xyhist + xyhist2)';
        
        rgbmapcell{Rctrind,Cctrind} = rgbmap;
        xyhist =  zeros(histres,histres);
        xyhist2 =  zeros(histres,histres);
    end
    disp([num2str(Rctrind) '/' num2str(length(Rctrarray))]);
    
end
toc


%todo put to be frawn into cdata
set (fig2, 'WindowButtonMotionFcn', @mouseMove);

    function f_out = fXYtohistT(f_in)
      
        
        f_out = ...        
        (ceil(...
    histres.*...
    ((((f_in.*(abs(f_in)<=histarea) + histarea.*(abs(f_in)>histarea))...
    ./histarea) + 1)./2)));

        
    end

    function moseClick(object, eventdata)
    end

    function mouseMove(object, eventdata)
        if object == gcf
            try
                if isempty(C2)
                    C2 = get (gca, 'CurrentPoint');
                    C2 = C2(:,[2 1 3:end]);
                end
                if isempty(hitcoodrstmp)
                    hitcoodrstmp = [-1 -1];
                end
                
                C1 = get (gca, 'CurrentPoint');
                C1 = C1(:,[2 1 3:end]);
                
                if ~isequal(C1,C2)%mouse moved to new point
                    
                    C2 = C1;
                    if all(C2(1,1:2)>=0) && all(C2(1,1:2)< histres )%within picture %refer to what instead of Children
                        
                        hitcoodrs = floor((C2(1,1:2)/(histres-1)).*(histscellsize-1)) + 1;
                        
                        if ~isequal(hitcoodrstmp,hitcoodrs)
                            
                            
                            %set(imagescptr,'CData',loghistscell{hitcoodrs(1),hitcoodrs(2)})
                            set(imagescptr,'CData',1-rgbmapcell{hitcoodrs(1),hitcoodrs(2)})
                            
                            radiusval = Rhistscell{hitcoodrs(1),hitcoodrs(2)};
                            centerval = Chistscell{hitcoodrs(1),hitcoodrs(2)};
                            
                            if exist('circoverlayplot')
                                delete(circoverlayplot);
                            end
                            
                            if exist('circoverlayplot0')
                                delete(circoverlayplot0);
                            end
                            if exist('xysplot')
                                delete(xysplot);
                            end
                            
                            hold on% %1st radius always 1
                            circoverlayplot0 = plot(ax2,ctmp*(histres/2)/histarea + histres/2 ,stmp*(histres/2)/histarea + histres/2,'w');
                            circoverlayplot = plot(ax2,ctmp*radiusval*(histres/2)/histarea + histres/2 + sc*centerval(1),stmp*radiusval*(histres/2)/histarea + histres/2 + sc*centerval(2),'w');
                         
                            hold off
                            
                            disp([radiusval centerval]);
                            hitcoodrstmp = hitcoodrs;
                            
                        end
                    end
                end
            catch me
                disp(me)
                set (object, 'WindowButtonMotionFcn', []);
                return
            end
        end
    end
end