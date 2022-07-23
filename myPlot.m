function myPlot(SD,valCurr,val,lowValue,highValue,n,T,arr,fig,videoObj)
        [nx,ny,nz] = size(SD);
        subplot 131
        slice(valCurr,ceil(nx/2),ceil(ny/2),ceil(nz/2));  %shading flat;
        title('Wave');
        %         caxis([lowValue highValue]);
        colorbar;
        axis equal;view(3);
        %view([0,0]);
      
        %         F = getframe(fig);
        %         writeVideo(vidObj, F);
        subplot 132
        slice(SD,ceil(nx/2),ceil(ny/2),ceil(nz/2)), axis equal; %shading flat;
        title('Model');view(3);
        %view([0,0]);
        subplot 133
        slice(val,ceil(nx/2),ceil(ny/2),ceil(nz/2));axis equal; %shading flat;
        caxis([-0.1    0.1]);
        title(['Time step :', num2str(n), ' of ', num2str(T)]);view(3);
        colorbar;
        axis equal;
        %view([0,0]);
        drawnow;
        F = getframe(fig);
        writeVideo(videoObj, F);
        
end