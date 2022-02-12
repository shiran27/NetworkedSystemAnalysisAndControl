% some praphics related to the system



timeSeries = 0:timeResolution:periodT;
figure()
for i = 1:1:numOfSubsystems
            
    ax = subplot(2,ceil(numOfSubsystems/2),i);
    ax.YColorMode = 'manual';
    
    yyaxis left
    plot(timeSeries, data(:,i,1),'-r','DisplayName','$\Vert x_i(t) \Vert$')
    hold on 
    plot(timeSeries, data(:,i,4),'-g','DisplayName','$\Vert y_i(t) \Vert$')
    ylabel('$\Vert x_i(t) \Vert$ and $\Vert y_i(t) \Vert$','Interpreter','Latex')
    title(['Subsystem - ',num2str(i)])
    ax.YColor = 'r';

    yyaxis right
    plot(timeSeries, data(:,i,2),'-b','DisplayName','$\Vert u_i(t) \Vert$')
    plot(timeSeries, data(:,i,3),'-','Color',[0,0,0,0.1],'DisplayName','$\Vert w_i(t) \Vert$')
    ylabel('$\Vert u_i(t) \Vert$ and $\Vert w_i(t) \Vert$','Interpreter','Latex')

    legend('Interpreter','Latex','Location','best')
    xlabel('Time - $t$','Interpreter','Latex')
    title(['Subsystem - ',num2str(i)])
    ax.YColor = 'b';
    
    
    
    grid on

      
        
end
















if videoMode
    % create the video writer with 1 fps
    writerObj = VideoWriter('myVideo.avi');
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(frameArray)
        % convert the image to a frame
        frame = frameArray(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end