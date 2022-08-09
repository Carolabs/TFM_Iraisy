function [] = animation(SYS, STORAGE)

% Displays the system motion

figure(1)
hold on

% TextBox with the time
mTextBox = uicontrol('style','text');
set(mTextBox,'String','Time = 0.0 s');
set(mTextBox,'Units','characters')

% Orientation of cameras
set(gca,'CameraTarget',[0 0 0])
set(gca,'CameraPosition',[0 0 2])
set(gca,'CameraUpVector',[0 1 0])

% Axis management
axis equal
axissize = SYS.L * 0.50 ;
axis([-axissize+SYS.L/2 axissize+SYS.L/2 ...
    -axissize+SYS.L/2 axissize+SYS.L/2 ...
    -axissize axissize])

% World axes
%displayFrame(eye(3), 0.0, 0.0, 0.0, 3.0);
%set(gca,'CLim',[0 100])
 
% Initialize handles
handles    = [];

% Initialization
arraysize = length(STORAGE.t);  % Number of stored points

for i=1:arraysize
    
    % Delete previously existing handles
    if ishandle(handles);      delete(handles);   end
    
    % Retrieve positions and time
    q = STORAGE.pos(:,i);
    t = STORAGE.t(i);
    
    % Change time of text box
    mystring = sprintf('Time = %0.2g s',  t);
    set(mTextBox,'String',mystring);
    
    % Draw
    handles(1) = line( [0,q(3)], [0,q(4)] );
    handles(2) = line( [SYS.xB,q(1)], [SYS.yB,q(2)], 'Color','red' );
    handles(3) = line( [q(3),q(5)], [q(4),q(6)], 'Color', 'black' ); 
    
    % Pause for visualization
    pause(0.001);
end