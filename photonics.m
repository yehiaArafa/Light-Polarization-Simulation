function varargout = photonics(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @photonics_OpeningFcn, ...
                   'gui_OutputFcn',  @photonics_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end

%You have to choose an element that have a \theta degree 

% --- Executes just before photonics is made visible.
function photonics_OpeningFcn(hObject, eventdata, handles, varargin)

handles.outputText = hObject;

guidata(hObject, handles);

% UIWAIT makes photonics wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.wavePlateGroup,'visible','off');
set(handles.customGroup,'visible','off');
set(handles.linearPolarizerGroup,'visible','off');
set(handles.LQGroup,'visible','off');
set(handles.axes1,'visible','off');
set(handles.jmGroup,'visible','off');
set(handles.axes2,'visible','off');
set(handles. mustHaveTheta,'visible','off');
set(handles.qwp, 'Value', 0);
set(handles.hwp, 'Value', 0);
set(handles.wavePlateGammaVal, 'Enable', 'off');
set(handles.addElementButton, 'Enable', 'off');
set(handles.runButton, 'Enable', 'off');
set(handles.intensityRadioButton , 'Enable', 'off');
set(handles.jmRadioButton, 'Enable', 'off');
set(handles.polarizationEllipseRadioButton, 'Enable', 'off');

clear global output;
clear global selection;
clear global elements;

end


% --- Outputs from this function are returned to the command line.
function varargout = photonics_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.outputText;

end


% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
% hObject    handle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.wavePlateGroup,'visible','off');
set(handles.customGroup,'visible','off');
set(handles.linearPolarizerGroup,'visible','off');
set(handles.LQGroup,'visible','off');
set(handles.addElementButton, 'Enable', 'off');
set(handles.linearPolarizerCheckBox,'Value', 0);
set(handles.generalWavePlateCheckBox,'Value', 0);
set(handles.liquidCrystalCheckBox,'Value', 0);
set(handles.customCheckBox,'Value', 0);
set(handles.mustHaveTheta,'visible','off');

if handles.horizontalPolarization.Value == 1
    polarization =[1;0];
elseif handles.verticalPolarization.Value == 1
    polarization =[0;1];
elseif handles.rightCircularPolarization.Value == 1
    polarization = 1/sqrt(2).*[1;j];
elseif handles.leftCircularPolarization.Value == 1
    polarization = 1/sqrt(2).*[1;-j];
elseif handles.linearPolarizationTheta.Value == 1 
    thetaVal =  degtorad(str2num(get(handles.LinearPolarizationThetaInput, 'String')));
    polarization = [cos(thetaVal);sin(thetaVal)];
end


global elements;


if handles.intensityRadioButton.Value == 1
  
  VariableThetaElement = elements(get(handles.outputListBox, 'Value'));
  
  if VariableThetaElement == "QWP" | VariableThetaElement == "HWP" | VariableThetaElement == "GWP" | VariableThetaElement == "LP"
   
    axes(handles.axes2);
    plot(0:0);
    set(handles.axes2,'visible','off'); 
    set(handles.jmGroup,'visible','off'); 
    set(handles.axes1,'visible','off');
    flipedElements = fliplr(elements);
   
    variableTheta=zeros(1,201);
    
    for i = 1: 201

    clear jonesMatrix;
    jonesMatrix = [1 0 ; 0 1];
    variableTheta(i)=(i-1)/200*360;

    for elm = flipedElements
 
        switch(elm)

            case "HP"
                jonesMatrix = jonesMatrix * [1 0; 0 0];

            case "VP"
                jonesMatrix = jonesMatrix * [0 0; 0 1];

            case "LP"
                if VariableThetaElement == "LP"
                    jonesMatrix = jonesMatrix * linpolarizer(variableTheta(i));
                else
                    LPTheta = str2num(get(handles.linearPolarizerTheta, 'String'));
                    jonesMatrix = jonesMatrix * linpolarizer(LPTheta);
                end
                
            case "QWP"
                WPGamma = pi/2;
                if VariableThetaElement == "QWP"
                    jonesMatrix =  jonesMatrix * waveplate(WPGamma,variableTheta(i));
                else
                    LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                    jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
                end
                
            case "HWP"
                WPGamma = pi;
                if VariableThetaElement == "HWP"
                    jonesMatrix =  jonesMatrix * waveplate(WPGamma,variableTheta(i));
                else
                    LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                    jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
                end

            case "GWP"
                WPGamma = str2num(get(handles.wavePlateGammaVal, 'String'));
                if VariableThetaElement == "GWP"
                    jonesMatrix = jonesMatrix * waveplate(WPGamma,variableTheta(i));
                else
                    LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                    jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
                 end
 
            case "LQ"
                vo = 1.3;
                vc = 1.25;
                no = 1.5;
                ne = 1.63;
              
                lqVolt = str2num(get(handles.lqVolt, 'string'));
                if lqVolt < vc || lqVolt == vc
                    lqTheta = 0;
                else
                    lqTheta = degtorad(180/2 - 2.*atand(exp(-((lqVolt-vc)./vo ))));
                end
                    ntheta = sqrt(1./( ((cos(lqTheta)^2)./ne^2) + ((sin(lqTheta)^2)./no^2) ));
                    lqGamma = (2*pi*(ntheta-no))./(ne-no);
                    jonesMatrix = jonesMatrix * waveplate(lqGamma,45);
               
            case "C"
                custom1 = str2num(get(handles.custom1, 'string'));
                custom2 = str2num(get(handles.custom2, 'string'));
                custom3 = str2num(get(handles.custom3, 'string'));
                custom4 = str2num(get(handles.custom4, 'string'));
                jonesMatrix = jonesMatrix * [custom1 custom2; custom3 custom4];      
        end
  
    end

    jonesVector = jonesMatrix * polarization;
    intensity(i) = abs(jonesVector(1)).^2+abs(jonesVector(2)).^2;

    end

    axes(handles.axes1);
    plot(variableTheta,intensity);
    axis([0 360 0 1]); set(gca,'XTick',[0 45 90 135 180 225 270 315 360]);
    xlabel('\theta deg'); ylabel('Intensity');
    switch(VariableThetaElement)
        case "LP"
             title('Rotating Linear Polarizer');
        case "QWP"
            title('Rotating Quarter Wave Plate');
        case "HWP"
            title('Rotating Half Wave Plate');
        case "GWP"
            title('Rotating Wave Plate');
    end
   
  else
    axes(handles.axes1);
    plot(0:0);
    set(handles.axes1,'visible','off');
    set(handles. mustHaveTheta,'visible','on');
  end

elseif handles.jmRadioButton.Value == 1
    
    axes(handles.axes1);
    plot(0:0);
    set(handles.axes1,'visible','off');
    
    axes(handles.axes2);
    plot(0:0);
    set(handles.axes2,'visible','off');
    
    set(handles.jmGroup,'visible','on');
    flipedElements2 = fliplr(elements);
    clear jonesMatrix;
    jonesMatrix = [1 0 ; 0 1];
    
    for elm = flipedElements2
        switch(elm)

            case "HP"
                jonesMatrix = jonesMatrix * [1 0; 0 0];

            case "VP"
                jonesMatrix = jonesMatrix * [0 0; 0 1];

            case "LP"
                LPTheta = str2num(get(handles.linearPolarizerTheta, 'String'));
                jonesMatrix = jonesMatrix * linpolarizer(LPTheta);
             
            case "QWP"
                WPGamma = pi/2;
                LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
                
            case "HWP"
                WPGamma = pi;
                LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
                
            case "GWP"
                WPGamma = str2num(get(handles.wavePlateGammaVal, 'String'));
                LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
             
            case "LQ"
                vo = 1.3;
                vc = 1.25;
                no = 1.5;
                ne = 1.63;
              
                lqVolt = str2num(get(handles.lqVolt, 'string'));
                if lqVolt < vc || lqVolt == vc
                    lqTheta = 0;
                else
                    lqTheta =  degtorad(180/2 - 2.*atand(exp(-((lqVolt-vc)./vo ))));

                end
                    ntheta = sqrt(1./( ((cos(lqTheta)^2)./ne^2) + ((sin(lqTheta)^2)./no^2) ));
                    lqGamma = (2*pi*(ntheta-no))./(ne-no);
                    jonesMatrix = jonesMatrix * waveplate(lqGamma,45);
                 
            case "C"
                custom1 = str2num(get(handles.custom1, 'string'));
                custom2 = str2num(get(handles.custom2, 'string'));
                custom3 = str2num(get(handles.custom3, 'string'));
                custom4 = str2num(get(handles.custom4, 'string'));
                jonesMatrix = jonesMatrix * [custom1 custom2; custom3 custom4];      
        end
  
    end
    
    if (real(jonesMatrix(1,1)) == 0 && imag(jonesMatrix(1,1)) == 0)
        jm1 = num2str(0);
    elseif(real(jonesMatrix(1,1)) == 0 && imag(jonesMatrix(1,1)) ~= 0)
        jm1 = num2str(imag(jonesMatrix(1,1)));
        jm1 = strcat(jm1, "j");
    elseif (imag(jonesMatrix(1,1)) == 0 && real(jonesMatrix(1,1)) ~= 0)
        jm1 = num2str(real(jonesMatrix(1,1)));   
    else
        jm1 = num2str(real(jonesMatrix(1,1)));
        jm1 = strcat(jm1, " +");
        jm1 = strcat(jm1, " (");
        jm1 = strcat(jm1, num2str(imag(jonesMatrix(1,1))));
        jm1 = strcat(jm1, "j");
        jm1 = strcat(jm1, " )");
    end
    
    if (real(jonesMatrix(1,2)) == 0 && imag(jonesMatrix(1,2)) == 0)
        jm2 = num2str(0);
    elseif(real(jonesMatrix(1,2)) == 0 && imag(jonesMatrix(1,2)) ~= 0)
        jm2 = num2str(imag(jonesMatrix(1,2)));
        jm2 = strcat(jm2, "j");
    elseif (imag(jonesMatrix(1,2)) == 0 && real(jonesMatrix(1,2)) ~= 0)
        jm2 = num2str(real(jonesMatrix(1,2)));    
    else
        jm2 = num2str(real(jonesMatrix(1,2)));
        jm2 = strcat(jm2, " +");
        jm2 = strcat(jm2, " (");
        jm2 = strcat(jm2, num2str(imag(jonesMatrix(1,2))));
        jm2 = strcat(jm2, "j");
        jm2 = strcat(jm2, " )");
    end
    
    if (real(jonesMatrix(2,1)) == 0 && imag(jonesMatrix(2,1)) == 0)
        jm3 = num2str(0);
    elseif(real(jonesMatrix(2,1)) == 0 && imag(jonesMatrix(2,1)) ~= 0)
        jm3 = num2str(imag(jonesMatrix(2,1)));
        jm3 = strcat(jm3, "j");
    elseif (imag(jonesMatrix(2,1)) == 0 && real(jonesMatrix(2,1)) ~= 0)
        jm3 = num2str(real(jonesMatrix(1,1)));
    else
        jm3 = num2str(real(jonesMatrix(2,1)));
        jm3 = strcat(jm3, " +");
        jm3 = strcat(jm3, " (");
        jm3 = strcat(jm3, num2str(imag(jonesMatrix(2,1))));
        jm3 = strcat(jm3, "j");
        jm3 = strcat(jm3, " )");
    end
    
    if (real(jonesMatrix(2,2)) == 0 && imag(jonesMatrix(2,2)) == 0)
        jm4 = num2str(0);
    elseif(real(jonesMatrix(2,2)) == 0  && imag(jonesMatrix(2,2)) ~= 0)
        jm4 = num2str(imag(jonesMatrix(2,2)));
        jm4 = strcat(jm4, "j");
    elseif (imag(jonesMatrix(2,2)) == 0 && real(jonesMatrix(2,2)) ~= 0)
        jm4 = num2str(real(jonesMatrix(2,2)));
    else
        jm4 = num2str(real(jonesMatrix(2,2)));
        jm4 = strcat(jm4, " +");
        jm4 = strcat(jm4, " (");
        jm4 = strcat(jm4, num2str(imag(jonesMatrix(2,2))));
        jm4 = strcat(jm4, "j");
        jm4 = strcat(jm4, " )");
    end
  
   
    set(handles.jm1,'string',jm1);
    set(handles.jm2,'string',jm2);
    set(handles.jm3,'string',jm3);
    set(handles.jm4,'string',jm4 );
    
    
    
elseif handles.polarizationEllipseRadioButton.Value ==1
    axes(handles.axes1);
    plot(0:0);
    set(handles.axes1,'visible','off');
    set(handles.jmGroup,'visible','off');
    set(handles.axes2,'visible','on');
    flipedElements2 = fliplr(elements);
    clear jonesMatrix;
    jonesMatrix = [1 0 ; 0 1];
    
    for elm = flipedElements2
        switch(elm)

            case "HP"
                jonesMatrix = jonesMatrix * [1 0; 0 0];

            case "VP"
                jonesMatrix = jonesMatrix * [0 0; 0 1];

            case "LP"
                LPTheta = str2num(get(handles.linearPolarizerTheta, 'String'));
                jonesMatrix = jonesMatrix * linpolarizer(LPTheta);
             
            case "QWP"
                WPGamma = pi/2;
                LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
                
            case "HWP"
                WPGamma = pi;
                LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
                
            case "GWP"
                WPGamma = str2num(get(handles.wavePlateGammaVal, 'String'));
                LPTheta = str2num(get(handles.wavePlateThetaVal, 'String'));
                jonesMatrix = jonesMatrix * waveplate(WPGamma,LPTheta);
             
            case "LQ"
                vo = 1.3;
                vc = 1.25;
                no = 1.5;
                ne = 1.63;
              
                lqVolt = str2num(get(handles.lqVolt, 'string'));
                if lqVolt < vc || lqVolt == vc
                    lqTheta = 0;
                else
                    lqTheta = degtorad(180/2 - 2.*atand(exp(-((lqVolt-vc)./vo ))));
                end
                    ntheta = sqrt(1./( ((cos(lqTheta)^2)./ne^2) + ((sin(lqTheta)^2)./no^2) ));
                    lqGamma = (2*pi*(ntheta-no))./(ne-no);
                    jonesMatrix = jonesMatrix * waveplate(lqGamma,45);
         
            case "C"
                custom1 = str2num(get(handles.custom1, 'string'));
                custom2 = str2num(get(handles.custom2, 'string'));
                custom3 = str2num(get(handles.custom3, 'string'));
                custom4 = str2num(get(handles.custom4, 'string'));
                jonesMatrix = jonesMatrix * [custom1 custom2; custom3 custom4];      
        end
  
    end
 
    jonesVector = jonesMatrix * polarization;
    
    ax=abs(jonesVector(1,1)); ay=abs(jonesVector(2,1));
    px=atan2( imag(jonesVector(1,1)), real(jonesVector(1,1)));
    py=atan2( imag(jonesVector(2,1)), real(jonesVector(2,1)));
    am=max(ax,ay);
    % angles
    R=ay/ax; p=py-px;
    psi=1/2*atan(2*R/(1-R^2)*cos(p)); %calculate major axis rotation
    chi=1/2*asin(2*R/(1+R^2)*sin(p)); %calculate ellipticity angle
    % vectors for plot
    t=linspace(0,1,100);
    ex=ax*cos(2*pi*t+px); %x field component
    ey=ay*cos(2*pi*t+py); %y field component
    axes(handles.axes2)
    plot(ex,ey); axis([-am am -am am]); axis square; xlabel('x');
    ylabel('y');
    title('Polarization Ellipse');
    grid on;

    
    
end



end


function[T]=linpolarizer(theta)
t=theta*pi/180; 
R=@(theta)[cos(theta) sin(theta); -sin(theta) cos(theta)]; 
T=R(-t)*[1 0;0 0].'*R(t); 
end


function[T] = waveplate(gamma,theta)
t=theta*pi/180;
R=@(theta)[cos(theta) sin(theta); -sin(theta) cos(theta)];
T=R(-t)*[1 0;0 exp(-1i*gamma)].'*R(t);
end


% --- Executes on button press in addElementButton.
function addElementButton_Callback(hObject, eventdata, handles)

global selection;
global output;
global elements;
set(handles.runButton, 'Enable', 'on');
set(handles.intensityRadioButton , 'Enable', 'on');
set(handles.jmRadioButton, 'Enable', 'on');
set(handles.polarizationEllipseRadioButton, 'Enable', 'on');

switch(selection)
    case "LP"
        if handles.horizontalPolarizerRadioButton.Value == 1
         temp = "Horizontal Polarizer";
         elements = [elements, "HP"];
       elseif handles.verticalPolarizerRadioButton.Value == 1
         temp = "Vertical Polarizer";
         elements = [elements, "VP"];
       elseif handles.linearPolarizerRadioButton.Value == 1
         LPTheta = get(handles.linearPolarizerTheta, 'String');
         temp = strcat("Linear Polarizer @(",LPTheta,")");
         elements = [elements, "LP"];
        end
       output = [output, temp];
       
       
    case "WP"
       if handles.qwp.Value == 1
         gamma = "90";
         elements = [elements, "QWP"];
       elseif handles.hwp.Value == 1
         gamma = "180";
         elements = [elements, "HWP"];
       else
          gamma = get(handles.wavePlateGammaVal, 'String');
          elements = [elements, "GWP"];
       end
       theta = get(handles.wavePlateThetaVal, 'String');
       temp = strcat("Wave Plate @(",gamma, ", ",theta,")");
       output = [output, temp];
     
    case "LQ"
         elements = [elements, "LQ"];
         volt = get(handles.lqVolt, 'string');
         temp = strcat("Liquid Crystal with ", volt, " volt");
         output = [output, temp];
       
    case "C"
       custom1 = get(handles.custom1, 'string');
       custom2 = get(handles.custom2, 'string');
       custom3 = get(handles.custom3, 'string');
       custom4 = get(handles.custom4, 'string');
       temp = strcat("Custom[", custom1, " ", custom2, " ; ", custom3, " ", custom4, "]");
       output = [output, temp];
       elements = [elements, "C"];      
end

set(handles.outputListBox,'string',output);

end


% --- Executes on slider movement.
function wavePlateThetaSlider_Callback(hObject, eventdata, handles)

sliderValue = get(handles.wavePlateThetaSlider,'Value');
set(handles.wavePlateThetaVal,'String',num2str(sliderValue));

end


% --- Executes during object creation, after setting all properties.
function wavePlateThetaSlider_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'min', 0);
set(hObject, 'max', 360);
set(hObject, 'SliderStep', [1/720, 0.1]);

end


% --- Executes on slider movement.
function wavePlateGammaSlider_Callback(hObject, eventdata, handles)

set(handles.qwp, 'Value', 0);
set(handles.hwp, 'Value', 0);

set(handles.wavePlateGammaVal, 'Enable', 'on');
sliderValue = get(handles.wavePlateGammaSlider,'Value');
set(handles.wavePlateGammaVal,'String',num2str(sliderValue));

end

% --- Executes during object creation, after setting all properties.
function wavePlateGammaSlider_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'min', 0);
set(hObject, 'max', 180);
set(hObject, 'SliderStep', [1/180, 0.1]);

end


% --- Executes on button press in linearPolarizerCheckBox.
function linearPolarizerCheckBox_Callback(hObject, eventdata, handles)

if(get(hObject,'Value') == 1)
    
set(handles.addElementButton, 'Enable', 'on');
set(handles.generalWavePlateCheckBox,'Value', 0);
set(handles.liquidCrystalCheckBox,'Value', 0);
set(handles.customCheckBox,'Value', 0);

set(handles.wavePlateGroup,'visible','off');
set(handles.customGroup,'visible','off');
set(handles.LQGroup,'visible','off');
set(handles.linearPolarizerGroup,'visible','on');
   
global selection;
selection = "LP";

else
    
set(handles.addElementButton, 'Enable', 'off');
set(handles.linearPolarizerGroup,'visible','off');

 
end

end


% --- Executes on button press in generalWavePlateCheckBox.
function generalWavePlateCheckBox_Callback(hObject, eventdata, handles)

if(get(hObject,'Value') == 1)
      
set(handles.addElementButton, 'Enable', 'on');
set(handles.linearPolarizerCheckBox,'Value', 0);
set(handles.liquidCrystalCheckBox,'Value', 0);
set(handles.customCheckBox,'Value', 0);

set(handles.customGroup,'visible','off');
set(handles.linearPolarizerGroup,'visible','off');
set(handles.LQGroup,'visible','off');
set(handles.wavePlateGroup,'visible','on');

global selection;
selection = "WP";

else
    
set(handles.addElementButton, 'Enable', 'off');
set(handles.wavePlateGroup,'visible','off');

end

end


% --- Executes on button press in liquidCrystalCheckBox.
function liquidCrystalCheckBox_Callback(hObject, eventdata, handles)

if(get(hObject,'Value') == 1)
    
set(handles.addElementButton, 'Enable', 'on');
set(handles.linearPolarizerCheckBox,'Value', 0);
set(handles.generalWavePlateCheckBox,'Value', 0);
set(handles.customCheckBox,'Value', 0);

set(handles.wavePlateGroup,'visible','off');
set(handles.customGroup,'visible','off');
set(handles.linearPolarizerGroup,'visible','off');
set(handles.LQGroup,'visible','on');
    
global selection;
selection = "LQ";


else
    
set(handles.addElementButton, 'Enable', 'off');
set(handles.LQGroup,'visible','off');

end
end


% --- Executes on button press in customCheckBox.
function customCheckBox_Callback(hObject, eventdata, handles)

if(get(hObject,'Value') == 1)
    
set(handles.addElementButton, 'Enable', 'on');
set(handles.linearPolarizerCheckBox,'Value', 0);
set(handles.generalWavePlateCheckBox,'Value', 0);
set(handles.liquidCrystalCheckBox,'Value', 0);

set(handles.wavePlateGroup,'visible','off');
set(handles.linearPolarizerGroup,'visible','off');
set(handles.LQGroup,'visible','off');
set(handles.customGroup,'visible','on');
    

global selection;
selection = "C";

else
    
set(handles.addElementButton, 'Enable', 'off');
set(handles.customGroup,'visible','off');

end

end


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)

axes(handles.axes1);
plot(0:0);
set(handles.axes1,'visible','off');
 
axes(handles.axes2);
plot(0:0);
set(handles.axes2,'visible','off');
 
set(handles. mustHaveTheta,'visible','off');
set(handles.wavePlateGroup,'visible','off');
set(handles.customGroup,'visible','off');
set(handles.LQGroup,'visible','off');
set(handles.linearPolarizerGroup,'visible','off');
set(handles.jmGroup,'visible','off');
set(handles. mustHaveTheta,'visible','off');
set(handles.qwp, 'Value', 0);
set(handles.hwp, 'Value', 0);
set(handles.wavePlateGammaVal, 'Enable', 'off');
set(handles.addElementButton, 'Enable', 'off');
set(handles.runButton, 'Enable', 'off');
set(handles.intensityRadioButton , 'Enable', 'off');
set(handles.jmRadioButton, 'Enable', 'off');
set(handles.polarizationEllipseRadioButton, 'Enable', 'off');

set(handles.intensityRadioButton,'Value', 1);
set(handles.jmRadioButton,'Value', 0);
set(handles.polarizationEllipseRadioButton,'Value', 0);

set(handles.linearPolarizerCheckBox,'Value', 0);
set(handles.generalWavePlateCheckBox,'Value', 0);
set(handles.liquidCrystalCheckBox,'Value', 0);
set(handles.customCheckBox,'Value', 0);

global output;
output = [];
global elements;
elements=[];

clear global output;
clear global selection;
clear global elements;



handles.outputListBox.String = [];

end


% --- Executes during object creation, after setting all properties.
function listBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function wavePlateGammaVal_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function wavePlateGammaVal_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function wavePlateThetaVal_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function wavePlateThetaVal_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function custom1_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function custom1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function custom3_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function custom3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function custom2_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function custom2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function custom4_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function custom4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function LinearPolarizationThetaInput_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function LinearPolarizationThetaInput_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function linearPolarizerTheta_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function linearPolarizerTheta_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in outputListBox.
function outputListBox_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function outputListBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function jm1_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function jm1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function jm2_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function jm2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function jm3_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function jm3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function jm4_Callback(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function jm4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function mustHaveTheta_Callback(hObject, eventdata, handles)

end


function mustHaveTheta_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function lqVolt_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function lqVolt_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function Vc_Callback(hObject, eventdata, handles)
% hObject    handle to Vc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vc as text
%        str2double(get(hObject,'String')) returns contents of Vc as a double
end

% --- Executes during object creation, after setting all properties.
function Vc_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function ne_Callback(hObject, eventdata, handles)
% hObject    handle to ne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ne as text
%        str2double(get(hObject,'String')) returns contents of ne as a double
end

% --- Executes during object creation, after setting all properties.
function ne_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function no_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function no_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Vo_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function Vo_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
