%Version 1.00 
%**************************************************************************
%***************************** Bispectrum GUI ******************************
%**************************************************************************
%---------------------------Credits----------------------------------------
%Wavelet Transform: Dmytro Iatsenko
%----------------------------Documentation---------------------------------
%Comnig Soon



function varargout = Bispectrum(varargin)
% BISPECTRUM MATLAB code for Bispectrum.fig
%      BISPECTRUM, by itself, creates a new BISPECTRUM or raises the existing
%      singleton*.
%
%      H = BISPECTRUM returns the handle to a new BISPECTRUM or the handle to
%      the existing singleton*.
%
%      BISPECTRUM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BISPECTRUM.M with the given input arguments.
%
%      BISPECTRUM('Property','Value',...) creates a new BISPECTRUM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bispectrum_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bispectrum_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help Bispectrum

% Last Modified by GUIDE v2.5 20-Jun-2017 14:52:46
%*************************************************************************%
%                BEGIN initialization code - DO NOT EDIT                  %
%                ----------------------------------------                 %
%*************************************************************************%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bispectrum_OpeningFcn, ...
                   'gui_OutputFcn',  @Bispectrum_OutputFcn, ...
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
%*************************************************************************%
%                END initialization code - DO NOT EDIT                    %
%*************************************************************************%


function Bispectrum_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = Bispectrum_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
function plot_type_CreateFcn(hObject, eventdata, handles)
function wavlet_transform_CreateFcn(hObject, eventdata, handles)
function surrogate_type_Callback(hObject, eventdata, handles)
function surrogate_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function surrogate_count_Callback(hObject, eventdata, handles)
function surrogate_count_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function surrogate_analysis_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function signal_name_Callback(hObject, eventdata, handles)
function signal_name_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function orientation_Callback(hObject, eventdata, handles)
function orientation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sampling_freq_Callback(hObject, eventdata, handles)
function sampling_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function max_freq_Callback(hObject, eventdata, handles)
function status_Callback(hObject, eventdata, handles)
function status_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Please Import Signal');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function max_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function min_freq_Callback(hObject, eventdata, handles)
function min_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function wavelet_type_Callback(hObject, eventdata, handles)
function wavelet_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function central_freq_Callback(hObject, eventdata, handles)
function central_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function preprocess_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cutedges_Callback(hObject, eventdata, handles)
function cutedges_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sampling_rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function length_Callback(hObject, eventdata, handles)
function length_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function intervals_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xlim_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ylim_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function detrend_signal_popup_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function display_type_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function surrogate_percentile_Callback(hObject, eventdata, handles)
set(hObject,'Enable','off');
function surrogate_percentile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sampling_rate_Callback(hObject, eventdata, handles)
function bisp_ord_Callback(hObject, eventdata, handles)
function bisp_ord_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function freq_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function freq_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%--------------------------------------------------Unused Callbacks--------

function intervals_Callback(hObject, eventdata, handles)
%Marking lines on the graphs    
    intervals = csv_to_mvar(get(handles.intervals,'String'));    
    display_selected = get(handles.display_type,'Value');
    if(size(intervals)>0)
        zval = 1;
        child_handles = allchild(handles.wt_pane);
        if(display_selected == 1 || display_selected == 2)
            for i = 1:size(child_handles,1)
                if(strcmp(get(child_handles(i),'Type'),'axes'))
                    set(child_handles(i),'Ytick',intervals);
                    hold(child_handles(i),'on');
                    warning('off');

                        for j = 1:size(intervals,2)
                            xl = get(child_handles(i),'xlim');
                            x = [xl(1) xl(2)];        
                            z = ones(1,size(x,2));
                            z = z.*zval;
                            y = intervals(j)*ones(1,size(x,2));
                            plot3(child_handles(i),x,y,z,'--k');
                        end
                    
                    warning('on');
                    hold(child_handles(i),'off');
                end
            end
        elseif display_selected >=3 
            for j = 1:size(intervals,2)
                yl = get(handles.bisp_amp_axis,'ylim');
                y = [yl(1) yl(2)];        
                x = intervals(j)*ones(1,size(y,2));
                plot(handles.bisp_amp_axis,x,y,'--k')
                yl = get(handles.bisp_phase_axis,'ylim');
                y = [yl(1) yl(2)];        
                x = intervals(j)*ones(1,size(y,2));
                plot(handles.bisp_phase_axis,x,y,'--k')
            end    
        end
    end
    
    
function preprocess_Callback(hObject, eventdata, handles)
%Detrending Part Visualisation
    data = guidata(hObject);
    sig = data.sig; 
    time_axis = data.time_axis;
    L = size(sig,2);
    fs = str2double(get(handles.sampling_freq,'String'));
    fmax = str2double(get(handles.max_freq,'String'));
    fmin = str2double(get(handles.min_freq,'String'));
    
    contents = cellstr(get(handles.detrend_signal_popup,'String'));
    i = contents{get(handles.detrend_signal_popup,'Value')};
    i = str2double(i);
    
    %Detrending
    cur_sig = sig(i,:);
    cur_sig = cur_sig(:);
    X=(1:length(cur_sig))'/fs; XM=ones(length(X),4); 

    for pn=1:3 
        CX=X.^pn; 
        XM(:,pn+1)=(CX-mean(CX))/std(CX); 
    end

    w=warning('off','all'); 
    new_signal=cur_sig-XM*(pinv(XM)*cur_sig); 
    warning(w);

    %Filtering
    fx=fft(new_signal,L); % Fourier transform of a signal

    Nq=ceil((L+1)/2); 
    ff=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L; 
    ff=ff(:); % frequencies in Fourier transform

    fx(abs(ff)<=max([fmin,fs/L]) | abs(ff)>=fmax)=0; % filter signal in a chosen frequency domain
    new_signal=ifft(fx);
    %Plotting
    
    plot(handles.plot_pp,time_axis,cur_sig);
    hold(handles.plot_pp,'on');
    plot(handles.plot_pp,time_axis,new_signal,'-r');
    
    
    %legend(handles.plot_pp,'Original','Pre-Processed','Location','Best');
    xlim(handles.plot_pp,[0,size(sig,2)./fs]);
%-------------------------------------------------------------------------    

function wavlet_transform_Callback(hObject, eventdata, handles)
%Does the wavelet transform 
% Get user input from GUI
    
    set(handles.status,'String','Calculating Wavelet Transform...');
    fs = str2double(get(handles.sampling_freq,'String'));
    fmax = str2double(get(handles.max_freq,'String'));
    fmin = str2double(get(handles.min_freq,'String'));
    fc =  str2double(get(handles.central_freq,'String'));
    ord = str2double(get(handles.bisp_ord,'String'));
    items = get(handles.wavelet_type,'String');
    index_selected = get(handles.wavelet_type,'Value');
    wavelet_type_selected = items{index_selected};
    
    items = get(handles.preprocess,'String');
    index_selected = get(handles.preprocess,'Value');
    preprocess_selected = items{index_selected};
    
    items = get(handles.cutedges,'String');
    index_selected = get(handles.cutedges,'Value');
    cutedges_selected = items{index_selected};
    
    sig = handles.sig;    
    
    n = size(sig,1) ;
    handles.WT = cell(n, 1);
    
    if isnan(ord)
      errordlg('Order must be specified','Parameter Error');
      return
    end
    
%Taking only selected part of the signal
    xl = get(handles.xlim,'String');
    xl = csv_to_mvar(xl);
    xl = xl.*fs;
    xl(2) = min(xl(2),size(sig,2));
    xl(1) = max(xl(1),1);
    sig = sig(:,xl(1):xl(2));
    xl = xl./fs;
      
    set(handles.status,'String','Calculating Wavelet Transform...');
    %Calculating wavelet transform and deciding parameter form
    if(isnan(fmax)&& isnan(fmin))
        if(isnan(fc))
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
            end
        else
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
            end
        end
    elseif(isnan(fmax))
        if(isnan(fc))
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);
            end
        else
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);
            end
        end
    elseif(isnan(fmin))
        if(isnan(fc))
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);
            end
        else
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);
            end
        end
    else
        if(isnan(fc))
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);
            end
        else
            for p = 1:n
            [handles.WT{p,1},handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);
            end
        end
    end
    
    %---------------
    
    handles.amp_WT = cell(n,1);
    handles.pow_WT = cell(n,1);
    handles.pow_arr = cell(n,1);
    handles.amp_arr = cell(n,1);
    
    %write the averaging function and the max function
    
    for p = 1:n
        handles.amp_WT{p,1} = abs(handles.WT{p,1});   
        handles.pow_WT{p,1} = abs(handles.WT{p,1}).^2;
        handles.pow_arr{p,1} = nanmean(handles.pow_WT{p,1}.');%Calculating Average Power
        handles.amp_arr{p,1} = nanmean(handles.amp_WT{p,1}.');%Calculating Average Amplitude  
    end
    
    %Calculating the Bispectrum
    [handles.bxxx, handles.bppp, handles.bpxx, handles.bxpp] = bispectrum(handles.sig(1,:), handles.sig(2,:),...
        handles.WT{1,1}, handles.WT{2,1}, handles.freqarr, .01, ord, fs, fc);
    
    guidata(hObject,handles);    
    display_type_Callback(hObject, eventdata, handles);
    guidata(hObject,handles);
    
    set(handles.display_type,'Enable','on');
    set(handles.intervals,'Enable','on');
  
function display_type_Callback(hObject, eventdata, handles)
tic
set(handles.freq_1,'String','');
set(handles.freq_2,'String','');
display_selection = get(hObject,'Value');

set(handles.status,'String','Plotting Data');

sig = handles.sig;
fs = str2double(get(handles.sampling_freq,'String'));
n = size(sig,2)/1000;
xl = csv_to_mvar(get(handles.xlim,'String'));
xl = xl.*fs;
xl(2) = min(xl(2),size(sig,2));
xl(1) = max(xl(1),1);
xl = xl./fs;
time_axis = xl(1):1/fs:xl(2);

if (display_selection == 1 || display_selection == 2) && isfield(handles,'WT')
    
    clear_pane_axes(handles.wt_pane);
    %Deciding the plot type
        
        %Actual Plotting    
        %-------------------------Surf Plot------------------------------------
        
        position = [0.06 0.122 0.6 0.849];
        handles.plot3d = subplot(1,3,[1 2],'Parent',handles.wt_pane,'position',position);
        position = [.75 .122 .196 .849];
        handles.plot_pow = subplot(1,3,3,'Parent',handles.wt_pane,'position',position);
        
        if(handles.plot_type == 1)           
            handles.peak_value = max(handles.pow_WT{display_selection,1}(:))+.1;
            pcolor(handles.plot3d, time_axis(1:n:end) ,handles.freqarr, handles.pow_WT{display_selection,1}(1:end,1:n:end)); 
            zlabel(handles.plot3d,'Power','fontweight','b','fontsize',12);
            plot(handles.plot_pow ,handles.pow_arr{display_selection,1}, handles.freqarr,'-k','LineWidth',3 );
        else            
            handles.peak_value = max(handles.amp_WT{display_selection,1}(:))+.1;
            pcolor(handles.plot3d, time_axis(1:n:end) ,handles.freqarr, handles.amp_WT{display_selection,1}(1:end,1:n:end)); 
            zlabel(handles.plot3d,'Amplitude','fontweight','b','fontsize',12);
            plot(handles.plot_pow ,handles.amp_arr{display_selection,1}, handles.freqarr,'-k','LineWidth',3 );
        end 
        
        shading(handles.plot3d,'interp');        
        set(handles.plot3d,'yscale','log');
        set(handles.plot3d,'ylim',[min(handles.freqarr) max(handles.freqarr)]);%making the axes tight
        set(handles.plot3d,'xlim',[time_axis(1) time_axis(end)]);%making the axes tight
        xlabel(handles.plot3d,'Time (s)','fontweight','b','fontsize',12);
        ylabel(handles.plot3d,'Frequency (Hz)','fontweight','b','fontsize',12);
        
        %------------------------Power Plot------------------------------------      
                        
        set(handles.plot_pow,'yscale','log');        
        ylim(handles.plot_pow,[min(handles.freqarr) max(handles.freqarr)]);
        set(handles.status,'String','Done Plotting');
        if(handles.plot_type == 1)       
            xlabel(handles.plot_pow,'Average Power','fontweight','b','fontsize',12);
        else   
            xlabel(handles.plot_pow,'Average Amplitude','fontweight','b','fontsize',12);
        end
        ylabel(handles.plot_pow,'Frequency (Hz)','fontweight','b','fontsize',12);
        
elseif display_selection == 3 || display_selection == 4 || display_selection == 5 || display_selection == 6
     %Plotting bispectrum   
        clear_pane_axes(handles.wt_pane);    
        position = [0.06 0.12 0.432 0.8];
        handles.bisp = axes('Parent',handles.wt_pane,'position',position,'fontsize',0.03);            
        position = [.56 .55 .42 .42];
        handles.bisp_amp_axis = axes('Parent',handles.wt_pane,'position',position,'fontsize',0.05);
        position = [.56 .09 .42 .42];
        handles.bisp_phase_axis = axes('Parent',handles.wt_pane,'position',position,'fontsize',0.05);
        if display_selection == 3
            pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bxxx) 
        elseif display_selection == 4
            pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bppp)
        elseif display_selection == 5
            pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bxpp)
        elseif display_selection == 6
            pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bpxx)
        end
        set(handles.bisp,'fontsize',0.03, 'yscale','log','xscale','log');
        idx_first = find(sum(~isnan(handles.bxxx),1) > 0, 1 ,'first');
        idx_last = find(sum(~isnan(handles.bxxx),1) > 0, 1 , 'last');      
        xlim(handles.bisp,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
        ylim(handles.bisp,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);        
        title(handles.bisp,'Bispectrum','fontsize',14,'fontweight','normal');
        xlabel(handles.bisp,'Frequency (Hz)');
        ylabel(handles.bisp,'Frequency (Hz)');
        shading(handles.bisp,'interp');
        
elseif display_selection == 7 && isfield(handles,'WT')
    %Plotting all plots
%         clear_pane_axes(handles.wt_pane);    
%         position = [0.06 0.55 0.22 0.40];
%         handles.bispxxx = axes('Parent',handles.wt_pane,'position',position);
%         position = [.331 .55 .22 .40];
%         handles.bispxpp = axes('Parent',handles.wt_pane,'position',position);
%         position = [.487 .07 .22 .40];
%         handles.bisppxx = axes('Parent',handles.wt_pane,'position',position);
%         position = [.758 .07 .22 .40];
%         handles.bispppp = axes('Parent',handles.wt_pane,'position',position);
%         position = [.64 .55 .34 .40];
%         handles.bisp_amp_axis = axes('Parent',handles.wt_pane,'position',position);
%         
%         position = [.06 .07 .34 .40];
%         handles.bisp_phase_axis = axes('Parent',handles.wt_pane,'position',position);
%         
%         pcolor(handles.bispxxx, handles.freqarr, handles.freqarr, handles.bxxx)
%         pcolor(handles.bispppp, handles.freqarr, handles.freqarr, handles.bppp)
%         pcolor(handles.bisppxx, handles.freqarr, handles.freqarr, handles.bpxx)
%         pcolor(handles.bispxpp, handles.freqarr, handles.freqarr, handles.bxpp)
%         
%         child_handles = allchild(handles.wt_pane);
%         for i = 1:size(child_handles,1)
%             if(strcmp(get(child_handles(i),'Type'),'axes'))     
%                 xlabel(child_handles(i),'Frequency (Hz)');
%                 shading(child_handles(i),'interp');
%                 set(child_handles(i),'yscale','log');
%                 set(child_handles(i),'xscale','log');
%             end
%         end               
else 
    error('Calculate Before Plotting');
end
guidata(hObject,handles);

function calculate_bisp_Callback(hObject, eventdata, handles)
%Plotting the Biphase and Biamp
display_selection = get(handles.display_type,'Value');
    if(display_selection>=3 && display_selection<=6) 
        hold(handles.bisp_amp_axis,'on');
        hold(handles.bisp_phase_axis,'on');
        f1 = str2double(get(handles.freq_1,'String'));
        f2 = str2double(get(handles.freq_2,'String'));
        fs = str2double(get(handles.sampling_freq,'String'));
        fc =  str2double(get(handles.central_freq,'String'));
        xl = csv_to_mvar(get(handles.xlim,'String'));
        xl = xl.*fs;
        xl(2) = min(xl(2),size(handles.sig,2));
        xl(1) = max(xl(1),1);
        xl = xl./fs;
        time_axis = xl(1):1/fs:xl(2);
        %Verify Calculation from Ola
        if display_selection == 3
            [handles.biamp, handles.biphase] = biphaseWav(handles.sig(1,:), handles.WT{1,1}, handles.WT{1,1}, handles.freqarr, f1, f2, fs, fc);
        elseif display_selection == 4
            [handles.biamp, handles.biphase] = biphaseWav(handles.sig(2,:), handles.WT{2,1}, handles.WT{2,1}, handles.freqarr, f1, f2, fs, fc);
        elseif display_selection == 5
            [handles.biamp, handles.biphase] = biphaseWav(handles.sig(2,:), handles.WT{1,1}, handles.WT{2,1}, handles.freqarr, f1, f2, fs, fc);
        elseif display_selection == 6
            [handles.biamp, handles.biphase] = biphaseWav(handles.sig(1,:), handles.WT{2,1}, handles.WT{1,1}, handles.freqarr, f1, f2, fs, fc);
        end
        plot(handles.bisp_amp_axis, time_axis, handles.biamp);
        plot(handles.bisp_phase_axis, time_axis, handles.biphase);        
        ylabel(handles.bisp_amp_axis,'Biamplitdue');
        ylabel(handles.bisp_phase_axis,'Biphase');
        xlabel(handles.bisp_phase_axis,'Time (s)');        
        guidata(hObject, handles);
    end
%Marking the point 

child_handles = allchild(handles.bisp);
for i = 1:size(child_handles,1)    
    if(strcmp(get(child_handles(i),'Type'),'line'))
        xdat = get(child_handles(i),'XData');
        mark = get(child_handles(i),'Marker');
        if(length(xdat) == 1 && strcmp(mark,'*'))
            set(child_handles(i),'Marker','o');
            set(child_handles(i),'MarkerEdgeColor','y');
        end
    end
end

function freq_2_Callback(hObject, eventdata, handles)
x = str2double(get(handles.freq_1,'String'));
y = str2double(get(handles.freq_2,'String'));
child_handles = allchild(handles.bisp);
for i = 1:size(child_handles,1)    
    if(strcmp(get(child_handles(i),'Type'),'line'))
        xdat = get(child_handles(i),'XData');
        mark = get(child_handles(i),'Marker');
        if(length(xdat) == 1 && strcmp(mark,'*'))
            delete(child_handles(i))
        end
    end
end
hold(handles.bisp,'on');
plot(handles.bisp, x, y, 'r*')

function freq_1_Callback(hObject, eventdata, handles)
x = str2double(get(handles.freq_1,'String'));
y = str2double(get(handles.freq_2,'String'));
child_handles = allchild(handles.bisp);
for i = 1:size(child_handles,1)    
    if(strcmp(get(child_handles(i),'Type'),'line'))
        xdat = get(child_handles(i),'XData');
        mark = get(child_handles(i),'Marker');
        if(length(xdat) == 1 && strcmp(mark,'*'))
            delete(child_handles(i))
        end
    end
end
hold(handles.bisp,'on');
plot(handles.bisp, x, y, 'r*')

function select_freq_Callback(hObject, eventdata, handles)
[x, y] = ginput(1);
child_handles = allchild(handles.bisp);
for i = 1:size(child_handles,1)    
    if(strcmp(get(child_handles(i),'Type'),'line'))
        xdat = get(child_handles(i),'XData');
        mark = get(child_handles(i),'Marker');
        if(length(xdat) == 1 && strcmp(mark,'*'))
            delete(child_handles(i))
        end
    end
end

hold(handles.bisp,'on');
plot(handles.bisp, x, y, 'r*')
set(handles.freq_1,'String',x);
set(handles.freq_2,'String',y);

% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
%Loading data

% --------------------------------------------------------------------
function csv_read_Callback(hObject, eventdata, handles)
%Read csv file
    set(handles.status,'String','Importing Signal...');
    
    sig = read_from_csv();
    fs = str2double(get(handles.sampling_freq,'String')); 
    
    data = guidata(hObject);
    data.sig = sig;   
    time = 1:size(sig,2);
    time = time./fs;
    data.time_axis = time;
    guidata(hObject,data);    
    
    
    time = 1:size(handles.sig,2);
    time = time./fs;          
    handles.time_axis = time;
    plot(handles.time_series_1, time, handles.sig(1,:));%Plotting the time_series part afte calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(handles.sig,2)./fs]);
    ylabel(handles.time_series_1,'Signal 1');
    
    
    plot(handles.time_series_2, time, handles.sig(2,:));%Plotting the time_series part afte calculation of appropriate limits
    xlim(handles.time_series_2,[0,size( handles.sig,2)./fs]);
    xlabel(handles.time_series_2,'Time (s)');
    ylabel(handles.time_series_2,'Signal 2');  
    linkaxes([handles.time_series_1 handles.time_series_2],'x');
        
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    
    guidata(hObject,handles);
    preprocess_Callback(hObject, eventdata, handles);%plots the detrended curve
    guidata(hObject,handles);
    get(handles.time_series_2,'fontsize')
    get(handles.time_series_2,'fontunits')
    set(handles.status,'String','Select Data And Continue With Wavelet Transform'); 

% --------------------------------------------------------------------
function mat_read_Callback(hObject, eventdata, handles)
%Read mat file    
    set(handles.status,'String','Importing Signal...');
    handles.sig = read_from_mat(); 
    handles.sig = struct2cell(handles.sig);
    handles.sig = cell2mat(handles.sig);
    fs = str2double(get(handles.sampling_freq,'String')); 
 
    time = 1:size(handles.sig,2);
    time = time./fs;          
    handles.time_axis = time;
    plot(handles.time_series_1, time, handles.sig(1,:));%Plotting the time_series part afte calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(handles.sig,2)./fs]);
    ylabel(handles.time_series_1,'Signal 1');
    
    
    plot(handles.time_series_2, time, handles.sig(2,:));%Plotting the time_series part afte calculation of appropriate limits
    xlim(handles.time_series_2,[0,size( handles.sig,2)./fs]);
    xlabel(handles.time_series_2,'Time (s)');
    ylabel(handles.time_series_2,'Signal 2');  
    linkaxes([handles.time_series_1 handles.time_series_2],'x');
        
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    
    guidata(hObject,handles);
    preprocess_Callback(hObject, eventdata, handles);%plots the detrended curve
    guidata(hObject,handles);
    get(handles.time_series_2,'fontsize')
    get(handles.time_series_2,'fontunits')
    set(handles.status,'String','Select Data And Continue With Wavelet Transform'); 
    
%---------------------------Limits-----------------------------
function xlim_Callback(hObject, eventdata, handles)
%When the values of xlim are changed the graphs are updated
    xl = csv_to_mvar(get(handles.xlim,'String'));
    xlim(handles.time_series_1,xl);
    xlim(handles.time_series_2,xl);
    xlim(handles.plot_pp,xl);
    t = xl(2) - xl(1);
    set(handles.length,'String',t);

function ylim_Callback(hObject, eventdata, handles)
%When the values of ylim are changed the graphs are updated  
    yl = csv_to_mvar(get(handles.ylim,'String'));
    ylim(handles.time_series_1,yl);
    ylim(handles.time_series_2,yl);


%---------------------------Updating Value of limits Limits-----------------------------
function refresh_limits_Callback(hObject, eventdata, handles)
%Calcualtes limits of the plot    
    
    x = get(handles.time_series_1,'xlim');
    xlim(handles.plot_pp,x);
    t = x(2) - x(1);
    x = strcat(num2str(x(1)),' , ',num2str(x(2)));    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat(num2str(y(1)),' , ',num2str(y(2)));
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);
    
% ---------------------------Zoom Updating--------------------------
function zoom_in_OffCallback(hObject, eventdata, handles)
%Refreshes the limit values right after the tool is deselected
    x = get(handles.time_series_1,'xlim');
    xlim(handles.plot_pp,x);
    t = x(2) - x(1);
    x = strcat(num2str(x(1)),' , ',num2str(x(2)));    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat(num2str(y(1)),' , ',num2str(y(2)));
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);

% -----------------------------Zoom Updating--------------------------
function zoom_out_OffCallback(hObject, eventdata, handles)
%Refreshes the limit values right after the tool is deselected
    x = get(handles.time_series_1,'xlim');
    xlim(handles.plot_pp,x);
    t = x(2) - x(1);
    x = strcat(num2str(x(1)),' , ',num2str(x(2)));    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat(num2str(y(1)),' , ',num2str(y(2)));
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);
    
function plot_type_SelectionChangeFcn(hObject, eventdata, handles)
%deciding which plot
    switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
        case 'power'
            plot_type = 1;
        case 'amp'
            plot_type = 2;
    end

    data = guidata(hObject);
    data.plot_type = plot_type;
    guidata(hObject,data); 

% ----------------------------------------Saving Files---------------
function save_Callback(hObject, eventdata, handles)
%Honestly you're just here because I don't know how to get rid of you

function save_3dplot_Callback(hObject, eventdata, handles)
%Saves the 3d plot
    Fig = figure;
    copyobj(handles.plot3d, Fig);
    Fig = tightfig(Fig);

function save_power_plot_Callback(hObject, eventdata, handles)
%Saves the power plot
    Fig = figure;
    copyobj(handles.plot_pow, Fig);
    Fig = tightfig(Fig);

function save_pow_arr_Callback(hObject, eventdata, handles)
%Saves the avg power array
    [FileName,PathName] = uiputfile
    save_location = strcat(PathName,FileName)
    data = guidata(hObject);
    pow_arr = data.pow_arr;
    save(save_location,'pow_arr');

function save_wt_Callback(hObject, eventdata, handles)
%Saves the wavelet transform
    [FileName,PathName] = uiputfile
    save_location = strcat(PathName,FileName)
    WT = handles.WT;
    freqarr = handles.freqarr;
    save(save_location,'freqarr','-v7.3');
    save(save_location,'WT','-v7.3');%Sometimes the compression is faulty

function detrend_signal_popup_Callback(hObject, eventdata, handles)
%Detrends the signal plots the chosen one
    cla(handles.plot_pp,'reset');
    preprocess_Callback(hObject, eventdata, handles);

function surrogate_analysis_Callback(hObject, eventdata, handles)
%To enable and disable the Percentile box
surrogate_analysis = get(hObject,'Value');
if surrogate_analysis == 1
     set(handles.surrogate_percentile,'Enable','off');
else
     set(handles.surrogate_percentile,'Enable','on');
end

function bisp_clear_Callback(hObject, eventdata, handles)
cla(handles.bisp_amp_axis,'reset');
cla(handles.bisp_phase_axis,'reset');
set(handles.freq_1,'String','');
set(handles.freq_2,'String','');
clear_axes_points(handles.bisp);

