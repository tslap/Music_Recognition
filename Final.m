filename = 'lovesong.wav'
[clip,sr] = audioread(filename);


t = audioinfo(filename).Duration;
N = audioinfo(filename).TotalSamples;

function Transcriber(clip,t,dt,N,sr)
  %Makes everything two track
  if size(clip)(2) ~= 2
    clip = [clip,clip];
  endif
  
  
  Clips = ClipWindower(clip,t,dt,sr);
  Flux = Flux(Clips);
  ts = SpecFlux(Flux,t,dt);
  Tempo = Rhythm(Flux,t,dt)
  FluxFinder(Clips,ts,N,sr)
  
 
endfunction

function FFTs = ClipWindower(clip,t,dt,sr)
  %This will split the clip into smaller clips of length dt
  windowLength = dt*sr;
  
  
  %slices inputed clip into subclips
  Windows = [];
  Clips = [clip(1:windowLength,:)];
  for n = 1: floor(2*t/dt)-2
   Windows = cat(3, Windows, clip( (n*windowLength)/2 : ((n+2)*windowLength)/2-1 , : ));
  endfor
  
  %Creates Hanning function
  X = [0:windowLength-1];
  Hann = (1/2)*(1-cos(2*pi*X/windowLength));
  Hann = transpose(Hann);
  Hann = cat(2, Hann, Hann);
  
  %Implements hanninng funnction
  for n = 1:size(Windows)(3)
    Windows(:,:,n) = Windows(:,:,n).*Hann;
  endfor
  
  %Takes the FFT of each window
  FFTs = [fft(Windows(:,:,1))];
  for n = 2:size(Windows)(3)
    FFTs = cat(3, FFTs, fft(Windows(:,:,n)));
  endfor
  
  
  %Takes only the positive frequencies 
  FFTs = FFTs(1:size(FFTs)(1)/2,:,:);
  
  %averages two channels into one
  FFTs = (FFTs(:,1,:) + FFTs(:,2,:)) / 2;
  FFTs = squeeze(abs(FFTs));
endfunction

function Flux = Flux(FFTs)
  sz = size(FFTs);
  Flux = [];
  
  %remove values that are only 0s

  for n = 1:sz(2)-1
    
    %Stops division by 0 for silennt clips
    if sum(FFTs(:,n)) == 0
      FFTs(1,n) = 1;
    endif
    
    %normalization
    FFTs(:,n) = FFTs(:,n) / sum(FFTs(:,n));
    
    %computes Spectral Flux
    Sum =0;    
    Sum = sum(FFTs(:,n) - FFTs(:,n+1))^2;
    Flux = [Flux; Sum];
    
  endfor
  

endfunction

function ts = SpecFlux(Flux, t, dt)
  
  sz = size(Flux);
  avg = sum(Flux)/size(Flux)(1);
  
  %Finds large peaks of flux for onset detection
  [mag, ts,w] = findpeaks(Flux, "MinPeakDistance",2/(200/60*dt), "MinPeakHeight" , avg/100);
  
  ts = [ts] ;
  ts = [1;ts;size(Flux)(1)]

  
endfunction

function Tempo = Rhythm(Flux, t, dt)

  
  szF = size(Flux);

  %Autocorrelation of Flux to find frequency
  [R, lag] = xcorr(Flux);
  szR = size(R);
  figure(2)
  
  %Finds large peaks which correspond to tempo and overtones 
  [peak,spaces] = findpeaks(R,"MinPeakDistance", 2/(300/60*dt));
  szS = size(spaces)(1);
  
  %Finds the average distance between peaks to find fundemental frequency
  Kappa = 0;
  for n=1:szS-1
    Kappa = Kappa + spaces(n+1) - spaces(n);
  endfor
  avgSpace = Kappa/(szS-1);


  %Unit conversion from clips -> seconds
  Tempo= 1/(avgSpace*(dt/2))*60;
  
  %Create Comb at tempo to produce plot
  comb = zeros(1,round(avgSpace));
  comb(1) = max(Flux);
  comb = repmat(comb, 1,round(szF(1) / avgSpace));
  
  
  figure(1)
  plot(Flux, 'linewidth', 0.75)
  title('Spectral Flux over time')
  xlabel('x * 0.16s')
  ylabel('Flux')
  set(gca,'YTickLabel',[]);
  xlim([0 200])
  grid minor

  plot(lag,R, 'linewidth', 0.75)
  grid minor
  title('Autocorrelation of Flux')
  set(gca, 'YTickLabel',[])
  ylabel('Magnitude')
  xlabel('X')
  xlim([-200 200])

  
  
  figure(3)
  stem(comb, 'linewidth', 0.5)
  hold on 
  plot(Flux, 'linewidth', 0.75)
  xlim([0 200])
  box on;
  title('Spectral Flux and Predicted Downbeats')
  xlabel('x * 0.16s')
  ylabel('Flux')
  set(gca, 'YTickLabel', []);
  hold off
  
endfunction

function FluxFinder(FFTs, ts,N,sr)
  
  deltaT = [];
  Notes = [];
  
  
  %Averages all of the FFTs for each note i.e. between two onsets
  AvgSum = 0;
  for n=1:size(ts)(1)-1
    for i = ts(n):ts(n+1)
      AvgSum = AvgSum + FFTs(:,i);
    endfor
    AvgFFT = AvgSum; 
    
    [mag, index] = max(AvgFFT);
    index/0.32;
    
    
    %Sends averaged FFTs to pitch function for detection
    note = pitch(AvgFFT,N,sr);
    Notes = [Notes; note];
    
    
  AvgSum = 0;
  AvgFFT = 0;
    
  endfor
 
  Notes
endfunction

function note = pitch(FFT,N,sr)
  
  %Autocorrelation of the fourier transform to pick out overtones  
  [R , lag] = xcorr(FFT);
  
  %limits range to common frequencies to limit unnecessary computation
  szR = size(R);
  R = R(round(szR(1)/2 -1000):round(szR(1)/2 +1000));
  szLag = size(lag);
  lag = (szLag(2)/2 - 1000:szLag(2)/2 + 1000);
  
  
  %Picks out large peaks from the Autocorrelation
  [val , pitch] = findpeaks(R, "MinPeakDistance", 20, "MinPeakHeight",max(R)/10 );
  
  
  %Finds average spacing inbetween peaks 
  Kappa = 0;
  for n=1:size(pitch)-1
    Kappa = Kappa + pitch(n+1) - pitch(n);
  endfor
  avgPitch = Kappa/ (size(pitch)(1)-1) ;
  
  %Unit conversion to hz before sending the note to be named
  freq = avgPitch * (sr/2)/size(FFT)(1);
  note = NoteIdentifyer(freq);
  
  
  figure(1)
  plot(FFT)
  xlim([0 400])
  title('Fourier Transform')
  xlabel('Frequency')
  ylabel('magnitude')
  
  figure(2)
  plot(lag,R)
  title('AutoCorrelation of Fourier Transform')
  ylabel('magnitude')
  xlabel('X')
  set(gca, 'XTickLabel', [])
  xlim([6000 8200])

endfunction





function Note = NoteIdentifyer(f)
  
  
  
  %finds interval from C4
  n = 12*log2(f/261.63);
  
  %Finds the difference in octave and difference from C within said octave
  DelOctave = floor(n/12);
  DelNote = round(n - 12 * DelOctave);
  
  Octave = 4 + DelOctave;
  
  if DelNote == 0 || DelNote == 12
    Note = 'C';
    elseif DelNote == 1
      Note = 'C#/Db';
    elseif DelNote == 2
      Note = 'D';
    elseif DelNote == 3
      Note = 'D#/Eb';
    elseif DelNote == 4
      Note = 'E';
    elseif DelNote == 5
      Note = 'F';
    elseif DelNote == 6
      Note = 'F#/Gb';
    elseif DelNote == 7
      Note = 'G';
    elseif DelNote == 8
      Note = 'G#/Ab';
    elseif DelNote == 9
      Note = 'A';
    elseif DelNote == 10
      Note = 'A#/Bb';
    elseif DelNote == 11
      Note = 'B';
   end
    
  Note = strcat('(',Note,num2str(Octave),')');
  
endfunction


