%% Modified DASH algorithm to correct hologram-to-hologram random phase offset
% "T" is the uncorrected transmission matrix.

function [T phscorrect signals Cp]=DASH_THG_FOV(T,Diter)



% note: this script assumes that the FOV is square. For example, the FOV
% in the original work is composed of a 26x26 holograms. 
N=sqrt(size(T,2));
%fieldsize=[N N];
k=fftfreq(N,1/N);
[kxx kyy]=meshgrid(k,k);
indd=1;
M=zeros(N,N,N*N);
for ii=1:N
    for jj=1:N
        M(:,:,indd)=fftshift(((-pi+ii*2*pi/N)*kxx+(-pi+jj*2*pi/N)*kyy));
        indd=indd+1;
    end
end



Cp=zeros(N,N);
phasesteps=[0 2*pi/5 4*pi/5 6*pi/5 8*pi/5];
f=.3;

%progress bar is another function, used to track the status/duration of the
%calcluation. 

upd = textprogressbar(size(T,2)*Diter, 'barlength', 20, ...
                         'updatestep', 1, ...
                         'startmsg', 'Correcting pupil phase forward... ',...
                         'endmsg', ' Finally!', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '+', ...
                         'emptybarsymbol', '-');


ind=1;
for mm=1:Diter
for ii=1:size(M,3)
    for jj=1:length(phasesteps)
        ESLM=exp(1i*angle(sqrt(f)*exp(1i*M(:,:,ii)+1i*phasesteps(jj))+sqrt(1-f)*exp(1i*angle(Cp))));
        % Calc Intensity
        Intensity(jj)=sum(abs(T*ESLM(:)).^2,'all');
        
    end
    %weighting function
    a(ii)=sum(sqrt(Intensity).*exp(1i.*phasesteps),'all')./length(phasesteps);
    Cp=Cp+a(ii).*exp(1i*M(:,:,ii));
    signals(ind) = mean(Intensity);
    ind=ind+1;
    upd(ind)
    %disp(ii);
end
end

figure; plot(signals)
xlabel('Iterations')
ylabel('Mean Intensity')
title('DASH Correction')

disp('Generating corrected R-Matrix...')
phscorrect=exp(1i*angle(Cp(:)));
for ii=1:size(T,2)
T(:,ii)=T(:,ii).*phscorrect(ii);
end
disp('Done!!!')





function f=fftfreq(npts,dt,alias_dt)
% returns a vector of the frequencies corresponding to the length
% of the signal and the time step.
% specifying alias_dt > dt returns the frequencies that would
% result from subsampling the raw signal at alias_dt rather than
% dt.
    
    
  if (nargin < 3)
      alias_dt = dt;
  end
  fmin = -1/(2*dt);
  df = 1/(npts*dt);
  f0 = -fmin;
  alias_fmin = -1/(2*alias_dt);
  f0a = -alias_fmin;
  
  ff = mod(linspace(0, 2*f0-df, npts)+f0,  2*f0)  - f0;
  fa = mod(                        ff+f0a, 2*f0a) - f0a;
  %  return the aliased frequencies
  f = fa;
end
end


