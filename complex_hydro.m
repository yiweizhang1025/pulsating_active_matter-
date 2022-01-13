% Parameters
N = 128;
L = 100;
rho = 1.3;
ep = 10;
D2 = 1;
om = 1e1;
lam = 5e-2;
C2_n = 12;
D = 4000;
dt = 1e-3;
interval_record = 1e-1;
initial_time = 4e2;
tmax = .5e0;
C2 = rho^C2_n;
%C2 = 0;
Nframe = floor(tmax/dt);

% Initial condition%
%u = rand(N,N)-0.5;
%v = rand(N,N)-0.5;
f0 = rho*ones(N,N)+1e-1*(rand(N,N)-0.5);
f = complex(u,v);
fbar = conj(f);
f2 = ep*f.^2/complex(4*D2,2*om);
f2bar = conj(f2);
t = -initial_time;
x = linspace(0,2*pi*L*(1-1/N),N);
y = linspace(0,2*pi*L*(1-1/N),N);
[X,Y] = meshgrid(x,y);
tlist = [];
Alist = [];
thlist = [];
Avarlist = [];
thvarlist = [];
videoisopen = 0;
figure1=figure('Position', [1024, 1024, 1024, 1024]);

% Check time first
Nx = N;
Ny = N;
tic; fft2(randn(3*Nx/2,3*Ny/2)); toc;
tic; fft2(randn(Nx,Ny)); toc;
tic; ifft2(fft2(randn(3*Nx/2,3*Ny/2)),'symmetric'); toc;
tic; ifft2(fft2(randn(Nx,Ny)),'symmetric'); toc;

% Optimize
fftw('dwisdom',[]);
fftw('planner','exhaustive');
fft2(randn(3*Nx/2,3*Ny/2));
fft2(randn(Nx,Ny));
ifft2(fft2(randn(3*Nx/2,3*Ny/2)),'symmetric');
ifft2(fft2(randn(Nx,Ny)),'symmetric');
fftinfo = fftw('dwisdom');
fftw('dwisdom',fftinfo);

% Then check again
tic; fft2(randn(3*Nx/2,3*Ny/2)); toc;
tic; fft2(randn(Nx,Ny)); toc;
tic; ifft2(fft2(randn(3*Nx/2,3*Ny/2)),'symmetric'); toc;
tic; ifft2(fft2(randn(Nx,Ny)),'symmetric'); toc;

% Iterations
while t<tmax
    % (0) Fourier Transform
    fk = fft2(f);
    f0k = fft2(f0);
    % (1) \tilde{f}(t+dt/2) from tilde{f}(t)
    qx = 1i * [0:N/2 (-N/2+1):-1].' /L;
    qy = 1i * [0:N/2 (-N/2+1):-1].'/L;
    q2 = -qx.*qx - qy.*qy;
    for m = 1:N
        kx = qx(m);
        for n=1:N
            ky = qy(n);
            k2 = -kx*kx-ky*ky;
            % convolution
            %f0k(m,n) = f0k(m,n)*(2-D*k2*dt/2)/(2+D*k2*dt/2)- C2/(1+lam)^2*((1+lam^2/2)*(cf0f0(m,n)+cf0f02(m,n))-lam*2*(imag(cf1f0(m,n)+cf0f12(m,n))))*dt/(2*(1+D*k2*dt/4));
            %f1k(m,n) = f1k(m,n)*(2-D*k2*dt/2)/(2+D*k2*dt/2)- C2/(1+lam)^2*((1+lam^2/2)*(cf1f0(m,n)+cf1f02(m,n)))*dt/(2*(1+D*k2*dt/4));
            % product
            f0k(m,n) = f0k(m,n)*(2-D*k2*dt/2)/(2+D*k2*dt/2);
            fk(m,n) = fk(m,n)*(2-D*k2*dt/2)/(2+D*k2*dt/2);
        end
    end
    t = t+dt/2;
    % (2) f(t+dt/2) from \tilde{f}(t+dt/2)
    f = ifft2(fk);
    fbar = conj(f);
    f0 = ifft2(f0k);
    f2 = ep*f.^2/complex(4*D2,2*om);
    f2bar = conj(f2);
    A2 = f.*fbar;
  
    % (3) f(t+dt) from f(t+dt/2)
    df = complex(-D2,-om)*f+ep/2*(-ep*A2.*f/complex(4*D2,2*om)+f0.*f)+complex(0,1)*C2*lam*f/(1+lam)^2.*(real(f)+lam/2*imag(f2))+sqrt(2*D2*rho)*complex(normrnd(0,1,[N,N]),normrnd(0,1,[N,N]));
    f = f + df*dt/2;
    %for m = 1:N
    %    kx = qx(m);
     %   for n=1:N
     %       ky = qy(n);
     %       k2 = -kx*kx-ky*ky;
     %       f0k(m,n) = f0k(m,n)*(2-D*k2*dt/2)/(2+D*k2*dt/2);
     %   end
    %end
    t = t+dt/2;
    fbar = conj(f);
    %f0 = ifft2(f0k);
    % (4) Output
    if t > 0
        % Initialise videos
        if videoisopen == 0
            myVideo = VideoWriter(sprintf('noisy_complex_hydro_N_%i_L_%i_rho_%g_ep_%i_D_%i_lam_%g_C_%i', N,L,rho,ep,D,lam,C2));%open video file
            myVideo.FrameRate = 20;  %can adjust this, 5 - 10 works well for me
            open(myVideo);
            n_record = 1;
            videoisopen = 1;
        end
        if t > n_record*interval_record
            %u = real(f);
            %v = imag(f);
            %A = (u.^2+v.^2).^(0.5);
            %theta = atan(v./A);
            A = abs(f);
            theta = angle(f);
            t
            % 1st plot for amp
           % subplot(2,2,1)
           % pcolor(A);
           % title('|f_1(x,y)|');
           % axis square;
           % colormap jet;
           % caxis([0 4]);
           % colorbar('eastoutside');
            
            % 2nd plot for phase
            %subplot(2,2,2)
            h = pcolor(theta);
            set(h, 'EdgeColor', 'none');
           % title('\Theta(x,y)');
            axis square;
            xticks([]);
            yticks([]);
            colormap jet;
            caxis([-pi,pi]);
            colorbar('eastoutside','XTick', -pi:2*pi:pi,'XTickLabel',{'-\pi','\pi'},'Fontsize',36)
            grid off
            
            % 3rd plot
            %A0 = mean(mean(A));
            %Alist(end+1) = A0;
            %tlist(end+1) = t;
            %subplot(2,2,3)
            %plot(tlist,Alist,'LineWidth',1.5);
            %axis square;
            %xlabel('t');
            %ylabel('|f_1|');
            %ylim([0 4])
            
            % 4th plot
            %th0 = mean(mean(theta));
            %thlist(end+1) = th0;
            % subplot(2,2,4)
            %plot(tlist,thlist,'LineWidth',1.5);
            %axis square;
            %xlabel('t');
            %ylabel('\Theta');
            %ylim([-pi pi])
            % 5rd plot
           % Avar = sqrt(var(reshape(A,[1,N*N])));
           % Avarlist(end+1) = Avar/A0;
           % subplot(3,2,3)
           % plot(tlist,Avarlist,'LineWidth',1.5);
           % axis square;
           % xlabel('t');
           % ylabel('std|f_1|/|f_1|');

            % 6th plot
           % thvar = sqrt(var(reshape(theta,[1,N*N])));
           % thvarlist(end+1) = thvar;
           % subplot(3,2,4)
           % plot(tlist,thvarlist,'LineWidth',1.5);
           % axis square;
           % xlabel('t');
           % ylabel('std(\Theta)');
            
            % Animation
            frame = getframe(gcf);
            writeVideo(myVideo,frame);
            n_record = n_record + 1;
        end
    end
end
saveas(gcf,'defect.png');
fileID = fopen('defect.dat','w');
for i=1:N
    x = i*L/N;
    for j=1:N
        y=j*L/N;
        fprintf(fileID,'%g\t%g\t%g\n',x,y,theta(i,j));
    end
end
close(myVideo);