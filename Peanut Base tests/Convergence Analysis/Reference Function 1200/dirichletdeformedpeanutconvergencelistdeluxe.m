% Computes Dirichlet eigenvalues and normed eigenfunctions for polygons.
% (NOTE: Output are wave numbers, i.e. have to be squared for eigenvalues) Input parameters:
% c = Number of corners (5 = Pentagon, 6 = hexagon etc.)
% Computing wavenumbers and their eigenfunctions from xstart to xend,
% covering the interval in adjacent circles of radius intervallength.
% Example: dirichletmultigon(5,6,10,1)

function eigenvaluearray=dirichletdeformedpeanutconvergencelistdeluxe(K) % K is the number of the eigenfunction tested
    global alpha
    global l
    global m
    global N
    global tol_rank
    global mp
    global R
    global faces
    global n
    global vi
    global vhi
    global wavenumber
    format long


    %% Generating/Importing eigenfunctions Bump function coefficients 
    mainfolder = fileparts(fileparts(pwd));
    datafolder = strcat(mainfolder,'/Data');
    programsfolder = strcat(mainfolder,'/Programs');
    addpath(genpath(programsfolder)); % Adding the programs folder
    wavenum = importdata(strcat(datafolder,'/PeanutWavenumbers.mat'));
    wavenum = real(wavenum);
    ev = wavenum.^2;     
    
    %% Defining checked interval and variables
    evlen = length(ev);
    mp = ev(K);
    if K==1
        R = (ev(2)-ev(1))/2;
    elseif K==evlen
        R = (ev(end)-ev(end-1))/2;
    else
        R = min(ev(K+1)-ev(K),ev(K)-ev(K-1))/2;
    end
    leftend = mp - R;
    rightend = mp + R;
    fprintf("Checking interval: [%.3f,%.3f]\n",mp-R,mp+R)
    plotarray = [];
    eigenvaluearray = [];
    nvector = [2400,2000,1800,1500,1200,1000,800,600,400,300,200,150,100];
    reso=81;
    nlen = length(nvector);
    disp(nvector)
    facevector = nvector/2;

    arrays = zeros(reso,reso,nlen);  %Arrays containing all the plot arrays for different n (ready to integrate)
    eigenvaluevector = zeros(1,nlen);
    n = 100;
    xx = linspace(0,2*pi,n+1);
    xfun = @(x) 0.06*((cos(x)+2).*(cos(x+0.6)+2).*(2+0.1*(cos(3*x))))-0.1;
    vix = xfun(xx);
    yfun = @(y) 0.06*((sin(y)+2).*(sin(y-0.5)+2).*(2+0.4*(cos(2*y))).*(1+0.1*(sin(4*y))))-0.06;
    viy = yfun(xx);
    figure(1)
    scatter(vix,viy)

    for nn = 1:nlen
        n=nvector(nn);
        fprintf("n=%d\n",n)
        alpha=(1-sqrt(3/5))/2;
        faces=n/2;
        fprintf('Number of faces: %d\n',faces);
        fprintf('Number of collocation nodes: %d\n',3*faces);
        N = 24; 
        l=8;
        tol_rank = 10^(-4);
        xx = linspace(0,2*pi,n+1);
        xfun = @(x) 0.06*((cos(x)+2).*(cos(x+0.6)+2).*(2+0.1*(cos(3*x))))-0.1;
        vix = xfun(xx);
        yfun = @(y) 0.06*((sin(y)+2).*(sin(y-0.5)+2).*(2+0.4*(cos(2*y))).*(1+0.1*(sin(4*y))))-0.06;
        viy = yfun(xx);
        vi = [vix;viy];
        outer = vi;
        indi1 = [1:2:n];
        indi2 = [2:2:n];
        indi3 = [3:2:n+1];
        vhi = vi;
        k=1;
        L=1;
        for i=1:faces
            x1=vi(1,indi1(i));
            y1=vi(2,indi1(i));
            x2=vi(1,indi2(i));
            y2=vi(2,indi2(i));
            x3=vi(1,indi3(i));
            y3=vi(2,indi3(i));
            s=alpha;
            vhi(1,L)=l1(s)*x1+l2(s)*x2+l3(s)*x3;
            vhi(2,L)=l1(s)*y1+l2(s)*y2+l3(s)*y3; 
            s=0.5;
            vhi(1,L+1)=l1(s)*x1+l2(s)*x2+l3(s)*x3;
            vhi(2,L+1)=l1(s)*y1+l2(s)*y2+l3(s)*y3; 
            s=1-alpha;
            vhi(1,L+2)=l1(s)*x1+l2(s)*x2+l3(s)*x3;
            vhi(2,L+2)=l1(s)*y1+l2(s)*y2+l3(s)*y3; 
            k=k+2;
            L=L+3;
        end
        m=3*faces;
        fprintf("Computing eigenvalues.\n")
        tic    
        [eigenvalues,V0] = IntegralAlgo(@matr);
        len = length(eigenvalues);
        start = 0;
        fprintf("Eigenvalue array before elimination:\n")
        disp(eigenvalues)
        i=1;
        while start<len
            if real(eigenvalues(i))<leftend || real(eigenvalues(i))>rightend || abs(imag(eigenvalues(i)))>0.1*abs(real(eigenvalues(i)))
                eigenvalues(i)=[];
            else
            i = i+1;
            end
            start = start+1;
        end
        s = size(eigenvalues);
        if s(1)==0 || s(2)==0
            eigenvalues=[];
        end
        eigenvalues = sort(eigenvalues);
        fprintf("Eigenvalue array after elimination:\n")
        disp(eigenvalues)
        fprintf("Computing eigenfunction.\n")

        %Create plot arrays for the eigenfunctions
        eigenvaluearray = [eigenvaluearray;eigenvalues];
        if isempty(eigenvalues)==0
            looplength=0;
        else
            looplength=1;
        end
        for iii=1:looplength
            V01=V0(:,iii);
            xvec=linspace(0,1,reso);
            yvec=linspace(0,1,reso);
            [xx,yy]=meshgrid(xvec,yvec);
            z = zeros(reso);
            zint = zeros(reso);
            zreal = zeros(reso);
            for i=1:reso
                for j=1:reso
                    x=xx(i,j);
                    y=yy(i,j);
                    [In,On]=inpolygon(x,y,outer(1,:),outer(2,:));
                    if On==1
                        z(i,j)=0;
                        zint(i,j)=0;
                    elseif In==1
                        P=[x;y];
                        summ=0;
                        k=1;
                        L=1;
                        for jj=1:n/2    % number of faces
                        % extract the first, middle, and endpoint
                            x1=vi(1,indi1(jj));
                            y1=vi(2,indi1(jj));
                            x2=vi(1,indi2(jj));
                            y2=vi(2,indi2(jj));
                            x3=vi(1,indi3(jj));
                            y3=vi(2,indi3(jj));
                            k=k+2;
                            E1=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),eigenvalues(iii),1);
                            E2=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),eigenvalues(iii),2);
                            E3=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),eigenvalues(iii),3);
                            summ=summ+E1*V01(L)+E2*V01(L+1)+E3*V01(L+2);
                            L=L+3;
                        end
                        z(i,j) = summ;
                        zint(i,j) = real(summ)^2;
                        zreal(i,j) = real(summ);
                    else
                        z(i,j)=NaN;
                        zint(i,j)=0;
                        zreal(i,j) = 0;
                    end
                end % end loop
            end % end loop
        xcoord=linspace(0,1,reso);
        ycoord=linspace(0,1,reso);
        norm = sqrt(trapz(ycoord,trapz(xcoord,zint,1),2));  %Double integral with trapezoidal rule
        sprintf("Norm = %d",norm)
        z = z/norm;
        zreal = zreal/norm;
          
        %%%% Sign check with the most accurate function
        % The sign is checked for three points in the array. The check
        % triggers (i.e. sign is flipped) if it is different in one of the
        % points in which the absolute value of the functions is greater
        % than 0.1 (to avoid floating point mistakes)
        mid = round(reso/2);
        if nn==1
            sign1 = sign(zreal(mid,mid));             %Point in the middle
            sign2 = sign(zreal(mid,round(1.5*mid)));  %Point on the x-axis to the right
            sign3 = sign(zreal(round(1.3*mid),round(1.3*mid))); %Point diagonally up right
        end
        if nn>1
            if (sign(zreal(mid,mid))~=sign1 && abs(zreal(mid,mid))>0.1) || ...
               (sign(zreal(mid,round(1.5*mid)))~=sign2 && abs(zreal(mid,round(1.5*mid)))>0.1) || ...
               (sign(zreal(round(1.3*mid),round(1.3*mid)))~=sign3 && abs(zreal(round(1.3*mid),round(1.3*mid)))>0.1)
               zreal = -zreal;
            end  % The sign is checked at three different points. If it is different at one of them AND 
        end      % the modulus of the function there is >0.1, the sign is flipped
          
        plotarray = cat(3,plotarray,z);
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        set(fig,'Visible', 'off');
        hold on
        imagesc(xvec,yvec,real(z))
        plot(vix,viy,'k','LineWidth',3.5);
        axis equal
        axis([0 1 0 1])
        hold off
        set(gca,'Ydir','normal')
        xlabel('x')
        ylabel('y')
        minimum = min(min(real(z)));
        maximum = max(max(real(z)));
        caxis( [minimum-0.1 maximum] )
        cmap = [1 1 1; parula(10000)];
        colormap(cmap);
        colorbar;
        ticks = 12;
        set(gca,'fontsize',ticks,'fontweight','bold','LineWidth',2);
        title(sprintf('Eigenfunction %.9f',eigenvalues(iii)))
        saveas(gcf,sprintf('Peanut Eigenfunction %03d %.8f, %d faces.png',K,eigenvalues(iii),faces))
        end 
        if isempty(eigenvalues)~=0
            eigenvaluevector(nn) = eigenvalues(1);     
        end
        arrays(:,:,nn) = zreal;
    end

    differencearray = zeros(reso,reso,nlen);
    normdifferences = zeros(1,nlen);
    eigenvaluedifferences = zeros(1,nlen);

    for i1=2:nlen
        differencearray(:,:,i1) = arrays(:,:,i1)-arrays(:,:,1);
    end

    for i1=2:nlen
        intarray = differencearray(:,:,i1).^2;
        normdifferences(i1) = sqrt(trapz(ycoord,trapz(xcoord,intarray,1),2)); 
    end

    for i1=2:nlen
        eigenvaluedifferences(i1) = abs(eigenvaluevector(i1)-eigenvaluevector(1));
    end

    normdifferences(1)=[];
    eigenvaluedifferences(1)=[];
    facevector(1)=[];

    save(sprintf('Peanut L2 %03d %.8f facemax=%d reso=%d.mat',K,eigenvaluevector(1),facevector(1),reso),'normdifferences','facevector')
    save(sprintf('Peanut ev %03d %.8f facemax=%d reso=%d.mat',K,eigenvaluevector(1),facevector(1),reso),'eigenvaluedifferences','facevector')
end

function A=matr(k)
  global wavenumber

  wavenumber=k;
  A = createDmatrix();
  A=-0.5*eye(size(A))+A;
end



function [eigenvalues,Z] = IntegralAlgo(func)
    global l
    global N
    global tol_rank
    global m

    ii = 0:1:N;
    tj = 2*pi*ii/N;

    % zufaellige Matrix VDach anlegen
    vD = eye(m,m);
    % A0,N und A1,N bestimmen
    A0 = 0;
    A1 = 0;
    for i=1:N
        p = phi(tj(i));
        ablP = ablPhi(tj(i));
        T=func(p);
        %cond(T)
        A0 = A0 + T\(vD * ablP); % i-ter Summand von A0
        A1 = A1 + T\(vD * p * ablP); % i-ter Summand von A1
    end
    A0 = A0 * 1/(1i*N);
    A1 = A1 * 1/(1i*N);
    % SVD von A0,N
    [V,SW,WH] = svd(A0);
    s = diag(SW);
    k = 0;
    % Rang Test von s
    for i=1:m
        if s(i)> tol_rank
            k=k+1;
        else
            break;
        end
    end
    % falls l=k l erhoehen und wiederholen
    %if l==k && ~(l==m)
    %    l=l+1;
    V0 = V(1:m,1:k);
    W0 = WH(1:m,1:k);
    s0 = diag(s(1:k)); 
    % B berechnen mit B = V0
    B = (V0' * A1 * W0) /s0;
    % Eigenwert-Problem fuer B berechnen
    [V,D] = eig(B);
    eigenvalues=diag(D);
    Z=V0*V;
end

% Funktion fuer Phi
function f=phi(t)
    global R
    global mp
    f =  mp + R * exp(1i*t);
    %f = a * cos(t) + b* sin(t);
end

% Funktion fuer die Ableitung von phi
function f = ablPhi(t)
    global R
    f = 1i * R * exp(1i*t);
end

function indi=piminus(j,n)
    if j==n
        indi=1;
    else
        indi=j+1;
    end
end

function A=createDmatrix()
    global wavenumber
    global vi
    global vhi
    global faces
    global n
    
    indi1=[1:2:n];
    indi2=[2:2:n];
    indi3=[3:2:n+1];
    % create the matrix of size n x n
    A=zeros(3*faces,3*faces);
    M0=zeros(3*faces,3*faces);
    for i=1:3*faces
        % extract collocation point of interval I_i
        xh=vhi(1,i);
        yh=vhi(2,i);
        k=1;
        L=1;
        for j=1:n/2    % number of faces
            % extract the first, middle, and endpoint
            x1=vi(1,indi1(j));
            y1=vi(2,indi1(j));
            x2=vi(1,indi2(j));
            y2=vi(2,indi2(j));
            x3=vi(1,indi3(j));
            y3=vi(2,indi3(j));
            k=k+2;
            if i==L
              E1=0;
            else
              E1=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,1);
            end
            if i==L+1
              E2=0;
            else
              E2=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,2);
            end
            if i==L+2
              E3=0;
            else
              E3=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,3);
            end
            
            A(i,L  )=E1;
            A(i,L+1)=E2;
            A(i,L+2)=E3;
            if i==L
              E1=0;
            else
              E1=intDL0(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,1);
            end
            if i==L+1
              E2=0;
            else
              E2=intDL0(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,2);
            end
            if i==L+2
              E3=0;
            else
              E3=intDL0(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,3);
            end
            
            M0(i,L  )=E1;
            M0(i,L+1)=E2;
            M0(i,L+2)=E3;
            L=L+3;
        end
    end
    for i=1:3*faces
        A(i,i)=-0.5-sum(M0(i,:),2);
    end
end

% integrate elastic single layer boundary function
function E=intDL(x1s,x2s,x3s,y1s,y2s,y3s,xhs,yhs,ks,basisf)
    global xh
    global yh
    global x1
    global x2
    global x3
    global y1
    global y2
    global y3
    global basis
    global k

    xh=xhs;
    yh=yhs;
    x1=x1s;
    x2=x2s;
    x3=x3s;
    y1=y1s;
    y2=y2s;
    y3=y3s;
    k=ks;
    basis=basisf;

    E=quadgk(@c11_1,0,1);
end

function y=c11_1(s)
    global xh
    global yh
    global x1
    global x2
    global x3
    global y1
    global y2
    global y3
    global k
    global basis
  
    % Jacobian
    dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
    dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
    J=sqrt(dxds.^2+dyds.^2);

    % m_k~
    a=xh-(l1(s)*x1+l2(s)*x2+l3(s)*x3);
    b=yh-(l1(s)*y1+l2(s)*y2+l3(s)*y3);
    r=sqrt(a.^2+b.^2);

    if basis==1
        f=l1a(s);
    elseif basis==2
        f=l2a(s);
    else
        f=l3a(s);
    end

    ny1=dyds./J;
    ny2=-dxds./J;
    
    y=1i*k*besselh(1,k*r)./(4*r).*(a.*ny1+b.*ny2).*J.*f;  % J would cancel out
end
   
% integrate elastic single layer boundary function
function E=intDL0(x1s,x2s,x3s,y1s,y2s,y3s,xhs,yhs,ks,basisf)
    global xh
    global yh
    global x1
    global x2
    global x3
    global y1
    global y2
    global y3
    global basis
    global k

    xh=xhs;
    yh=yhs;
    x1=x1s;
    x2=x2s;
    x3=x3s;
    y1=y1s;
    y2=y2s;
    y3=y3s;
    k=ks;
    basis=basisf;

    E=quadgk(@c11_2,0,1);
end

function y=c11_2(s)
    global xh
    global yh
    global x1
    global x2
    global x3
    global y1
    global y2
    global y3
    global basis
  
    % Jacobian
    dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
    dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
    J=sqrt(dxds.^2+dyds.^2);

    % m_k~
    a=xh-(l1(s)*x1+l2(s)*x2+l3(s)*x3);
    b=yh-(l1(s)*y1+l2(s)*y2+l3(s)*y3);
    r=sqrt(a.^2+b.^2);

    if basis==1
        f=l1a(s);
    elseif basis==2
        f=l2a(s);
    else
        f=l3a(s);
    end

    ny1=dyds./J;
    ny2=-dxds./J;
    
    y=1./(2*pi*r.*r).*(a.*ny1+b.*ny2).*J.*f;  % J would cancel out
end
    
    % integrate elastic single layer boundary function
function E=intSL(x1s,x2s,x3s,y1s,y2s,y3s,xhs,yhs,ks,basisf,check)
    global xh
    global yh
    global x1
    global x2
    global x3
    global y1
    global y2
    global y3
    global basis
    global k
    global alpha

    xh=xhs;
    yh=yhs;
    x1=x1s;
    x2=x2s;
    x3=x3s;
    y1=y1s;
    y2=y2s;
    y3=y3s;
    k=ks;
    basis=basisf;

    if check==0
        E=quadgk(@c11_3,0,1);
    elseif check==1
        E=quadgk(@c11_3,0,alpha)+quadgk(@c11_3,alpha,1);
    elseif check==2
        E=quadgk(@c11_3,0,0.5)+quadgk(@c11_3,0.5,1);
    elseif check==3
        E=quadgk(@c11_3,0,1-alpha)+quadgk(@c11_3,1-alpha,1);
    else
        fprintf('Something is wrong!\n');
    end
end


function y=c11_3(s)
    global xh
    global yh
    global x1
    global x2
    global x3
    global y1
    global y2
    global y3
    global k
    global basis

    % Jacobian
    dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
    dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
    J=sqrt(dxds.^2+dyds.^2);

    % m_k~
    a=xh-(l1(s)*x1+l2(s)*x2+l3(s)*x3);
    b=yh-(l1(s)*y1+l2(s)*y2+l3(s)*y3);
    r=sqrt(a.^2+b.^2);

    if basis==1
        f=l1a(s);
    elseif basis==2
        f=l2a(s);
    else
        f=l3a(s);
    end

    y=1i*besselh(0,k*r)/4.*J.*f;
end

% first Lagrange basis function
function y=l1(s)
    u=1-s;
    y=u.*(2*u-1);
end

function y=l1a(s)
    global alpha
    
    u=1-s;
    y=(u-alpha)/(1-2*alpha).*(1-2*s)/(1-2*alpha); % check ok
end

% second Lagrange basis function
function y=l2(s)
    u=1-s;
    y=4*s.*u;
end

function y=l2a(s)
    global alpha
    
    u=1-s;
    y=4*(s-alpha)/(1-2*alpha).*(u-alpha)/(1-2*alpha); % check ok
end

% third Lagrange basis function
function y=l3(s)
    y=s.*(2*s-1);
end

function y=l3a(s)
    global alpha
    
    y=(s-alpha)/(1-2*alpha).*(2*s-1)/(1-2*alpha); % check ok
end
