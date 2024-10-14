%creating crystal
clear
% clc
%run('../startup.m')
%run('../../chunkie/startup.m')

fprintf('Setting up domain...\n')

Nx = 20;
hx = 1/Nx;
Ny = 10;
hy = 1/Ny;
x = -1:hx:1;
y = -1:hy:1;

fsave = ['result_mscat_volume_Nx' num2str(Nx) '_Ny' num2str(Ny) ...
      '_zk6piNx_gau_ref0_nref4.mat'];

a = hx/3;% this means the separation is hx/3
b = hy/3;% this means the separation is hy/3

% the rectangle should have then hal-size [hx/3+hx/9,hy/3+hy/9]

shift=[0,0];
foscx = @(x,y,c) 0.1*cos(3*(x-c(1))./a/2*pi) + 0.5*cos(10*(x-c(1))./a/2*pi) + sin(17*(x-c(1))./a/2*pi);
foscy = @(x,y,c) 0.1*sin(2*(y-c(2))./b/2*pi) + 0.5*cos(9*(y-c(2))./b/2*pi) + sin(19*(y-c(2))./b/2*pi);

qfun1 = @(x,y,c) -(2 + foscx(x,y,c).*foscy(x,y,c)).*exp(-50*((x-c(1)).^2/4/a.^2+ (y-c(2)).^2/4/b.^2));

ifplot = 0;
n = 16;
nref = 4;

if ifplot
    rect_use = [-a, a, -b, b];
    uv = surfacemesh.square(n, nref, rect_use);
    u = uv.x;
    v = uv.y;
    xdom = cell(size(u));
    ydom = cell(size(u));
    zdom = cell(size(u));
    for k=1:length(uv)
        xdom{k} = u{k};
        ydom{k} = v{k};
        zdom{k} = 0*u{k};
    end
    dom = surfacemesh(xdom, ydom, zdom);
    q = surfacefun(@(x,y,z) qfun1(x,y,[0,0]), dom);
end


%%

zk = 6*pi*Nx;
fprintf('Setting up kernels...\n')
Sk = kernel('helm', 's', zk);
Skp = kernel('helm', 'sprime', zk);
    
Dk = kernel('helm', 'd', zk);
Dkp = kernel('helm', 'dprime', zk);

%loop for self convergence
fprintf('Doing self-convergence...\n')
for iself=1:1

    
    fprintf('Wavenumber=%d\n',zk)

    %set-up kernel and matrix A

    
    fprintf('Self convergence %d...\n',iself)
    fprintf('Setting up proxy rectangles...\n')
    tic
    % this has to be fixed. ppw isnot the idea here.
    npt_rect= 30*(2^(iself-1));
    Nppw_hor = npt_rect;
    Nppw_ver = 2*npt_rect;
    npxys = [Nppw_hor, Nppw_ver];    

    %fmm eps and gmres
    % eps_gmres = 10^(-(6+(iself-1)*3));
    % eps_fmm = 1e-14;
    eps_gmres = 1e-9;
    eps_fmm = 1e-9;

    opts = [];
    opts.iflege = 1;
    ls = [a+hx/9,b+hy/9];
    [pr,ptau,pw,pin] = chnk.flam.proxy_rect_pts(ls, npxys, opts);
    pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
    sqpw = sqrt(pw);

    plot(pr(1,:), pr(2,:), 'k.')
    

    %rectangle base
    pinfo = [];
    pinfo.r = pr;
    pinfo.n = pn;
    pinfo.wts = pw;
    Time_base_r(iself) = toc;
    fprintf('Time to set up  base rect=%d\n',Time_base_r(iself))  

    %setting up scattering matrices
    fprintf('Time scattering matrix\n')
    tic
    opts_L = [];
    opts_L.n = n;
    opts_L.nref = nref;
    xylim = [-a, -b; a, b];
    S = get_volume_scattering_matrix(zk, @(x,y) qfun1(x,y,[0,0]), ...
      xylim, pinfo, opts_L);
    Time_scat_mat(iself) = toc;
    fprintf('Time to calculate scat mat=%d\n',Time_scat_mat(iself))    

    fprintf('Finding correction Tself for Tmat\n')
    tic
    Dkmat = sqpw'.*Dk.eval(pinfo, pinfo).*sqpw;
    Skmat = sqpw'.*Sk.eval(pinfo, pinfo).*sqpw;
    Dkpmat = sqpw'.*Dkp.eval(pinfo, pinfo).*sqpw;
    Skpmat = sqpw'.*Skp.eval(pinfo, pinfo).*sqpw;
        
    Tself = [Dkmat, -Skmat; ...
             Dkpmat, -Skpmat];
    Tself(isnan(Tself)) = 0;

    Time_corrmat(iself) = toc;
    fprintf('Time to calculate Tmat correction mat=%d\n',Time_corrmat(iself))    

    %find crystal
    %construct particles and rectangles
    fprintf('Building domain - particles and rectangles...\n')
    tic
    num_obj = 1;
    cryst = cell(length(x), length(y));
    rect = cell(length(x), length(y));
    for ii=1:length(x)
        for jj=1:length(y)

            if mod(ii,2) 
               shift=[x(ii),y(jj)];
            else
                shift=[x(ii),y(jj)+hy/2];
            end

           
            cryst{ii,jj}.include = 1;

            % for 20 and 10
            if  (ii>=15) && (ii <=16) && (jj>12)
                cryst{ii,jj}.include = 0;
            end
             
            if  (jj==12) && (ii<=16)% && (jj>=4)
                cryst{ii,jj}.include = 0;
            end

            % for 10 and 5
            %if  (ii>7) && (ii <=8) && (jj>6)
            %    cryst{ii,jj}.include = 0;
            %end

            %if  (jj==6) && (ii<=8)% && (jj>=4)
            %    cryst{ii,jj}.include = 0;
            %end

            pinfo = [];
            pinfo.r = pr + shift';
            pinfo.n = pn;
            pinfo.wts = pw;
                
            rect{ii,jj}.pinfo = pinfo;        

            if cryst{ii,jj}.include == 1
                list_rect{num_obj} = rect{ii,jj}.pinfo;
                num_obj = num_obj + 1;
            end

        end
    end
    Time_dom(iself) = toc;
    fprintf('Time to set-up domain=%d\n',Time_dom(iself))

    % set plane wave data
    % Set up plane wave data as sources;
    fprintf('Setting up plane wave data...\n')
    tic
    thet = 0;%rand*2*pi;
    ct = cos(thet); sd = sin(thet);

    udata_pxy = [];
    for ii=1:length(list_rect)
    
        pinfo = list_rect{ii};

        px = pinfo.r(1,:).';
        py = pinfo.r(2,:).';

        pnx = pinfo.n(1,:).';
        pny = pinfo.n(2,:).';

        pw = pinfo.wts;
        pw = pw(:);
        sqpw = sqrt(pw);

        uincs = exp(1j * zk * (px*ct + py*sd) );
        gxuincs = 1j * zk * ct .* uincs;
        gyuincs = 1j * zk * sd .* uincs;

        dudnincs = pnx .* gxuincs + pny .* gyuincs;

        uin_pxy = [uincs .* sqpw(:); 
                   dudnincs .* sqpw(:)]; 
        
        udata_pxy =[udata_pxy;
                    S*uin_pxy];

    end
    Time_rhs(iself) = toc;
    fprintf('Time to calculate rhs=%d\n',Time_rhs(iself))    

    fprintf('Solving problem...\n')        
    
    afun = @(x) SSolve_function(x, zk, list_rect, eps_fmm, S, Tself);

    x_test = rand(length(udata_pxy),1);

    tic
    y_test = afun(x_test);
    toc

    Nit = ceil(0.1*length(list_rect)*length(pinfo.r(1,:)));
    Nit = 4000;
    tic
    uout_pxy_gmres = gmres(afun, udata_pxy, [], eps_gmres, Nit);        
    Time_dens(iself) = toc
    fprintf('Time to calculate densities=%d\n',Time_dens(iself))

    %calculating field
    tic
    targinfo = [];
    t = -2:0.01:2;
    [X,Y] = meshgrid(t);
    % targinfo.r = [5;5];
    % targinfo.n = [1;0];%this is unecessary
    targinfo.r = [X(:)';Y(:)'];
    targinfo.n = [1;0];%this is unecessary
    
    out = calculate_field_fmm(eps_fmm, zk, list_rect, targinfo, uout_pxy_gmres);
    soltarg2(iself,:) = out;
    Time_calcfield(iself) = toc;
    fprintf('Time to get field in the domain fmm=%d\n',Time_calcfield(iself))

    save(fsave,'-v7.3')
end %iself

%Self convergence
fprintf('Self convergence series...\n')
for ii = 1: size(soltarg2,1)-1
     norm(soltarg2(ii,1)-soltarg2(end,1))/norm(soltarg2(end,1))
end

save(fsave,'-v7.3')

%%
uin = exp(1j*zk*X);
uin = reshape(uin, size(soltarg2(end,:)));
pdata = soltarg2(end, :) + uin;
xuse = X(:);
yuse = Y(:);
for ii=1:length(list_rect)
    pinfo = list_rect{ii};
    xr = pinfo.r(1,:);
    yr = pinfo.r(2,:);
    
    xm = mean(xr);
    ym = mean(yr);
    
    xx = [xm - ls(1)/2, xm + ls(1)/2, xm + ls(1)/2, xm - ls(1)/2];
    yy = [ym - ls(2)/2, ym - ls(2)/2, ym + ls(2)/2, ym + ls(2)/2];
    in = inpolygon(xuse, yuse, xx, yy);
    pdata(in) = NaN;
   
end
pdata = reshape(pdata, size(X));
pcolor(X, Y, abs(pdata)); shading interp;

function out = SSolve_function(x,zk,list_rect,eps,S,Tself)
%     Ssolve = eye(nn) + (Stotal-eye(nn))*Tmat;
%     Stotal = [S1, zeros(mm);
%              zeros(mm),  S2];
% Si -> scattering matrix

% fprintf('test')
    % eps = 1e-6;%1e-14;       
    % this is considering a unique rectangle
    % I am not sure if this will change with rotation
    % have to think
    pw = list_rect{1}.wts;
    pw = pw(:);
    sqpw = sqrt(pw);    
    npanels = length(list_rect);
    nnn = length(sqpw);
    
    % here we sue different Tmat from manas
    aux = Tmat_fast_matvec(eps, zk, npanels, nnn, list_rect, sqpw, Tself, x);    
    % aux = Tmat_fmm_multip(eps,zk,npanels,nnn,list_rect,sqpw,x);
    
    numobj = length(list_rect);
    aux1 = [];    
    for ii=1:numobj
        nS2 = size(S,2);
        aux1 = [aux1;
                S*aux(nS2*(ii-1)+1:nS2*ii) + aux(nS2*(ii-1)+1:nS2*ii)];
        
    end

    out = x - aux1;

end

function out = calculate_field_fmm(eps,zk,plist,targ,density)

    srcuse = [];
    srcuse.sources = [];
    srcuse.charges = [];
    srcuse.dipstr  = [];
    srcuse.dipvec  = [];
    
    target = targ.r(1:2,:); 

    sqpw = sqrt(plist{1}.wts(:));

    npols = length(plist); 
    
    dens = density .* repmat(sqpw,2*npols,1);
    
    nnn = size(plist{1}.r,2);
    
    % eps = 1e-14;
    
    for ip = 1 : npols
        srcuse.sources = [srcuse.sources plist{ip}.r(1:2,:)];
        srcuse.charges = [srcuse.charges -transpose(dens(2*nnn*(ip-1)+nnn+1:2*nnn*ip))];
        srcuse.dipstr  = [srcuse.dipstr  transpose(dens(2*nnn*(ip-1)+1:2*nnn*(ip-1)+nnn))];
        srcuse.dipvec  = [srcuse.dipvec  plist{ip}.n(1:2,:)];
    end
    
    pg  = 0;
    pgt = 1;
    
    U = hfmm2d(eps, zk, srcuse, pg,target,pgt);
    
    out = transpose(U.pottarg);

end
