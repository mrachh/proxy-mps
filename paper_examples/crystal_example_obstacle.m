%creating crystal
clear
% clc
%run('../startup.m')
%run('../../chunkie/startup.m')

fprintf('Setting up domain...\n')

Nx = 2;%20
hx = 1/Nx;
Ny = 1;%10
hy = 1/Ny;
x = -1:hx:1;
y = -1:hy:1;

fsave='result_mscat_wiggleN8_Nx20_Ny10_zk6piNy_new.mat';

a = hx/3;% this means the separation is hx/3
b = hy/3;% this means the separation is hy/3

% the rectangle should have then hal-size [hx/3+hx/9,hy/3+hy/9]
cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;

shift=[0,0];
%coefs for general ellipse
coefs=zeros(1,8);
coefs(1)  = 1;
coefs(end) = 0.1;

% chnkr0 = chunkerfunc(@(t)ellipse(t, a, b), cparams, pref); 
% chnkr0 = sort(chnkr0);

chnkr0 = chunkerfunc(@(t)wigglellipse(t, a, b, shift, coefs), cparams, pref); 
chnkr0 = sort(chnkr0);

%for specific ellipse of convo
% chnkr0 = chunkerfunc(@(t)wigglellipse1(t, a, b, shift, [1,0.1]), cparams, pref); 
% chnkr0 = sort(chnkr0);

zk = 6*pi*Nx;
fprintf('Setting up kernels...\n')
Sk = kernel('helm', 's', zk);
Skp = kernel('helm', 'sprime', zk);
    
Dk = kernel('helm', 'd', zk);
Dkp = kernel('helm', 'dprime', zk);
    
Ck = kernel('helm', 'c', zk, [1.0, -1j*zk]);

%loop for self convergence
fprintf('Doing self-convergence...\n')
for iself=1:4

    
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
    n = chnkr0.npt;
    A = chunkermat(chnkr0, Ck) + 0.5*eye(n);
    [S, C, Cp] = get_scattering_matrices(zk, chnkr0, pinfo, A);    
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
    for ii=1:length(x)
        for jj=1:length(y)

            if mod(ii,2) 
               shift=[x(ii),y(jj)];
            else
                shift=[x(ii),y(jj)+hy/2];
            end

            % cryst{ii,jj}.chnkr = chunkerfunc(@(t)ellipse(t, a, b, shift), cparams, pref); 
            cryst{ii,jj}.chnkr = chunkerfunc(@(t)wigglellipse(t, a, b, shift, coefs), cparams, pref); 
            % cryst{ii,jj}.chnkr = chunkerfunc(@(t)wigglellipse1(t, a, b, shift, [1,0.1]), cparams, pref); 
            cryst{ii,jj}.chnkr = sort(cryst{ii,jj}.chnkr);
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
                list_obj{num_obj}.chnkr = cryst{ii,jj}.chnkr;
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
    for ii=1:length(list_obj)
    
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
    
    afun = @(x)SSolve_function(x,zk,list_rect,eps_fmm,S,Tself);

    x_test = rand(length(udata_pxy),1);

    tic
    y_test = afun(x_test);
    toc

    Nit = ceil(0.1*length(list_rect)*length(pinfo.r(1,:)));
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
    
    out = calculate_field_fmm(eps_fmm,zk,list_rect,targinfo,uout_pxy_gmres);
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

% plot crystal with rectangles
% fprintf('Plot domain...\n')
%figure
%pcolor(X,Y,reshape(real(soltarg2(end,:)),size(X,1),size(X,1))); shading interp
%hold on;
%axis equal
%for ii=1:num_obj-1
%
%    plot(list_obj{ii}.chnkr,'k')
%    plot(list_rect{ii}.r(1,:),list_rect{ii}.r(2,:),'r')
%
%end
%
%uin = exp(1i*zk*X);



function [r,d,d2] = ellipse(t,a,b, varargin)
    shift = [0 0];
    if (nargin == 4)
        shift = varargin{1};
    end
    x0 = shift(1);
    y0 = shift(2);
    xs = x0 + a*cos(t);
    ys = y0 + b*sin(t);
    dxs = -a*sin(t);
    dys = b*cos(t);

    d2xs = -a*cos(t);
    d2ys = -b*sin(t);

    r = [(xs(:)).'; (ys(:)).'];
    d = [(dxs(:)).'; (dys(:)).'];
    d2 = [(d2xs(:)).'; (d2ys(:)).'];
end

function [r,d,d2] = wigglellipse(t,a,b,varargin)
    shift = [0 0];
    rt   = 1;
    drt  = 0;
    d2rt = 0;
    if (nargin == 4)
        shift = varargin{1};
    end

    if (nargin == 5)
        shift = varargin{1};
        coefs = varargin{2};
        N = length(coefs);
        rt   = coefs(1)+cos(bsxfun(@times,t,1:N-1))*coefs(2:N)';
        drt  = bsxfun(@times,(1:N-1)',-sin(bsxfun(@times,(1:N-1)',t')))'*coefs(2:N)';
        d2rt = bsxfun(@times,((1:N-1).*(1:N-1))',-cos(bsxfun(@times,(1:N-1)',t')))'*coefs(2:N)';
    end

    x0 = shift(1);
    y0 = shift(2);
    ct = cos(t);
    st = sin(t);

    xs = x0 + a*rt.*ct;
    ys = y0 + b*rt.*st;
    dxs = -a*rt.*st+a*drt.*ct;
    dys = b*rt.*ct+b*drt.*st;

    d2xs = - a*rt.*ct - 2 * a*drt.*st + a*d2rt.*ct;
    d2ys = - b*rt.*st + 2 * b*drt.*ct + b*d2rt.*st;

    r = [(xs(:)).'; (ys(:)).'];
    d = [(dxs(:)).'; (dys(:)).'];
    d2 = [(d2xs(:)).'; (d2ys(:)).'];
end

function [r,d,d2] = wigglellipse1(t,a,b,varargin)
    shift = [0 0];
    rt   = 1;
    drt  = 0;
    d2rt = 0;
    if (nargin == 4)
        shift = varargin{1};
    end

    if (nargin == 5)
        shift = varargin{1};
        coefs = varargin{2};       
        rt   = coefs(1)+cos(20*t)*coefs(2);
        drt  = -20*sin(20*t)*coefs(2);
        d2rt = -400*cos(20*t)*coefs(2);
    end

    x0 = shift(1);
    y0 = shift(2);
    ct = cos(t);
    st = sin(t);

    xs = x0 + a*rt.*ct;
    ys = y0 + b*rt.*st;
    dxs = a * (-rt.*st + drt.*ct);
    dys = b * ( rt.*ct + drt.*st);

    d2xs = a * ( - rt.*ct - 2 * drt.*st + d2rt.*ct);
    d2ys = b * ( - rt.*st + 2 * drt.*ct + d2rt.*st);

    r = [(xs(:)).'; (ys(:)).'];
    d = [(dxs(:)).'; (dys(:)).'];
    d2 = [(d2xs(:)).'; (d2ys(:)).'];
end

function [S, C_chnkr_to_pxy, Cp_chnkr_to_pxy] = get_scattering_matrices(zk, chnkr, pinfo, A)
    wts = chnkr.wts(:);
    Ck = kernel('helm', 'c', zk, [1.0,-1j*zk]);
    Ckp = kernel('helm', 'cprime', zk, [1.0,-1j*zk]);
    Sk = kernel('helm', 's', zk);
    Dk = kernel('helm', 'd', zk);

    pr = pinfo.r;
    pw = pinfo.wts;
    
    [~,npxy] = size(pr);
    sqpw = sqrt(pw);
    S = zeros(2*npxy);
    
    
    cinfo = [];
    cinfo.r = chnkr.r(:,:);
    cinfo.n = chnkr.n(:,:);
    
    D_pxy_to_chnkr = Dk.eval(pinfo, cinfo).*(sqpw(:).');
    S_pxy_to_chnkr = Sk.eval(pinfo, cinfo).*(sqpw(:).');


    C_chnkr_to_pxy = sqpw(:).*Ck.eval(cinfo, pinfo).*(wts(:).');
    Cp_chnkr_to_pxy = sqpw(:).*Ckp.eval(cinfo, pinfo).*(wts(:).');

    AtmpD = A \ D_pxy_to_chnkr;
    AtmpS = A \ S_pxy_to_chnkr;

    iind1 = 1:npxy;
    iind2 = (npxy+1):(2*npxy);
    S(iind1, iind1) = -C_chnkr_to_pxy*AtmpD;
    S(iind1, iind2) = C_chnkr_to_pxy*AtmpS;

    S(iind2, iind1) = -Cp_chnkr_to_pxy*AtmpD;
    S(iind2, iind2) = Cp_chnkr_to_pxy*AtmpS;
end

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

        % chnkr = list_obj{ii}.chnkr;
        % pinfo = list_rect{ii};    
        % n = chnkr.npt;
        % A = chunkermat(chnkr, Ck) + 0.5*eye(n);
        % [S, C, Cp] = get_scattering_matrices(zk, chnkr, pinfo, A);
        nS2 = size(S,2);
        aux1 = [aux1;
                S*aux(nS2*(ii-1)+1:nS2*ii) - aux(nS2*(ii-1)+1:nS2*ii)];
        
    end

    out = x + aux1;

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
