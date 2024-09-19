clear
clc

zk = 1.1;

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 16;
narms = 0;
amp = 0.25;
chnkr = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref); 
chnkr = sort(chnkr);
wts = chnkr.wts; wts = wts(:);

% Set up kernel
Ck = kernel('helm', 'c', zk, [1.0, -1j*zk]);
Ckp = kernel('helm', 'cprime', zk, [1.0, -1j*zk]);

Sk = kernel('helm', 's', zk);
Skp = kernel('helm', 'sprime', zk);

Dk = kernel('helm', 'd', zk);
Dkp = kernel('helm', 'dprime', zk);

n = chnkr.npt;
A = chunkermat(chnkr, Ck) + 0.5*eye(n);

% Analytic solution test
ns = 10;
ns_out = ns; 
nt = 1000;
[sources, sources_out, targets] = get_sources(narms, amp, ns, ns_out, nt);
strengths = randn(ns, 1);

err1 = test_analytic_solution(zk, chnkr, sources, targets, strengths, A);
fprintf('error in analytic solution test=%d\n',err1);

% Setup proxy points 
npxy0 = 30;      % Number of proxy points per edge
npxy = 4*npxy0;
opts = [];
opts.iflege = 1;
bsize = 1.25;
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy, opts);


pn = [ptau(2,:); -ptau(1,:)] ./ sqrt(sum(ptau.^2,1));


pr = pr * bsize;
pw = pw * bsize;

% change proxy to circle
ifcircle = 0;
if ifcircle
    tt = 0:2*pi/npxy:2*pi - 2*pi/npxy;
    ct = cos(tt);
    st = sin(tt);

    rr = npxy/2/pi;
    pr = zeros(2, npxy);
    pn = zeros(2, npxy);
    pw = zeros(npxy, 1);
    
    pr(1,:) = rr*ct;
    pr(2,:) = rr*st;
    
    pn(1,:) = ct;
    pn(2,:) = st;
    
    pw(:) = 2*pi/npxy*rr;
end

pinfo = [];
pinfo.r = pr;
pinfo.n = pn;
pinfo.wts = pw;
sqpw = sqrt(pw);

% plot boundary, proxy and sources
figure(1)
clf
plot(chnkr,'kx'), hold on;
plot(pr(1,:),pr(2,:),'b.');
plot(sources_out(1,:),sources_out(2,:),'ro');
axis equal


%% Now test scattering matrix for 2 objects

chnkr1 = chnkr;
shift = [10 0];
chnkr2 = chunkerfunc(@(t) starfish(t, narms, amp, shift), cparams, pref); 
chnkr2 = sort(chnkr2);


A11 = A;
pinfo1 = pinfo;
[S1, Q1, C1, Cp1] = get_scattering_matrices(zk, chnkr1, pinfo1, A11);

pinfo2 = pinfo;
pinfo2.r = pinfo2.r + shift';

A22 = chunkermat(chnkr2, Ck) + 0.5*eye(n);
[S2, Q2, C2, Cp2] = get_scattering_matrices(zk, chnkr2, pinfo2, A22);

figure(2)
clf
plot(chnkr1,'k.'); hold on;
plot(chnkr2, 'r.');
plot(pinfo1.r(1,:),pinfo1.r(2,:),'ko');
plot(pinfo2.r(1,:),pinfo2.r(2,:),'ro');

axis equal

% set up extended matrix
A21 = chunkerkernevalmat(chnkr1, Ck, chnkr2);
A12 = chunkerkernevalmat(chnkr2, Ck, chnkr1);

Afull = [A11, A12; 
         A21, A22];

% Set up plane wave data as sources;
thet = rand*2*pi*0;
cd = cos(thet); sd = sin(thet);

px1 = pinfo1.r(1,:).';
py1 = pinfo1.r(2,:).';

pnx1 = pinfo1.n(1,:).';
pny1 = pinfo1.n(2,:).';
pw1 = pinfo1.wts;
pw1 = pw1(:);
sqpw1 = sqrt(pw1);

uincs1 = exp(1j*zk*(px1*cd + py1*sd));
gxuincs1 = 1j*zk*cd.*uincs1;
gyuincs1 = 1j*zk*sd.*uincs1;

dudnincs1 = pnx1.*gxuincs1 + pny1.*gyuincs1;


px2 = pinfo2.r(1,:).';
py2 = pinfo2.r(2,:).';

pnx2 = pinfo2.n(1,:).';
pny2 = pinfo2.n(2,:).';
pw2 = pinfo2.wts;
pw2 = pw2(:);
sqpw2 = sqrt(pw2);

uincs2 = exp(1j*zk*(px2*cd + py2*sd));
gxuincs2 = 1j * zk * cd .* uincs2;
gyuincs2 = 1j * zk * sd .* uincs2;

dudnincs2 = pnx2 .* gxuincs2 + pny2 .* gyuincs2;

uin_pxy = [uincs1.*sqpw1; 
          dudnincs1.*sqpw1(:); 
          uincs2.*sqpw2(:); 
          dudnincs2.*sqpw2(:)]; 

[mm,~] = size(S1);

Stotal = [S1*Q1', zeros(mm);
          zeros(mm),  S2*Q2'];

% Test Stotal by solving two decoupled problems on the individual
% scatterers

uin1_pxy = [uincs1.*sqpw1; 
            dudnincs1.*sqpw1(:)];
   
uin2_pxy = [uincs2.*sqpw2(:); 
            dudnincs2.*sqpw2(:)];    
   
cj1 = Q1'*uin1_pxy;
dat1 = Q1*cj1;


cj2 = Q2'*uin2_pxy;
dat2 = Q2*cj2;

% Build the translation operators
D_p1_to_p2 = sqpw.' .* Dk.eval(pinfo1, pinfo2) .* sqpw;
S_p1_to_p2 = sqpw.' .* Sk.eval(pinfo1, pinfo2) .* sqpw; 
Dp_p1_to_p2 = sqpw.' .* Dkp.eval(pinfo1, pinfo2) .* sqpw;
Sp_p1_to_p2 = sqpw.' .* Skp.eval(pinfo1, pinfo2) .* sqpw;

D_p2_to_p1 = sqpw.' .* Dk.eval(pinfo2, pinfo1) .* sqpw; 
Dp_p2_to_p1 = sqpw.' .* Dkp.eval(pinfo2, pinfo1) .* sqpw; 
S_p2_to_p1 = sqpw.' .* Sk.eval(pinfo2, pinfo1) .* sqpw;
Sp_p2_to_p1 = sqpw.' .* Skp.eval(pinfo2, pinfo1) .* sqpw;


[nnn, ~] = size(D_p1_to_p2);
zmat = zeros(nnn);
zeye = eye(nnn);
Tmat = [zmat,        zmat,          D_p2_to_p1,  -S_p2_to_p1;
        zmat,        zmat,          Dp_p2_to_p1, -Sp_p2_to_p1;
        D_p1_to_p2,  -S_p1_to_p2,   zmat,        zmat;
        Dp_p1_to_p2, -Sp_p1_to_p2   zmat,        zmat];
    
errs = test_translation_operators(zk, sources, strengths, pinfo1, pinfo2, Tmat);    
fprintf('errors in greens identities and translation operators = %d\n',errs)


% Build the matrix

[nn, ~] = size(Stotal);
Ssolve = eye(nn) +(Stotal-eye(nn))*Tmat;




udata_pxy = Stotal * uin_pxy;

% Test accuracy of udata_pxy
ubdry1 = exp(1j * zk * (chnkr1.r(1,:)*cd + chnkr1.r(2,:)*sd) ).';
sig1 = A11 \ ubdry1;
uex1 = [sqpw1(:).*(C1*(sig1.*chnkr1.wts(:)));
        sqpw1(:).*(Cp1*(sig1.*chnkr1.wts(:)))];

ubdry2 = exp(1j * zk * (chnkr2.r(1,:)*cd + chnkr2.r(2,:)*sd) ).';
sig2 = A22 \ ubdry2;
uex2 = [sqpw2(:) .* (C2*(sig2.*chnkr2.wts(:)));
        sqpw2(:) .* (Cp2*(sig2.*chnkr2.wts(:))) ];
    
uex_data = [uex1; 
            uex2];    
        
% Test 1 step of translation operators
uin_2 = Tmat * udata_pxy;

% Confirm that uin_2 is the same as solving the pde on 
% object1, and then evaluating the field and it's derivative
% on the proxy of object 2 and vice vrsa

cinfo1 = [];
cinfo1.r = chnkr1.r(:,:);
cinfo1.n = chnkr1.n(:,:);
ueval2 = sqpw2(:) .* (Ck.eval(cinfo1, pinfo2)*(sig1.*chnkr1.wts(:)));
dudneval2 = sqpw2(:) .* (Ckp.eval(cinfo1, pinfo2)*(sig1.*chnkr1.wts(:)));


cinfo2 = [];
cinfo2.r = chnkr2.r(:,:);
cinfo2.n = chnkr2.n(:,:);
ueval1 = sqpw1(:) .* (Ck.eval(cinfo2, pinfo1)*(sig2.*chnkr2.wts(:)));
dudneval1 = sqpw1(:) .* (Ckp.eval(cinfo2, pinfo1)*(sig2.*chnkr2.wts(:)));

uin_2_test = [ueval1; 
              dudneval1;
              ueval2;
              dudneval2];
fprintf('Error in evaluting scattered field at other proxy=%d\n',norm(uin_2_test - uin_2));


% Now test that the solution obtained by using this as data
% on the second obstacle and solving the boundary value problem is the same
% as applying the scattering matrix 

ubdry1_new = Ck.eval(cinfo2, cinfo1)*(sig2.*chnkr2.wts(:));
sig1_new = A11 \ ubdry1_new;

uout1 = [sqpw1(:) .* (C1*(sig1_new.*chnkr1.wts(:)));
        sqpw1(:) .* (Cp1*(sig1_new.*chnkr1.wts(:))) ];


ubdry2_new = Ck.eval(cinfo1, cinfo2)*(sig1.*chnkr1.wts(:));
sig2_new = A22 \ ubdry2_new;

uout2 = [sqpw2(:) .* (C2*(sig2_new.*chnkr2.wts(:)));
        sqpw2(:) .* (Cp2*(sig2_new.*chnkr2.wts(:))) ];

uout_tot = [uout1;
            uout2];
        
uout_1bounce = Stotal*uin_2;
fprintf('Error in solve after 1 scattering bounce = %d\n',norm(uout_tot - uout_1bounce));
uout2 = udata_pxy + uout_1bounce;


uout_pxy = Ssolve \ udata_pxy;

[nn, ~] = size(Stotal);

afun = @(x) x + (Stotal-eye(nn))*(Tmat*x);
uout_pxy_gmres = gmres(afun, udata_pxy, [], 1e-14);

% Construct solution via direct computation
chnkrs(1,2) = chunker();
chnkrs(1) = chnkr1;
chnkrs(2) = chnkr2;

chnkrtotal = merge(chnkrs);
x = chnkrtotal.r(1,:).';
y = chnkrtotal.r(2,:).';
ubdry = exp(1j * zk * (x*cd + y*sd) ); 

Afull2 = chunkermat(chnkrtotal, Ck) + 0.5*eye(2*n);
sig = Afull2 \ ubdry;

cinfo_use = [];
cinfo_use.r = chnkrtotal.r(:,:);
cinfo_use.n = chnkrtotal.n(:,:);

pinfo_use = [];
npxy1 = size(px1,1);
npxy2 = size(px2,1);
pinfo_use.r = zeros(2,npxy1+npxy2);
pinfo_use.n = zeros(2,npxy1+npxy2);
pinfo_use.r(1,1:npxy1) = px1;
pinfo_use.r(1,npxy1+1:end) = px2;

pinfo_use.r(2,1:npxy1) = py1;
pinfo_use.r(2,npxy1+1:end) = py2;

pinfo_use.n(1,1:npxy1) = pnx1;
pinfo_use.n(1,npxy1+1:end) = pnx2;

pinfo_use.n(2,1:npxy1) = pny1;
pinfo_use.n(2,npxy1+1:end) = pny2;
sqpw_total = [sqpw1(:); sqpw2(:)];


u_ex = sqpw_total(:).*(Ck.eval(cinfo_use, pinfo_use)*(sig.*chnkrtotal.wts(:)));
dudn_ex = sqpw_total(:).*(Ckp.eval(cinfo_use, pinfo_use)*(sig.*chnkrtotal.wts(:)));


uout_pxy_ex = [u_ex(1:npxy1);
               dudn_ex(1:npxy1);
               u_ex(npxy1+1:end);
               dudn_ex(npxy1+1:end)];
           
err1 = norm(uout_pxy_ex - uout_pxy);
fprintf('Error in final solution = %d\n',err1);

err2 = norm(uout_pxy_ex - udata_pxy);
err3 = norm(uout_pxy_ex - uout2);

%% test solution at arbitrary target

targinfo = [];
targinfo.r = pinfo1.r + shift.'/2;
targinfo.n = pinfo1.n;

i1ind = 1:npxy1;
i2ind = (npxy1+1):(2*npxy1);
i3ind = (2*npxy1+1):(2*npxy1+npxy2);
i4ind = (2*npxy1+npxy2+1):(2*(npxy1+npxy2));
uuse = uout_pxy;
u_ex2 = (Ck.eval(cinfo_use, targinfo)*(sig.*chnkrtotal.wts(:)));
D1targ = Dk.eval(pinfo1, targinfo)*(uuse(i1ind).*sqpw1(:));
S1targ = Sk.eval(pinfo1, targinfo)*(uuse(i2ind).*sqpw1(:));

D2targ = Dk.eval(pinfo2, targinfo)*(uuse(i3ind).*sqpw2(:));
S2targ = Sk.eval(pinfo2, targinfo)*(uuse(i4ind).*sqpw2(:));

utest2 = D1targ - S1targ + D2targ - S2targ;
fprintf('error in arbitrary point=%d\n',norm(utest2 - u_ex2))

figure(3)
clf
plot(pinfo1.r(1,:), pinfo1.r(2,:),'k.'); hold on;
plot(pinfo2.r(1,:), pinfo2.r(2,:),'k.');
plot(targinfo.r(1,:), targinfo.r(2,:),'b.');
axis equal






function errs = test_translation_operators(zk, sources, strengths, pinfo1, pinfo2, Tmat)
% This subrotuine tests blocks of the outgoing to incoming translation operator 
    Sk = kernel('helm', 's', zk);
    Skp = kernel('helm', 'sprime', zk);

    pw1 = pinfo1.wts;
    pw2 = pinfo2.wts;
    errs = zeros(2,3);
    sqpw1 = sqrt(pw1);
    sqpw2 = sqrt(pw2);
    
    [~, n1] = size(pinfo1.r);
    [~, n2] = size(pinfo2.r);
    
    
    i1 = 1;
    i2 = n1+1;
    i3 = 2*n1 + 1;
    i4 = 2*n1 + n2 + 1;
    i5 = 2*n1 + 2*n2 + 1;
    
    D_p2_to_p1 = Tmat(i1:(i2-1),i3:(i4-1));
    S_p2_to_p1 = -Tmat(i1:(i2-1),i4:(i5-1));
    
    Dp_p2_to_p1 = Tmat(i2:(i3-1),i3:(i4-1));
    Sp_p2_to_p1 = -Tmat(i2:(i3-1),i4:(i5-1));
    
    D_p1_to_p2 = Tmat(i3:(i4-1),i1:(i2-1));
    S_p1_to_p2 = -Tmat(i3:(i4-1),i2:(i3-1));

    Dp_p1_to_p2 = Tmat(i4:(i5-1),i1:(i2-1));
    Sp_p1_to_p2 = -Tmat(i4:(i5-1),i2:(i3-1));

    
    % Test green's identity part of Tmat
             
    srcinfo = [];
    srcinfo.r = sources;
    utest1 = (Sk.eval(srcinfo, pinfo1)*strengths).*sqpw1;
    dudntest1 = (Skp.eval(srcinfo, pinfo1)*strengths).*sqpw1;

    uout = D_p1_to_p2*utest1 - S_p1_to_p2*dudntest1;
    dudnout = Dp_p1_to_p2*utest1 - Sp_p1_to_p2*dudntest1;

    uex_out2 = (Sk.eval(srcinfo, pinfo2)*strengths).*sqpw2;
    dudnex_out2 = (Skp.eval(srcinfo, pinfo2)*strengths).*sqpw2;

    err1 = norm(uout - uex_out2)/norm(uex_out2);
    err2 = norm(dudnout - dudnex_out2)/norm(uex_out2);
    
    errs(1,1) = err1;
    errs(2,1) = err2;


    shift = pinfo2.r(:,1) - pinfo1.r(:,1);
    srcinfo = [];
    srcinfo.r = sources + shift;
    utest2 = (Sk.eval(srcinfo, pinfo2)*strengths).*sqpw2;
    dudntest2 = (Skp.eval(srcinfo, pinfo2)*strengths).*sqpw2;

    uout = D_p2_to_p1*utest2 - S_p2_to_p1*dudntest2;
    dudnout = Dp_p2_to_p1*utest2 - Sp_p2_to_p1*dudntest2;

    uex_out1 = (Sk.eval(srcinfo, pinfo1)*strengths).*sqpw1;
    dudnex_out1 = (Skp.eval(srcinfo, pinfo1)*strengths).*sqpw1;

    err1 = norm(uout - uex_out1)/norm(uex_out1);
    err2 = norm(dudnout - dudnex_out1)/norm(uex_out1);
    
    errs(1,2) = err1;
    errs(2,2) = err2;
    
    uin = [utest1; 
           dudntest1;
           utest2;
           dudntest2];
       
    uout = Tmat*uin;
    
    uout_ex = [uex_out1;
               dudnex_out1;
               uex_out2;
               dudnex_out2];

    
    errs(1,3) = norm(uout-uout_ex);
    errs(2,3) = 0;


end



function [S, Q, C_chnkr_to_pxy, Cp_chnkr_to_pxy] = get_scattering_matrices(zk, chnkr, pinfo, A)
    wts = chnkr.wts(:);
    Ck = kernel('helm', 'c', zk, [1.0,-1j*zk]);
    Ckp = kernel('helm', 'cprime', zk, [1.0,-1j*zk]);
    Sk = kernel('helm', 's', zk);
    Dk = kernel('helm', 'd', zk);

    pr = pinfo.r;
    pw = pinfo.wts;
    pn = pinfo.n;
    [~,npxy] = size(pr);

    nthet = npxy+20;
    dhats = 0:2*pi/nthet:2*pi-2*pi/nthet;
    
    cd = cos(dhats); sd = sin(dhats);
    
    uincs = exp(1j * zk * (pr(1,:).'*cd + pr(2,:).'*sd) );
    gxuincs = 1j * zk * cd .* uincs;
    gyuincs = 1j * zk * sd .* uincs;

    sqpw = sqrt(pw);
    dudnincs = pn(1,:)'.*gxuincs + pn(2,:)'.*gyuincs;
    data_mat = [uincs.*sqpw(:); 
                dudnincs.*sqpw(:)];


    [Q,~] = qr(data_mat,0);
    [~,nrhs] = size(Q);

    S = zeros(size(Q));

    cinfo = [];
    cinfo.r = chnkr.r(:,:);
    cinfo.n = chnkr.n(:,:);
    D_pxy_to_chnkr = Dk.eval(pinfo, cinfo);
    S_pxy_to_chnkr = Sk.eval(pinfo, cinfo);


    C_chnkr_to_pxy = Ck.eval(cinfo, pinfo);
    Cp_chnkr_to_pxy = Ckp.eval(cinfo, pinfo);

    for i=1:nrhs
        uin = Q(1:npxy,i).*sqpw(:);
        dudnin = Q(npxy+1:end,i).*sqpw(:); 

        ubdry = -D_pxy_to_chnkr*uin + S_pxy_to_chnkr*dudnin;
        sig = A \ ubdry;

        uout = C_chnkr_to_pxy*(sig.*wts(:));
        dudnout = Cp_chnkr_to_pxy*(sig.*wts(:));

        S(1:npxy,i) = uout.*sqpw(:);
        S(npxy+1:end,i) = dudnout.*sqpw(:);  
    end

end


function [sources, sources_out, targets] = get_sources(narms, amp, ns, ns_out, nt)


    % sources

    
    ts = 2*pi*rand(ns, 1);
    sources = starfish(ts, narms, amp);
    sources = 0.5*sources;
    

    % exterior sources
    ts = 2*pi*rand(ns_out, 1);
    sources_out = starfish(ts, narms, amp);
    sources_out = sources_out .* (3.5 + 3*repmat(rand(1, ns_out), 2, 1));


    % targets
    ts = 2*pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1 + 3*repmat(rand(1, nt), 2, 1));
end


function err1 = test_analytic_solution(zk, chnkr, sources, targets, strengths, A)
    
    Ck = kernel('helm', 'c', zk, [1.0,-1j*zk]);
    Sk = kernel('helm', 's', zk);
    
   
    srcinfo  = []; srcinfo.r = sources;
    targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
    kernmats = Sk.eval(srcinfo, targinfo);
    ubdry_analytic = kernmats*strengths;

    % Solve
    sol_analytic = A \ ubdry_analytic;


    % Compute exact solution
    srcinfo  = []; srcinfo.r  = sources;
    targinfo = []; targinfo.r = targets;
    kernmatstarg = Sk.eval(srcinfo, targinfo);
    utarg = kernmatstarg*strengths;

    % Compute solution using chunkerkerneval
    % evaluate at targets and compare

    opts.usesmooth = false;
    opts.verb = false;
    opts.quadkgparams = {'RelTol', 1e-16, 'AbsTol', 1.0e-16};


    Dsol = chunkerkerneval(chnkr, Ck, sol_analytic, targets, opts);


    wchnkr = chnkr.wts;
    err1 = norm(utarg-Dsol, 'inf') / dot(abs(sol_analytic(:)), wchnkr(:));

end
