%% MAE 598: Project 3 - Topology Optimization
% Aishwarya Ledalla


%% Topology Optimization of a Random Structure
clear;
clc;

top88(10,10,0.5,3,2,10)

%%%% Modified by Max Yi Ren (ASU) %%%%%%%%%%%%%%%%%%%%%%%
%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%

function top88(nelx,nely,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES: 6061-T6 Aluminum
E0 = 68.9; % GPa
Emin = 1e-9;
nu = 0.33;

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];

KE = 1/(1-nu^2)/24*([A11 A12;A12 A11]+nu*[B11 B12;B12 B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1)),64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8)),64*nelx*nely,1);

%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(3,1,-1,3*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+3)],[2*(nelx+1)*(nely+3)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;

for i1 = 1:nelx
    for j1 = 1:nely

        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;

%% START ITERATION
while change > 0.01
    loop = loop + 1;

    %% FE-ANALYSIS
    K = (KE')*(Emin + (xPhys^penal)*(E0 - Emin));
    U(freedofs) = linsolve(K,F(freedofs));

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = sum(( dot(U(:),KE) * U(:) )); % element-wise strain energy
    c = sum( (Emin + (xPhys^penal)*(Emax-Emin))*ce ) ; % total strain energy
    dc = (-penal * (xPhys^(penal-1))*(Emax-Emin))*ce ; % design sensitivity
    dv = ones(nely,nelx);

    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end

    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    [xnew,g] = oc(nelx,nely,x,volfrac,dc,dv,g);

    if ft == 0
        xPhys = x;
    elseif ft == 1
        xPhys = H*x'/Hs;
    end

    change = max(abs(xnew(:)-x(:)));
    x = xnew;

    %% PRINT RESULTS
    fprintf('It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xPhys(:)),change);

    %% PLOT DENSITIES
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end

%% Optimality criterion
    function [xnew, g] = oc(nelx,nely,x,volfrac,dc,dv,g)
    l1 = 0;
    l2 = 1e9;
    move = 0.2;
    % reshape to perform vector operations
    xnew = zeros(nelx,nely);
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0, max(x - move, min(1, min(x + move, x*sqrt(-dc/dv/lmid)))));
        gt = g + sum(dv*(xnew-x));
        if gt > 0
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund, Department of Solid Mechanics, %
% Technical University of Denmark, %
% DK-2800 Lyngby, Denmark. %
% Please sent your comments to: sigmund@fam.dtu.dk %
% %
% The code is intended for educational purposes and theoretical details %
% are discussed in the paper %
% "Efficient topology optimization in MATLAB using 88 lines of code, %
% E. Andreassen, A. Clausen, M. Schevenels, %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010 %
% This version is based on earlier 99-line code %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization, %
% Vol 21, pp. 120--127. %
% %
% The code as well as a postscript version of the paper can be %
% downloaded from the web-site: http://www.topopt.dtu.dk %
% %
% Disclaimer: %
% The authors reserves all rights but do not guaranty that the code is %
% free from errors. Furthermore, we shall not be liable in any event %
% caused by the use of the program. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%