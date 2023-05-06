function Multiphase_mfse_MBB(nptx,npty,refine,volfra,E,corlencoe)
%% THE INPUT PARAMETER MEANS
% nptx      : the observation points in x-axis;
% npty      : the observation points in y-axis;
% refine    : the ratio of adjacent observation point distance to the FE mesh size; 
% volfra    : the vector of different materials maximum volume fraction;
% E         : the vector of different materials Young¡¯s modulus;
% corlencoe : the correlation length ratio in material field.
%% PARAMETER PREPARATION
Num_materials = length(volfra); % number of materials
located = 0:1/(Num_materials):1; % material demarcation point
nu = 0.3; ptdist = 2; E = [1e-9,E]; elelen = ptdist/refine;
nelx = refine*nptx; nely = refine*npty;
tolne = nelx*nely; tolnd = (nelx+1)*(nely+1);
tolnpt = nptx*npty; tolvol = tolne*elelen^2;
fprintf([' Number of material-field  points:%10i \n'...
    ' Number of finite elements:%10i\n'],tolnpt,tolne);
%% PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:tolnd,1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,tolne,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], tolne,1);
iK = reshape(kron(edofMat,ones(8,1))', 64*tolne,1);
jK = reshape(kron(edofMat,ones(1,8))', 64*tolne,1);
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = 1/(1-nu^2)/24*([A11 A12; A12' A11]+nu*[B11 B12; B12' B11]);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*(nelx*(nely+1)+1), 1, -100, 2*tolnd, 1); U = zeros(2*tolnd,1);
fixeddofs = [(1:2)*2*(nely+1),2*(nelx)*(nely+1)+1:2:2*(nelx+1)*(nely+1)];
freedofs = setdiff(1:2*tolnd,fixeddofs);
%% MATERIAL FIELD SERIES EXPANSION
[eIntopMat,ptIntopMat,Xe,Ye] = MFSE2D(nptx,npty,ptdist,refine,corlencoe);
beta = 0; penal = 3;
%% INITIALIZE ITERATION
x = (-(1-volfra(1))*ones(1,tolnpt))/ptIntopMat; x = x';
ePhi = eIntopMat'*x;
loop = 0; obj = 0.;
change = 1.; ichange = 1; neig = length(x); n = neig;
xmin = -1000*ones(n,1); xmax = 1000*ones(n,1);
low = xmin; upp = xmax; xold1 = x;  xold2 = x; clf;
Obj = []; Volf = cell(Num_materials,1);
%% START ITERATION
while (change>=0.005 || beta<12) && loop<200
    loop = loop + 1; objold = obj;
    %% MULTI-HEAVISIDE PROJECTION
    [ePhiProj, edproj] = Multi_Heaviside(ePhi, beta, Num_materials);
    %% ELASTIC MODULUS INTERPOLATION
    Esk = zeros(64,tolne); dEdPhiProj = zeros(tolne,1);
    for i = 1 : Num_materials
        [temp_1,~] = find(ePhiProj <= located(i+1));
        [temp_2,~] = find(ePhiProj >=  located(i));
        locations = intersect(temp_1, temp_2);
        temp_ePhiProj = ePhiProj(locations) / (located(i+1)-located(i))....
            - (located(i+1) / (located(i+1)-located(i)) - 1);     
        Esk(:,locations) = KE(:)*(E(i) + (temp_ePhiProj'.^penal)*(E(i+1) - E(i)));
        dEdPhiProj(locations) = -penal * (E(i+1) - E(i))....
            * (temp_ePhiProj').^(penal-1);
    end
    %% FE-ANALYSIS
    sK = reshape(Esk,64*tolne, 1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    obj = F'*U; ce = sum((U(edofMat)*KE).*U(edofMat), 2);
    dcdx = eIntopMat*(dEdPhiProj.*ce.*edproj);
    vol = zeros(Num_materials,1);dvdphi = zeros(Num_materials,tolne);
    for i = 1 : Num_materials
        if i ~= Num_materials
            [temp_1,~] = find(ePhiProj < located(i+1)); [temp_2,~] = find(ePhiProj >=  located(i));
            locations = intersect(temp_1, temp_2);
            temp_ePhiProj = ePhiProj(locations) / (located(i+1)-located(i))....
                - (located(i+1) / (located(i+1)-located(i)) - 1);
            vol(i) = sum(temp_ePhiProj*elelen^2); dvdphi(i,locations) = 1;
            [temp_1,~] = find(ePhiProj < located(i+2)); [temp_2,~] = find(ePhiProj >=  located(i+1));
            locations = intersect(temp_1, temp_2);
            temp_ePhiProj = ePhiProj(locations) / (located(i+2)-located(i+1))....
                - (located(i+2) / (located(i+2)-located(i+1)) - 1);
            vol(i) = vol(i) + sum((1-temp_ePhiProj)*elelen^2); dvdphi(i,locations) = -1;                  
        else
            [temp_1,~] = find(ePhiProj <= located(end)); [temp_2,~] = find(ePhiProj >=  located(i));
            locations = intersect(temp_1, temp_2);
            temp_ePhiProj = ePhiProj(locations) / (located(i+1)-located(i))....
                - (located(i+1) / (located(i+1)-located(i)) - 1);
            vol(i) = sum(temp_ePhiProj*elelen^2); dvdphi(i,locations) = 1;
        end
    end
    voldgdx = zeros(neig,Num_materials);
    for i = 1 : Num_materials
        voldgdx(:,i) = eIntopMat*(edproj(:).*dvdphi(i,:)'*elelen^2);
    end    
    %% UPDATE DESIGN VARIABLES
    m = Num_materials; cc = 10000*ones(m,1);
    d = zeros(m,1); a0 = 1; a = zeros(m,1);
    fval = zeros(m, 1); dfdx = zeros(m, n);
    for i = 1 : Num_materials
        fval(i) = 100*(vol(i)/tolvol-volfra(i));
        dfdx(i,:) = 100*voldgdx(:,i)/tolvol;
    end
    xmax = x + 0.3; xmin = x - 0.3;
    [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(m,n,loop,x,xmin,xmax,xold1,xold2, ...
        obj,dcdx,fval,dfdx,low,upp,a0,a,cc,d);
    xold2 = xold1; xold1 = x; x = xmma;  ePhi = eIntopMat'*x;
    %% TUNE PROJECTION PARAMETER
    change = abs(obj-objold)/obj;
    if change < 0.005 && loop > 30
        ichange = ichange + 1;
    else
        ichange = 1;
    end
    if mod(ichange,3) == 0
        beta = min(beta + 1,12);
    end
    %% PRINT RESULTS
    fprintf([' It.:%5i Obj.:%9.4f vol1.:%6.3f vol2.:%6.3f numdesvars :%5i' ...
        ' beta:%5.1f ch.:%6.3f\n'],...
        loop,obj,vol(1)/tolvol,vol(2)/tolvol,neig,beta,change);
    %% PLOT MATERIALS
    figure(1); displayx = zeros(nely, 2*nelx);
    displayx(:, 1:nelx) = reshape(ePhiProj, nely, nelx); displayx(:, nelx+1:end) = displayx(:, nelx:-1:1);
    colormap(flipud(jet)); clims = [-1 -0];
    imagesc(-displayx,clims); axis equal; axis tight;axis off
    title('Materials distribution'); pause(1e-6);
    %% PLOT ITERATION HISTORY
    Obj = cat(2, Obj, obj);
    for i = 1 : Num_materials
        Volf{i,1} = cat(2, Volf{i,1}, vol(i)/tolvol);
    end
    figure(2), plotConvergence(1:loop, Obj, Volf);
end
%% EXTRACT DIFFERENT MATERIAL DISTRIBUTIONS
plot_material = zeros(Num_materials,tolne);
for i = 1 : Num_materials
    if i ~= Num_materials
        [temp_1,~] = find(ePhiProj <= located(i+1)); [temp_2,~] = find(ePhiProj >= located(i));
        locations = intersect(temp_1, temp_2);
        plot_material(i,locations) = ePhiProj(locations) / (located(i+1)-located(i))....
            - (located(i+1) / (located(i+1)-located(i)) - 1);
        [temp_1,~] = find(ePhiProj <= located(i+2)); [temp_2,~] = find(ePhiProj >= located(i+1));
        locations = intersect(temp_1, temp_2);
        plot_material(i,locations) = plot_material(i,locations)' + ....
            1 - (ePhiProj(locations) / (located(i+2)-located(i+1))....
            - (located(i+2) / (located(i+2)-located(i+1)) - 1));
    else
        [temp_1,~] = find(ePhiProj <= located(i+1)); [temp_2,~] = find(ePhiProj >= located(i));
        locations = intersect(temp_1, temp_2);
        plot_material(i,locations) = ePhiProj(locations) / (located(i+1)-located(i))....
            - (located(i+1) / (located(i+1)-located(i)) - 1);
    end
end
%% PLOT DIFFERENT MATERIAL DISTRIBUTIONS
Xe = reshape(Xe,nely,nelx); Ye = reshape(Ye,size(Xe));
for i = 1 : Num_materials
    figure(i + 2)
    contourf([Xe nptx*ptdist+Xe],[Ye Ye],....
        [reshape(plot_material(i,:),nely,nelx) fliplr(reshape(plot_material(i,:),nely,nelx))],....
        [0.5 0.5]);
    axis equal; axis tight; axis off
    if i == 1
        title('Weak material');
    else
        title('Strong material');
    end
end
end

function [eIntopMat,ptIntopMat,Xe,Ye] = MFSE2D(nptx,npty,ptdist,refine,corlencoe)
corlen = corlencoe*min(nptx,npty)*ptdist;
elelen = ptdist/refine; nelx = refine*nptx; nely = refine*npty;
tolne = nelx*nely;  tolnpt = nptx*npty;
%% BUILD CORRELATION MATRIX
[Xpt, Ypt] = meshgrid((0.5:1:nptx)*ptdist, (npty-0.5:-1:0.5)*ptdist);
Xpt = Xpt(:); Ypt = Ypt(:);
corMat = zeros(tolnpt,tolnpt);
for i = 1:size(corMat,1)
    for j = i+1:size(corMat,2)
        corMat(i,j) = exp(-(((Xpt(j)-Xpt(i))^2+(Ypt(j)-Ypt(i))^2)/corlen^2));
    end
end
corMat = corMat + corMat';
for i = 1:size(corMat, 1)
    corMat(i,i) = 1;
end
%% DO SERIES EXPANSION OF THE MATERIAL FIELD
if size(corMat,1) < 1e4
    [eigfunMat, eigvalMat] = eig(corMat);
else
    [eigfunMat, eigvalMat] = eigs(corMat,1500);
end
eigvalVec = diag(eigvalMat);
[eigvalVec, eigsortind] = sort(eigvalVec, 'descend');
neig = 0; tmpsum = 0.;
while tmpsum<(1-1e-6)*sum(abs(eigvalVec))
    neig = neig+1;
    tmpsum = tmpsum+eigvalVec(neig);
end
EXPANMat = sparse(1:neig, 1:neig, eigvalVec(1:neig).^(-1/2), neig, neig)...
    *eigfunMat(:,eigsortind(1:neig))'; clear eigfunMat;
%% COMPUTE PHI ON ELEMENTS AND MATERIAL-FIELD POINTS
[Xe, Ye] = meshgrid((0.5:1:nelx)*elelen, (nely-0.5:-1:0.5)*elelen);
Xe = Xe(:); Ye = Ye(:);
eIntopMat = zeros(neig, tolne);
grsize = min(round(tolnpt/20), tolne); ngr = ceil(tolne/grsize);
for igr = 1:ngr
    eind = (igr-1)*grsize + 1:min(igr*grsize, tolne);
    Xe_sub = Xe(eind); Ye_sub = Ye(eind);
    eptvals = exp(-(((repmat(Xpt',length(eind),1)-repmat(Xe_sub, 1, tolnpt)).^2 ...
        +(repmat(Ypt',length(eind),1)-repmat(Ye_sub, 1, tolnpt)).^2)/corlen^2))';
    eptvals(abs(eptvals)<1e-9) = 0;
    eIntopMat(:,eind) = EXPANMat*eptvals;
end
ptIntopMat = EXPANMat*corMat';
end

function [ ePhiProj, edproj] = Multi_Heaviside(ePhi, beta, Num_materials)
%% EXTEND HEAVISIDE PROJECTION
located = -1:2/(Num_materials):1;
ePhiProj = zeros(size(ePhi)); edproj = zeros(size(ePhi));
for i = 1:length(ePhi)
    for j = 1 : length(located) - 1
        if ePhi(i,1) > located(j) && ePhi(i,1) <= located(j+1)
            temp_ePhi = ePhi(i,1)*Num_materials + Num_materials+1-2*j;
            [ePhiProj(i,1),edproj(i,1)] = Threshold(temp_ePhi, beta);
            ePhiProj(i,1) = (ePhiProj(i,1) + 1)/2;
            ePhiProj(i,1) = 1/(Num_materials) * ePhiProj(i,1) + (j-1)/Num_materials;
            edproj(i,1) = edproj(i,1);
        elseif ePhi(i,1) < located(1)
            ePhiProj(i,1) = 0;
            edproj(i,1) = 0;
        elseif ePhi(i,1) > located(end)
            ePhiProj(i,1) = 1;
            edproj(i,1) = 0;
        end
    end
end
end

function [ePhiProj, edproj] = Threshold(ePhi, beta)
%% GENERAL PROJECTION FUNCTION
ePhiProj = sign(ePhi).*(1-exp(-beta*abs(ePhi))+abs(ePhi)*exp(-beta));
ePhiProj(ePhi>1) = 1; ePhiProj(ePhi<-1) = -1;
edproj = beta*exp(-beta*abs(ePhi)) + exp(-beta);
edproj(ePhi>1) = 0; edproj(ePhi<-1) = 0;
end

function plotConvergence(iter, values , volfrac)
%% PLOT OBJECTIVE
axes1 = gca; yyaxis(axes1,'left');
plot(iter, values,'black-','LineWidth',1.5);
ylabel('Structural compliance','FontSize',14,'FontName','Times New Roman','Color',[0 0 1]);
set(axes1,'YColor','black','FontSize',14,'FontName','Times New Roman');
%% PLOT VOLUME FRACTION
yyaxis(axes1,'right');
for i = 1 : length(volfrac)
    volfrac_temp = volfrac{i,1};
    if i == 1
    plot(iter,volfrac_temp,'b-.','LineWidth',1.5);hold on;
    else
    plot(iter,volfrac_temp,'r-.','LineWidth',1.5);hold on;
    end   
end
set(axes1,'ylim',[0 0.5],'ytick',0:.1:0.5,'YColor','black','FontSize',14,'FontName','Times New Roman');
ylabel('Volume constraint','FontSize',14,'FontName','Times New Roman');
xlabel('Number of iterations','FontSize',14,'FontName','Times New Roman');
drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is supplementary to the corresponding paper:                    %
% An effective topological representation and dimensional reduction        %
% approach for multi-material topology optimization,                       %
% Jianwen Bao, Zhaoyou Sun, Pai Liu, Yangjun Luo                           %
% Submitted to               , 2022                                        %
%                                                                          %
% This code is based on                                                    %
% Efficient topology optimization in MATLAB using 88 lines of code.        %
% Andreassen, Erik, et al. Struct Multidiscipl Optim (2021)                %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%