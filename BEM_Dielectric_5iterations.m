%% ========================================================
%% 3D BEM WITH DIELECTRICS + ITERATIVE REFINEMENT + VISUALIZATION
%% ========================================================

clc;
clear;
close all;

%% ========================================================
%% Load geometry
%% ========================================================
run('/home/kainat/untitled3.m');

%% ========================================================
%% Parameters
%% ========================================================
eps0     = 8.854187817e-12;
eps_r_E1 = 3;
eps_r_E2 = 4;

%% ========================================================
%% Detect conductors
%% ========================================================
conductors = {};
i = 1;
while evalin('base', sprintf('exist(''C%d'',''var'')', i))
    conductors{i} = eval(sprintf('C%d', i));
    i = i + 1;
end
num_conductors = length(conductors);

%% ========================================================
%% Detect dielectrics
%% ========================================================
dielectrics = {};
j = 1;
while evalin('base', sprintf('exist(''D%d'',''var'')', j))
    dielectrics{j} = eval(sprintf('D%d', j));
    j = j + 1;
end
num_dielectrics = length(dielectrics);

%% ========================================================
%% Iteration control
%% ========================================================
max_iter = 5;
eps_conv = 1e-6;

C_prev = [];

%% ========================================================
%% PRINT HEADER
%% ========================================================
fprintf('Detected %d conductors: ', num_conductors);
for k = 1:num_conductors
    fprintf('C%d ', k);
end
fprintf('\n\n');

%% ========================================================
%% MAIN ITERATION LOOP
%% ========================================================
for iter = 1:max_iter

fprintf('\n====================================================\n');
fprintf('               ITERATION %d\n', iter);
fprintf('====================================================\n');

%% ========================================================
%% Combine panels
%% ========================================================
all_panels = [];
for k = 1:num_conductors
    all_panels = [all_panels; conductors{k}];
end
for k = 1:num_dielectrics
    all_panels = [all_panels; dielectrics{k}];
end

Nc = sum(cellfun(@(c) size(c,1), conductors));
Nd = sum(cellfun(@(d) size(d,1), dielectrics));
N  = Nc + Nd;

fprintf('Total panels: %d\n', N);

%% ========================================================
%% Geometry
%% ========================================================
panel_centers = zeros(N,3);
panel_areas   = zeros(N,1);
panel_normals = zeros(N,3);

for p = 1:N
    P1 = all_panels(p,1:3);
    P2 = all_panels(p,4:6);
    P3 = all_panels(p,7:9);
    P4 = all_panels(p,10:12);

    panel_centers(p,:) = (P1+P2+P3+P4)/4;

    v1 = P2-P1; v2 = P3-P1;
    v3 = P4-P1; v4 = P3-P4;

    panel_areas(p) = 0.5*norm(cross(v1,v2)) + 0.5*norm(cross(v3,v4));

    n_vec = cross(v1,v2);
    if norm(n_vec)>0
        panel_normals(p,:) = n_vec/norm(n_vec);
    end
end

%% ========================================================
%% FULL DEBUG OUTPUT ONLY ITERATION 1
%% ========================================================
if iter == 1
    fprintf('\n========== PANEL CENTERS ==========\n');
    disp(panel_centers);

    fprintf('\n========== PANEL NORMALS ==========\n');
    disp(panel_normals);

    fprintf('\n========== DISTANCE MATRIX ==========\n');
    dist_matrix = zeros(N,N);
    for i=1:N
        for j=1:N
            dist_matrix(i,j)=norm(panel_centers(i,:)-panel_centers(j,:));
        end
    end
    disp(dist_matrix);
end

%% ========================================================
%% Distance matrix (always needed)
%% ========================================================
dist_matrix = zeros(N,N);
for i=1:N
    for j=1:N
        dist_matrix(i,j)=norm(panel_centers(i,:)-panel_centers(j,:));
    end
end

%% ========================================================
%% Equivalent size
%% ========================================================
hx_equiv = sqrt(panel_areas);
hy_equiv = sqrt(panel_areas);

%% ========================================================
%% Build matrix
%% ========================================================
Full_Matrix = zeros(N,N);

for i=1:Nc
    for j=1:N
        if i==j
            hx=hx_equiv(i); hy=hy_equiv(i);
            Full_Matrix(i,j)=(2*(hx*log((hy+sqrt(hx^2+hy^2))/hx)+ ...
                                  hy*log((hx+sqrt(hx^2+hy^2))/hy))) ...
                                  /(4*pi*eps0);
        else
            Full_Matrix(i,j)=panel_areas(j)/(4*pi*eps0*dist_matrix(i,j));
        end
    end
end

for i=Nc+1:N
    for j=1:N
        if i==j
            Full_Matrix(i,j)=(eps_r_E1+eps_r_E2)/(2*eps0);
        else
            r=panel_centers(i,:)-panel_centers(j,:);
            R=norm(r);
            if R>0
                Full_Matrix(i,j)=(eps_r_E1-eps_r_E2)*panel_areas(j)* ...
                    dot(r,panel_normals(i,:))/(4*pi*eps0*R^3);
            end
        end
    end
end

%% ========================================================
%% Solve
%% ========================================================
C_matrix=zeros(num_conductors);

cond_start=zeros(1,num_conductors);
cond_end=zeros(1,num_conductors);

idx=1;
for k=1:num_conductors
    np=size(conductors{k},1);
    cond_start(k)=idx;
    cond_end(k)=idx+np-1;
    idx=idx+np;
end

for exc=1:num_conductors

    V=zeros(N,1);
    V(cond_start(exc):cond_end(exc))=1;

    sigma=Full_Matrix\V;

    if iter==1
        fprintf('\n========== SIGMA VECTOR ==========\n');
        disp(sigma);
    end

    Q=sigma.*panel_areas;

    for k=1:num_conductors
        idxr=cond_start(k):cond_end(k);
        C_matrix(k,exc)=sum(Q(idxr));
    end
end

%% ========================================================
%% OUTPUT FORMAT
%% ========================================================
fprintf('\n=== Area Analysis ===\n');
for k=1:num_conductors
    idxr=cond_start(k):cond_end(k);
    fprintf('Conductor %d: %d panels, total area = %.4e m^2\n', ...
        k,length(idxr),sum(panel_areas(idxr)));
end

fprintf('\n=== Physical Layout (3D) ===\n');
for k=1:num_conductors
    idxr=cond_start(k):cond_end(k);
    c=mean(panel_centers(idxr,:),1);
    fprintf('Conductor %d: Center at (%.4e, %.4e, %.4e)\n', ...
        k,c(1),c(2),c(3));
end

fprintf('\nCapacitance matrix (F):\n');
disp(C_matrix);

%% ========================================================
%% CONVERGENCE
%% ========================================================
if iter>1
    ratio=norm(C_matrix(:)-C_prev)/norm(C_matrix(:));
    fprintf('\nAuto-stop check: ratio = %.3e\n',ratio);
else
    ratio=Inf;
end
C_prev=C_matrix(:);

%% ========================================================
%% VISUALIZATION (EVERY ITERATION)
%% ========================================================
figure(1); clf; hold on;

colors={'r','g','b','m','c','y'};

for i=1:num_conductors
    curC=conductors{i};
    for p=1:size(curC,1)
        Verts=reshape(curC(p,:),3,4)';
        fill3(Verts(:,1),Verts(:,2),Verts(:,3),colors{mod(i-1,6)+1}, ...
            'FaceAlpha',0.4,'EdgeColor','k');
    end
end

for d=1:num_dielectrics
    curD=dielectrics{d};
    for p=1:size(curD,1)
        Verts=reshape(curD(p,:),3,4)';
        fill3(Verts(:,1),Verts(:,2),Verts(:,3),[0.6 0.6 0.6], ...
            'FaceAlpha',0.2,'EdgeColor','k');
    end
end

title(['Iteration ',num2str(iter)]);
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on; view(3);

drawnow;

%% ========================================================
%% SUBDIVISION
%% ========================================================
if iter<max_iter

    fprintf('\nSubdividing 3D panels into 4 sub-panels each...\n');

    new_conductors=cell(1,num_conductors);
    new_dielectrics=cell(1,num_dielectrics);

    for c=1:num_conductors
        curC=conductors{c};
        M=size(curC,1);
        newC=zeros(M*4,12);

        for p=1:M
            V1=curC(p,1:3); V2=curC(p,4:6);
            V3=curC(p,7:9); V4=curC(p,10:12);

            M12=(V1+V2)/2; M23=(V2+V3)/2;
            M34=(V3+V4)/2; M41=(V4+V1)/2;
            Fc=(V1+V2+V3+V4)/4;

            newC((p-1)*4+1,:)=[V1 M12 Fc M41];
            newC((p-1)*4+2,:)=[M12 V2 M23 Fc];
            newC((p-1)*4+3,:)=[Fc M23 V3 M34];
            newC((p-1)*4+4,:)=[M41 Fc M34 V4];
        end
        new_conductors{c}=newC;
    end

    for d=1:num_dielectrics
        curD=dielectrics{d};
        M=size(curD,1);
        newD=zeros(M*4,12);

        for p=1:M
            V1=curD(p,1:3); V2=curD(p,4:6);
            V3=curD(p,7:9); V4=curD(p,10:12);

            M12=(V1+V2)/2; M23=(V2+V3)/2;
            M34=(V3+V4)/2; M41=(V4+V1)/2;
            Fc=(V1+V2+V3+V4)/4;

            newD((p-1)*4+1,:)=[V1 M12 Fc M41];
            newD((p-1)*4+2,:)=[M12 V2 M23 Fc];
            newD((p-1)*4+3,:)=[Fc M23 V3 M34];
            newD((p-1)*4+4,:)=[M41 Fc M34 V4];
        end
        new_dielectrics{d}=newD;
    end

    conductors=new_conductors;
    dielectrics=new_dielectrics;
end

end
