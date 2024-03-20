	close all
	clear
	clc

	fid=fopen('P4119_1Blade.neu','w');


	mm=11;
	nn=26;
	Hub_flag=0;

	D=1;   % xlsread('data.xlsx',10,'D2');
	Z=1;   % xlsread('data.xlsx',10,'B2');
    
    sheet_no=1;
	r_R=xlsread('data.xlsx',sheet_no,'B6:B16');
	c_D=xlsread('data.xlsx',sheet_no,'C6:C16');
	P_D=xlsread('data.xlsx',sheet_no,'D6:D16');
	theta_skew=xlsread('data.xlsx',sheet_no,'E6:E16');
	rake_D=xlsread('data.xlsx',sheet_no,'F6:F16'); % this is total rake (total rake=apparent rake+skew inuced rake)
											% so if it was non zero, then you should not consider iS (skew induced rake)
	t0_D=xlsread('data.xlsx',sheet_no,'G6:G16');
	f0_c=xlsread('data.xlsx',sheet_no,'H6:H16');

	theta_skew=(pi/180).*theta_skew;

	xc=xlsread('data.xlsx',sheet_no,'J6:J31');
	yc0=xlsread('data.xlsx',sheet_no,'K6:K31');
% 	si=xlsread('data.xlsx',10,'L6:L31');
	y_tmax=xlsread('data.xlsx',sheet_no,'L6:L31');

	[MM,II] = max(yc0);

% 	si=atan(si);
	phi=0:2*pi/Z:2*pi;
	C=D.*c_D;
	tmax=D.*(t0_D);
	fmax=C.*(f0_c);
	r=(D/2).*r_R;

	for i=1:mm
		theta_nt(i)=atan(P_D(i)/r_R(i)/pi);
		tangent(i)=tan(theta_nt(i));
	end
	
	ti=tmax*y_tmax';

	
	iG=r.*rake_D;   % 2*R*RAKE_D
	iS=r.*theta_skew.*tangent';

for j=1:mm
    
    %%%
    for k=1:nn
        yc(k)=yc0(k)*fmax(j)/MM;        
    end
    si=first_derivative(xc,yc);
    si=atan(si);
    %%%
    
    for k=1:nn
        for m=1:Z
            Xc(j,k)=xc(k)*C(j);
            Yc(j,k)=yc(k);%*C(j);
            
            xu(j,k)=Xc(j,k)-ti(j,k)*sin(si(k));
            yu(j,k)=Yc(j,k)+ti(j,k)*cos(si(k));
            
            xl(j,k)=Xc(j,k)+ti(j,k)*sin(si(k));
            yl(j,k)=Yc(j,k)-ti(j,k)*cos(si(k));

            xpu(j,k)=-(iG(j)+iS(j))+(0.5*C(j)-xu(j,k))*sin(theta_nt(j))+yu(j,k)*cos(theta_nt(j));
            ypu(j,k)=r(j)*sin(theta_skew(j)-(1/r(j))*((0.5*C(j)-xu(j,k))*cos(theta_nt(j))-yu(j,k)*sin(theta_nt(j))));
            zpu(j,k)=r(j)*cos(theta_skew(j)-(1/r(j))*((0.5*C(j)-xu(j,k))*cos(theta_nt(j))-yu(j,k)*sin(theta_nt(j))));

            xpl(j,k)=-(iG(j)+iS(j))+(0.5*C(j)-xl(j,k))*sin(theta_nt(j))+yl(j,k)*cos(theta_nt(j));
            ypl(j,k)=r(j)*sin(theta_skew(j)-(1/r(j))*((0.5*C(j)-xl(j,k))*cos(theta_nt(j))-yl(j,k)*sin(theta_nt(j))));
            zpl(j,k)=r(j)*cos(theta_skew(j)-(1/r(j))*((0.5*C(j)-xl(j,k))*cos(theta_nt(j))-yl(j,k)*sin(theta_nt(j))));

            xpu_prime(j,k,m)=xpu(j,k);
            ypu_prime(j,k,m)=ypu(j,k)*cos(phi(m))-zpu(j,k)*sin(phi(m));
            zpu_prime(j,k,m)=ypu(j,k)*sin(phi(m))+zpu(j,k)*cos(phi(m));
            xpl_prime(j,k,m)=xpl(j,k);
            ypl_prime(j,k,m)=ypl(j,k)*cos(phi(m))-zpl(j,k)*sin(phi(m));
            zpl_prime(j,k,m)=ypl(j,k)*sin(phi(m))+zpl(j,k)*cos(phi(m));
        end
    end
end


for l=1:Z
    mesh(xpu_prime(:,:,l),ypu_prime(:,:,l),zpu_prime(:,:,l))
    hold on
    mesh(xpl_prime(:,:,l),ypl_prime(:,:,l),zpl_prime(:,:,l))
end

if Hub_flag == 1
    Dhub=D*r_R(1);
    Rhub=Dhub/2;
    Lhub = Dhub;

    xxx = 0:0.1:2;
    [yh0,zh0,xh0] = cylinder(Rhub*sqrt(1-xxx.^2/2^2),25);
    xh0a = -2*Rhub*xh0 - Rhub;
    mesh(xh0a,yh0,zh0);        

    xh0b = +2*Rhub*xh0 + Rhub;
    mesh(xh0b,yh0,zh0);
    
    lll=ones(1,20);
    [yh1,zh1,xh1] = cylinder(Rhub*lll,25); 
    xh1 = Lhub*xh1 - Rhub;           
    mesh(xh1,yh1,zh1);
end
plot3([-D/2,D/2],[0,0],[0,0],'r')
plot3([0,0],[-D/2,D/2],[0,0],'g')
plot3([0,0],[0,0],[-D/2,D/2],'b')
colormap(winter(1))
view([0 0])
axis equal

[l,m,n]=size(xpu_prime);
xnodes=zeros(l,2*(m-1));
ynodes=zeros(l,2*(m-1));
znodes=zeros(l,2*(m-1));

for i=1:l
    for j=1:m
        xnodes(i,j)=xpu_prime(i,j,1);
        ynodes(i,j)=ypu_prime(i,j,1);
        znodes(i,j)=zpu_prime(i,j,1);
    end
    for j=m+1:2*(m-1)
        xnodes(i,j)=xpl_prime(i,2*m-j,1);
        ynodes(i,j)=ypl_prime(i,2*m-j,1);
        znodes(i,j)=zpl_prime(i,2*m-j,1);
    end
end


[r, c]=size(xnodes);

no_of_foils=r; %mm
no_of_nodes_per_foils=c; %2*(nn-1)

no_of_swise_panels=r-1;
no_of_cwise_panels_per_face=c/2;

npoints=no_of_foils*no_of_nodes_per_foils;
npanels=no_of_swise_panels*no_of_nodes_per_foils;

Nodes=zeros(Z*npoints,3);
Elems=zeros(Z*npanels,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for zz=1:Z
    cnt=1;
    inv_rot_mat=[1 0 0;0 cos((zz-1)*2*pi/Z) -sin((zz-1)*2*pi/Z);0 sin((zz-1)*2*pi/Z) cos((zz-1)*2*pi/Z)];
    for i=1:no_of_foils
        for j=1:no_of_nodes_per_foils
            A=inv_rot_mat*[xnodes(i,j) ynodes(i,j) znodes(i,j)]';
            Nodes(cnt+(zz-1)*npoints,1)=A(1);
            Nodes(cnt+(zz-1)*npoints,2)=A(2);
            Nodes(cnt+(zz-1)*npoints,3)=A(3);
            cnt=cnt+1;
        end
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt=1;
for i=1:no_of_swise_panels
    for j=1:no_of_nodes_per_foils
        if mod(cnt,no_of_nodes_per_foils)==0
            Elems(cnt,1)=cnt;
            Elems(cnt,2)=cnt-no_of_nodes_per_foils+1;
            Elems(cnt,3)=cnt+1;
            Elems(cnt,4)=cnt+no_of_nodes_per_foils;
        else
            Elems(cnt,1)=cnt;
            Elems(cnt,2)=cnt+1;
            Elems(cnt,3)=cnt+1+no_of_nodes_per_foils;
            Elems(cnt,4)=cnt+no_of_nodes_per_foils;
        end
        cnt=cnt+1;
    end
end
for zz=1:Z-1
    cnt=1;
    for i=1:no_of_swise_panels
        for j=1:no_of_nodes_per_foils
            Elems(cnt+zz*npanels,1)=Elems(cnt+(zz-1)*npanels,1)+npoints;
            Elems(cnt+zz*npanels,2)=Elems(cnt+(zz-1)*npanels,2)+npoints;
            Elems(cnt+zz*npanels,3)=Elems(cnt+(zz-1)*npanels,3)+npoints;
            Elems(cnt+zz*npanels,4)=Elems(cnt+(zz-1)*npanels,4)+npoints;
            cnt=cnt+1;
        end
    end      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write headers
fprintf(fid,'        CONTROL INFO 2.4.6\n');
fprintf(fid,'** GAMBIT NEUTRAL FILE\n');
fprintf(fid,'PROPELLER\n');
fprintf(fid,'PROGRAM:                Gambit     VERSION:  2.4.6\n');
fprintf(fid,'\n');
fprintf(fid,'     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\n');
fprintf(fid,'%10d%10d         2         %d         2         3\n',Z*npoints,Z*npanels,Z);
fprintf(fid,'ENDOFSECTION\n');

%write nodal coordinates
fprintf(fid,'   NODAL COORDINATES 2.4.6\n');
for i=1:Z*npoints
    fprintf(fid,'%10d %19.11e %19.11e %19.11e\n',i,Nodes(i,1),Nodes(i,2),Nodes(i,3));
end
fprintf(fid,'ENDOFSECTION\n');

%write panel conectivity
fprintf(fid,'      ELEMENTS/CELLS 2.4.6\n');
for i=1:Z*npanels
    fprintf(fid,'%8d  2  4 %8d %8d %8d %8d\n',i,Elems(i,1),Elems(i,2),Elems(i,3),Elems(i,4));
end
fprintf(fid,'ENDOFSECTION\n');

%write element group
back=zeros(1,Z*npanels/2);
fore=zeros(1,Z*npanels/2);
cnt=1;
for zz=1:Z
    for i=1:no_of_swise_panels
        for j=1:no_of_nodes_per_foils/2
            back(1,cnt)=j+(i-1)*no_of_nodes_per_foils + (zz-1)*npanels;
            fore(1,cnt)=j+(i-1)*no_of_nodes_per_foils + (no_of_nodes_per_foils/2) + (zz-1)*npanels;
            cnt=cnt+1;
        end
    end
end
fprintf(fid,'       ELEMENT GROUP 2.4.6\n');
fprintf(fid,'GROUP: %10d ELEMENTS: %10d MATERIAL: %10d NFLAGS: %10d\n',1,Z*npanels/2,2,1);
fprintf(fid,'                            back\n');
fprintf(fid,'       0\n');
cc=mod(Z*npanels/2,10);
cnt=1;
for i=1:(Z*npanels/2-cc)/10
    for j=1:10
        fprintf(fid,'%8d',back(1,cnt));
        cnt=cnt+1;
    end
    fprintf(fid,'\n');
end
for i=1:cc
        fprintf(fid,'%8d',back(1,cnt));
        cnt=cnt+1;
end
fprintf(fid,'\nENDOFSECTION\n');
fprintf(fid,'       ELEMENT GROUP 2.4.6\n');
fprintf(fid,'GROUP: %10d ELEMENTS: %10d MATERIAL: %10d NFLAGS: %10d\n',2,Z*npanels/2,2,1);
fprintf(fid,'                            fore\n');
fprintf(fid,'       0\n');
cnt=1;
for i=1:(Z*npanels/2-cc)/10
    for j=1:10
        fprintf(fid,'%8d',fore(1,cnt));
        cnt=cnt+1;
    end
    fprintf(fid,'\n');
end
for i=1:cc
        fprintf(fid,'%8d',fore(1,cnt));
        cnt=cnt+1;
end
fprintf(fid,'\nENDOFSECTION\n');

%write boundary conditions
for zz=1:Z
    fprintf(fid,' BOUNDARY CONDITIONS 2.4.6\n');
    fprintf(fid,'                             te%d       1%8d       0       6\n',zz,2*no_of_swise_panels);
    cnt=no_of_nodes_per_foils/2 + (zz-1)*npanels;
    for i=1:no_of_swise_panels
        fprintf(fid,'%10d%5d%5d\n',cnt,2,2);
        cnt=cnt+no_of_nodes_per_foils;
    end
    cnt=1+(no_of_nodes_per_foils/2) + (zz-1)*npanels;
    for i=1:no_of_swise_panels
        fprintf(fid,'%10d%5d%5d\n',cnt,2,4);
        cnt=cnt+no_of_nodes_per_foils;
    end
    fprintf(fid,'ENDOFSECTION\n');
end

fclose(fid);



