void RoboTunaGenerator(struct Grid *Grids, char *FILENAME)
{
    int i,j,k;
    int n=41; ///number of stations, take it odd for better result!!!!
    int m=15; ///number of segments in each station
    int s=15; ///number of stations for Cadeual fin(necessary to be odd,cause always a section in the middle)
    double f=1.1; ///expansion ratio between stations
    double AR=1.5; ///Aspect ratio of Ovals
    double L=1.0; ///characteristic length of fish body
    double H=0.3*L; ///characteristic length of Cadeual fin
    double Dx[n-1];
    double Dzf[s-1];
    double x[n];
    double zf[s];
    double sum=0.0;
    double npow=(double) (n-1.0)/2.0;
    double nflr=(double) (n-2.0)/2.0;
    double r;
    double y[m][n],z[m][n];
    double xNACte,xNACle;
    double scale,centr;
    char path[80];

    int NACAlen=(int) sizeof(NACA16)/sizeof(double)/2;
    double xNACA[s][NACAlen],yNACA[s][NACAlen],zNACA[s][NACAlen];

    x[0]=-0.3*L;
    Dx[0]=(0.5*L*(f-1.0))/(pow(f,npow)-1.0);

    for(i=1;i<(n-1);i++)
    {
        if(i<=(int)nflr)
            Dx[i]=f*Dx[i-1];
        else
            Dx[i]=Dx[n-2-i];
    }
    for(i=1;i<n;i++)
    {
        sum+=Dx[i-1];
        x[i]=x[0]+sum;
    }

    for (i=1;i<(n-1);i++)
        for (j=0;j<m;j++)
        {
            if((x[i]/L)<=0.1)
            {
                r=0.152*L*tanh(6.0*(x[i]/L)+1.8)/AR;
                y[j][i]=r*sin(j*pi/m);
                z[j][i]=AR*r*cos(j*pi/m);
            }
            else
            {
                r=L*(0.075-0.076*tanh(7.0*(x[i]/L)-3.15))/AR;
                y[j][i]=r*sin(j*pi/m);
                z[j][i]=AR*r*cos(j*pi/m);
            }
        }
    /// in 2ta felan ye meghdar alloc beshan ta bad mehdare daghigh ro peida konam!!!
    Grids->NN=(2*m*(n-2)+2) + (NACAlen*s) + (2*NACAlen*15);  ///body nodes + fin nodes
    Grids->NC=(2*m*(n-1)) + (NACAlen*(s-1)) + (2*NACAlen*14);
    Grids->Nodes=(struct Node*) malloc(sizeof(struct Node)*Grids->NN);
    Grids->Cells=(struct Cell*) malloc(sizeof(struct Cell)*Grids->NC);


    Grids->Nodes[0].Pos.X=x[0];
    Grids->Nodes[0].Pos.Y=0.0;
    Grids->Nodes[0].Pos.Z=0.0;
    k=1;
    for (i=1;i<(n-1);i++)
    {
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=y[j][i];
            Grids->Nodes[k].Pos.Z=z[j][i];
            k++;
        }
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=-y[j][i];
            Grids->Nodes[k].Pos.Z=-z[j][i];
            k++;
        }
    }
    Grids->Nodes[k].Pos.X=x[n-1];
    Grids->Nodes[k].Pos.Y=0.0;
    Grids->Nodes[k].Pos.Z=0.0;

    ///triangular cells of head part
    for(i=0;i<(2*m-1);i++)
    {
        Grids->Cells[i].Type=3;
        Grids->Cells[i].NodePerCell=3;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*3);
        Grids->Cells[i].NodeList[0]=0;
        Grids->Cells[i].NodeList[1]=i+1;
        Grids->Cells[i].NodeList[2]=i+2;
    }
    Grids->Cells[2*m-1].Type=3;
    Grids->Cells[2*m-1].NodePerCell=3;
    Grids->Cells[2*m-1].NodeList=(int *)malloc(sizeof(int)*3);
    Grids->Cells[2*m-1].NodeList[0]=0;
    Grids->Cells[2*m-1].NodeList[1]=2*m;
    Grids->Cells[2*m-1].NodeList[2]=1;

    ///quadrilateral cells of mid part
    for(i=1;i<(n-2);i++)
    {
        for(j=(2*i*m);j<(2*(i+1)*m-1);j++)
        {
            Grids->Cells[j].Type=2;
            Grids->Cells[j].NodePerCell=4;
            Grids->Cells[j].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[j].NodeList[0]=j-2*m+1;
            Grids->Cells[j].NodeList[1]=j+1;
            Grids->Cells[j].NodeList[2]=j+2;
            Grids->Cells[j].NodeList[3]=j-2*m+2;
        }
        Grids->Cells[2*(i+1)*m-1].Type=2;
        Grids->Cells[2*(i+1)*m-1].NodePerCell=4;
        Grids->Cells[2*(i+1)*m-1].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[2*(i+1)*m-1].NodeList[0]=2*m*i;
        Grids->Cells[2*(i+1)*m-1].NodeList[1]=2*m*(i+1);
        Grids->Cells[2*(i+1)*m-1].NodeList[2]=2*m*i+1;
        Grids->Cells[2*(i+1)*m-1].NodeList[3]=2*m*(i-1)+1;
    }

    ///triangular cells of aft part
    for(i=2*m*(n-2);i<(2*m*(n-1)-1);i++)
    {
        Grids->Cells[i].Type=3;
        Grids->Cells[i].NodePerCell=3;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*3);
        Grids->Cells[i].NodeList[0]=i-2*m+1;
        Grids->Cells[i].NodeList[1]=2*m*(n-2)+1;
        Grids->Cells[i].NodeList[2]=i-2*m+2;
    }
    Grids->Cells[2*m*(n-1)-1].Type=3;
    Grids->Cells[2*m*(n-1)-1].NodePerCell=3;
    Grids->Cells[2*m*(n-1)-1].NodeList=(int *)malloc(sizeof(int)*3);
    Grids->Cells[2*m*(n-1)-1].NodeList[0]=2*m*(n-2);
    Grids->Cells[2*m*(n-1)-1].NodeList[1]=2*m*(n-2)+1;
    Grids->Cells[2*m*(n-1)-1].NodeList[2]=2*m*(n-3)+1;


    /// Generation of Caedual Fin
    zf[0]=-0.15*L;
    Dzf[0]=(0.5*H*(f-1.0))/(pow(f,(double)((s-1)/2.0))-1.0);
    for(i=1;i<(s-1);i++)
    {
        if(i<(int)(s/2.0))
            Dzf[i]=f*Dzf[i-1];
        else
            Dzf[i]=Dzf[s-2-i];
    }
    sum=0.0;
    for(i=1;i<s;i++)
    {
        sum+=Dzf[i-1];
        zf[i]=zf[0]+sum;
    }

    ///Fill Node Coordinates of Fins
    int nodeCNT=2*m*(n-2)+2;
    for(i=0;i<s;i++)
    {
        xNACte=L*(-40.74*pow(fabs(zf[i]/L),3.0)+9.666*pow(zf[i]/L,2.0)+0.77);
        xNACle=L*(39.543*pow(fabs(zf[i]/L),3.0)-3.685*pow(zf[i]/L,2.0)+0.636*fabs(zf[i]/L)+0.7);
        scale=(xNACte-xNACle)/(2.0*fabs(NACA16[0][0]));
        centr=(xNACte+xNACle)/2.0;
        for(j=0;j<NACAlen;j++)
        {
            xNACA[i][j]=(NACA16[j][0]*scale)+centr;
            yNACA[i][j]=(NACA16[j][1]*scale);
            zNACA[i][j]=zf[i];

            Grids->Nodes[nodeCNT].Pos.X=xNACA[i][j];
            Grids->Nodes[nodeCNT].Pos.Y=yNACA[i][j];
            Grids->Nodes[nodeCNT].Pos.Z=zNACA[i][j];
            nodeCNT++;
        }
    }

    ///Cell Connectivity for Fins
    nodeCNT=2*m*(n-2)+2;
    int cellCNT=2*m*(n-1);
    for(i=0;i<(s-1);i++)
        for(j=0;j<NACAlen;j++)
        {
            Grids->Cells[cellCNT].Type=2;
            Grids->Cells[cellCNT].NodePerCell=4;
            Grids->Cells[cellCNT].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[cellCNT].NodeList[0]=nodeCNT+i*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[1]=(j==(NACAlen-1))?(nodeCNT+i*NACAlen):(nodeCNT+i*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[2]=(j==(NACAlen-1))?(nodeCNT+(i+1)*NACAlen):(nodeCNT+(i+1)*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[3]=nodeCNT+(i+1)*NACAlen+j;
            cellCNT++;
        }

    ///Generation of Ventoral/Dorsal Fin
    double x_finlet[15];
    x_finlet[0]=0.24;
    for(i=1;i<15;i++)
        x_finlet[i]=x_finlet[i-1]+0.025714;

    ///Fill Node Coordinates of dorsal/ventoral fins
    nodeCNT=(2*m*(n-2)+2) + (NACAlen*s);
    double teta=45.*pi/180.;
    double xxx, yyy, zzz;
    for(i=0;i<15;i++)
    {
        scale=0.02/(2.0*fabs(NACA16[15][0]));//0.08
        centr=x_finlet[i]+0.01;
        for(j=0;j<NACAlen;j++)
        {
            xxx=((NACA16[j][0]*scale)+centr)-((NACA16[15][0]*scale)+centr);
            yyy=(NACA16[j][1]*scale);
            zzz=0.0;

            Grids->Nodes[nodeCNT].Pos.X=(xxx*cos(teta) - zzz*sin(teta)) + ((NACA16[15][0]*scale)+centr);
            Grids->Nodes[nodeCNT].Pos.Y=yyy;
            Grids->Nodes[nodeCNT].Pos.Z=(xxx*sin(teta) + zzz*cos(teta)) + L*(0.075-0.076*tanh(7.0*(x_finlet[i]/L)-3.15));

            nodeCNT++;
        }
    }

    for(i=0;i<15;i++)
    {
        scale=0.02/(2.0*fabs(NACA16[15][0]));
        centr=x_finlet[i]+0.01;
        for(j=0;j<NACAlen;j++)
        {
            xxx=((NACA16[j][0]*scale)+centr)-((NACA16[15][0]*scale)+centr);
            yyy=(NACA16[j][1]*scale);
            zzz=0.0;

            Grids->Nodes[nodeCNT].Pos.X=(xxx*cos(-teta) - zzz*sin(-teta)) + ((NACA16[15][0]*scale)+centr);
            Grids->Nodes[nodeCNT].Pos.Y=yyy;
            Grids->Nodes[nodeCNT].Pos.Z=(xxx*sin(-teta) + zzz*cos(-teta)) - L*(0.075-0.076*tanh(7.0*(x_finlet[i]/L)-3.15));

            nodeCNT++;
        }
    }

    ///Cell Connectivity for dorsal/ventoral fins
    nodeCNT=(2*m*(n-2)+2) + (NACAlen*s);
    cellCNT=(2*m*(n-1)) + (NACAlen*(s-1));
    for(i=0;i<14;i++)
        for(j=0;j<NACAlen;j++)
        {
            Grids->Cells[cellCNT].Type=2;
            Grids->Cells[cellCNT].NodePerCell=4;
            Grids->Cells[cellCNT].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[cellCNT].NodeList[0]=nodeCNT+i*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[1]=nodeCNT+(i+1)*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[2]=(j==(NACAlen-1))?(nodeCNT+(i+1)*NACAlen):(nodeCNT+(i+1)*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[3]=(j==(NACAlen-1))?(nodeCNT+i*NACAlen):(nodeCNT+i*NACAlen+j+1);
            cellCNT++;
        }

    nodeCNT=(2*m*(n-2)+2) + (NACAlen*s) + NACAlen*15;
    for(i=0;i<14;i++)
        for(j=0;j<NACAlen;j++)
        {
            Grids->Cells[cellCNT].Type=2;
            Grids->Cells[cellCNT].NodePerCell=4;
            Grids->Cells[cellCNT].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[cellCNT].NodeList[0]=nodeCNT+i*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[1]=(j==(NACAlen-1))?(nodeCNT+i*NACAlen):(nodeCNT+i*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[2]=(j==(NACAlen-1))?(nodeCNT+(i+1)*NACAlen):(nodeCNT+(i+1)*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[3]=nodeCNT+(i+1)*NACAlen+j;
            cellCNT++;
        }






    FILE *fid;
    sprintf(path,"./mesh/%s.neu",FILENAME);
    fid=fopen(path,"w");
    fprintf(fid,"%26s\n","CONTROL INFO 2.4.6");
    fprintf(fid,"%s\n","** GAMBIT NEUTRAL FILE");
    fprintf(fid,"%s\n","FishGeometry");
    fprintf(fid,"%s\n","PROGRAM:                Gambit     VERSION:  2.4.6");
    fprintf(fid,"%s\n"," ");
    fprintf(fid,"%10s%10s%10s%10s%10s%10s\n","NUMNP","NELEM","NGRPS","NBSETS","NDFCD","NDFVL");
    fprintf(fid,"%10d%10d%10d%10d%10d%10d\n",Grids->NN,Grids->NC,0,0,3,3);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","NODAL COORDINATES 2.4.6");

    for(i=0;i<Grids->NN;i++)
        fprintf(fid,"%10d%20.11e%20.11e%20.11e\n",i+1,Grids->Nodes[i].Pos.X,Grids->Nodes[i].Pos.Y,Grids->Nodes[i].Pos.Z);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","ELEMENTS/CELLS 2.4.6");
    for(i=0;i<Grids->NC;i++)
    {
        fprintf(fid,"%8d%3d%3d ",i+1,Grids->Cells[i].Type,Grids->Cells[i].NodePerCell);
        for(j=0;j<Grids->Cells[i].NodePerCell;j++)
            fprintf(fid,"%8d",1+Grids->Cells[i].NodeList[j]);
        fprintf(fid,"\n");
    }
    fprintf(fid,"%s\n","ENDOFSECTION");
    fclose(fid);



}







/// NOT COMPLETED ///
void GiantDanioGenerator(struct Grid *Grids, char *FILENAME)
{
    int i,j,k;
    int n=41; ///number of stations, take it odd for better result!!!!
    int m=15; ///number of segments in each station
    int s=15; ///number of stations for Cadeual fin(necessary to be odd,cause always a section in the middle)
    double f=1.1; ///expansion ratio between stations
    double AR=1.0; ///Aspect ratio of Ovals
    double L=1.0; ///characteristic length of fish body
    double H=0.3*L; ///characteristic length of Cadeual fin
    double Dx[n-1];
    double Dzf[s-1];
    double x[n];
    double zf[s];
    double sum=0.0;
    double npow=(double) (n-1.0)/2.0;
    double nflr=(double) (n-2.0)/2.0;
    double r;
    double y[m][n],z[m][n];
    double xNACte,xNACle;
    double scale,centr;
    char path[80];

    int NACAlen=(int) sizeof(NACA16)/sizeof(double)/2;
    double xNACA[s][NACAlen],yNACA[s][NACAlen],zNACA[s][NACAlen];

    x[0]=-0.3*L;
    Dx[0]=(0.5*L*(f-1.0))/(pow(f,npow)-1.0);

    for(i=1;i<(n-1);i++)
    {
        if(i<=(int)nflr)
            Dx[i]=f*Dx[i-1];
        else
            Dx[i]=Dx[n-2-i];
    }
    for(i=1;i<n;i++)
    {
        sum+=Dx[i-1];
        x[i]=x[0]+sum;
    }

    for (i=1;i<(n-1);i++)
    {
        for (j=0;j<m;j++)
        {
            if((x[i]/L)<=0.1)
            {
                r=L*(0.152*tanh(6.0*(x[i]/L)+1.8));
                y[j][i]=r*sin(j*pi/m)/AR;
                z[j][i]=r*cos(j*pi/m);
            }
            else
            {
                r=L*(0.075-0.076*tanh(6.3*(x[i]/L)-3.08));
                y[j][i]=r*sin(j*pi/m)/AR;
                z[j][i]=r*cos(j*pi/m);
            }
        }
    }
    Grids->NN=(2*m*(n-2)+2) + (NACAlen*s) + (2*NACAlen*8);  ///body nodes + fin nodes
    Grids->NC=(2*m*(n-1)) + (NACAlen*(s-1)) + (2*NACAlen*7);
    Grids->Nodes=(struct Node*) malloc(sizeof(struct Node)*Grids->NN);
    Grids->Cells=(struct Cell*) malloc(sizeof(struct Cell)*Grids->NC);


    Grids->Nodes[0].Pos.X=x[0];
    Grids->Nodes[0].Pos.Y=0.0;
    Grids->Nodes[0].Pos.Z=L*0.0975;
    k=1;
    for (i=1;i<(n-1);i++)
    {
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=y[j][i];// + L*sin(j*pi/m)*(0.0975*tanh(-(0.3+x[i]/L)/0.15)+0.0975)/AR;
            Grids->Nodes[k].Pos.Z=z[j][i] + L*(0.0975*tanh(-(0.3+x[i]/L)/0.15)+0.0975);
            k++;
        }
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=-y[j][i];// + L*sin(j*pi/m)*(0.0975*tanh(-(0.3+x[i]/L)/0.15)+0.0975)/AR;
            Grids->Nodes[k].Pos.Z=-z[j][i] + L*(0.0975*tanh(-(0.3+x[i]/L)/0.15)+0.0975);
            k++;
        }
    }
    Grids->Nodes[k].Pos.X=x[n-1];
    Grids->Nodes[k].Pos.Y=0.0;
    Grids->Nodes[k].Pos.Z=0.0;

    ///triangular cells of head part
    for(i=0;i<(2*m-1);i++)
    {
        Grids->Cells[i].Type=3;
        Grids->Cells[i].NodePerCell=3;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*3);
        Grids->Cells[i].NodeList[0]=0;
        Grids->Cells[i].NodeList[1]=i+1;
        Grids->Cells[i].NodeList[2]=i+2;
    }
    Grids->Cells[2*m-1].Type=3;
    Grids->Cells[2*m-1].NodePerCell=3;
    Grids->Cells[2*m-1].NodeList=(int *)malloc(sizeof(int)*3);
    Grids->Cells[2*m-1].NodeList[0]=0;
    Grids->Cells[2*m-1].NodeList[1]=2*m;
    Grids->Cells[2*m-1].NodeList[2]=1;

    ///quadrilateral cells of mid part
    for(i=1;i<(n-2);i++)
    {
        for(j=(2*i*m);j<(2*(i+1)*m-1);j++)
        {
            Grids->Cells[j].Type=2;
            Grids->Cells[j].NodePerCell=4;
            Grids->Cells[j].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[j].NodeList[0]=j-2*m+1;
            Grids->Cells[j].NodeList[1]=j+1;
            Grids->Cells[j].NodeList[2]=j+2;
            Grids->Cells[j].NodeList[3]=j-2*m+2;
        }
        Grids->Cells[2*(i+1)*m-1].Type=2;
        Grids->Cells[2*(i+1)*m-1].NodePerCell=4;
        Grids->Cells[2*(i+1)*m-1].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[2*(i+1)*m-1].NodeList[0]=2*m*i;
        Grids->Cells[2*(i+1)*m-1].NodeList[1]=2*m*(i+1);
        Grids->Cells[2*(i+1)*m-1].NodeList[2]=2*m*i+1;
        Grids->Cells[2*(i+1)*m-1].NodeList[3]=2*m*(i-1)+1;
    }

    ///triangular cells of aft part
    for(i=2*m*(n-2);i<(2*m*(n-1)-1);i++)
    {
        Grids->Cells[i].Type=3;
        Grids->Cells[i].NodePerCell=3;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*3);
        Grids->Cells[i].NodeList[0]=i-2*m+1;
        Grids->Cells[i].NodeList[1]=2*m*(n-2)+1;
        Grids->Cells[i].NodeList[2]=i-2*m+2;
    }
    Grids->Cells[2*m*(n-1)-1].Type=3;
    Grids->Cells[2*m*(n-1)-1].NodePerCell=3;
    Grids->Cells[2*m*(n-1)-1].NodeList=(int *)malloc(sizeof(int)*3);
    Grids->Cells[2*m*(n-1)-1].NodeList[0]=2*m*(n-2);
    Grids->Cells[2*m*(n-1)-1].NodeList[1]=2*m*(n-2)+1;
    Grids->Cells[2*m*(n-1)-1].NodeList[2]=2*m*(n-3)+1;


    /// Generation of Caedual Fin
    zf[0]=-0.15*L;
    Dzf[0]=(0.5*H*(f-1.0))/(pow(f,(double)((s-1)/2.0))-1.0);
    for(i=1;i<(s-1);i++)
    {
        if(i<(int)(s/2.0))
            Dzf[i]=f*Dzf[i-1];
        else
            Dzf[i]=Dzf[s-2-i];
    }
    sum=0.0;
    for(i=1;i<s;i++)
    {
        sum+=Dzf[i-1];
        zf[i]=zf[0]+sum;
    }

    ///Fill Node Coordinates of Fins
    int nodeCNT=2*m*(n-2)+2;
    for(i=0;i<s;i++)
    {
        xNACte=L*(-40.74*pow(fabs(zf[i]/L),3.0)+9.666*pow(zf[i]/L,2.0)-0.15*fabs(zf[i]/L)+0.8075);
        xNACle=L*(39.543*pow(fabs(zf[i]/L),3.0)-3.685*pow(zf[i]/L,2.0)+0.636*fabs(zf[i]/L)+0.7);
        scale=(xNACte-xNACle)/(fabs(NACA16[0][0]-NACA16[15][0]));
        centr=(xNACte+xNACle)/2.0;
        for(j=0;j<NACAlen;j++)
        {
            xNACA[i][j]=(NACA16[j][0]-(NACA16[0][0]+NACA16[15][0])/2.0)*scale+centr;
            yNACA[i][j]=(NACA16[j][1]*scale);
            zNACA[i][j]=zf[i];

            Grids->Nodes[nodeCNT].Pos.X=xNACA[i][j];
            Grids->Nodes[nodeCNT].Pos.Y=yNACA[i][j];
            Grids->Nodes[nodeCNT].Pos.Z=zNACA[i][j];
            nodeCNT++;
        }
    }

    ///Cell Connectivity for Fins
    nodeCNT=2*m*(n-2)+2;
    int cellCNT=2*m*(n-1);
    for(i=0;i<(s-1);i++)
        for(j=0;j<NACAlen;j++)
        {
            Grids->Cells[cellCNT].Type=2;
            Grids->Cells[cellCNT].NodePerCell=4;
            Grids->Cells[cellCNT].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[cellCNT].NodeList[0]=nodeCNT+i*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[1]=(j==(NACAlen-1))?(nodeCNT+i*NACAlen):(nodeCNT+i*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[2]=(j==(NACAlen-1))?(nodeCNT+(i+1)*NACAlen):(nodeCNT+(i+1)*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[3]=nodeCNT+(i+1)*NACAlen+j;
            cellCNT++;
        }

    /// ///////////////////////////////
    /// ///////// IN PROGRESS /////////
    /// ///////////////////////////////
    /// ///////////////////////////////
    ///Generation of Ventoral/Dorsal Fin
    double x_finlet_u[8];
    double x_finlet_l[8];
    x_finlet_u[0]=0.3436;
    x_finlet_l[0]=0.4400;
    for(i=1;i<8;i++)
    {
        x_finlet_u[i]=x_finlet_u[i-1]+0.018817; /// ghablan: x_finlet_u[i]=x_finlet_u[i-1]+ 0.018314;
        x_finlet_l[i]=x_finlet_l[i-1]+0.010322; /// ghablan: x_finlet_l[i]=x_finlet_l[i-1]+ 0.012686;
    }

    ///Fill Node Coordinates of dorsal/ventoral fins
    nodeCNT=(2*m*(n-2)+2) + (NACAlen*s);
//    double teta_u[15]={0.0850,0.1323,0.1917,0.2679,0.3679,0.5014,0.6807,0.9154,1.1980,1.4921,1.7523,1.9575,2.1113,2.2258,2.3124};
//    double leng_u[15]={0.1352,0.1212,0.1076,0.0945,0.0821,0.0709,0.0613,0.0543,0.0511,0.0523,0.0576,0.066,0.0766,0.0885,0.1013};
//    double teta_l[15]={0.5676,0.6537,0.7538,0.8689,0.9991,1.1420,1.2934,1.4467,1.5949,1.7323,1.8555,1.9636,2.0570,2.1372,2.2060};
//    double leng_l[15]={0.0859,0.0796,0.074,0.0692,0.0655,0.063,0.0619,0.0622,0.064,0.0671,0.0713,0.0765,0.0824,0.089,0.0961};
    double teta_u[8]={-0.3735,-0.3810,-0.3907,-0.4037,-0.4220,-0.4495,-0.4955,-0.5878};///ok
    double leng_u[8]={0.2706,0.2403,0.2100,0.1797,0.1494,0.1193,0.0892,0.0596};///ok
    double teta_l[8]={-0.3888,-0.4067,-0.4288,-0.4567,-0.4930,-0.5420,-0.6115,-0.7160};///ok
    double leng_l[8]={0.1901,0.1720,0.1540,0.1361,0.1184,0.1008,0.0836,0.0670};///ok

    double xxx, yyy, zzz;
    for(i=0;i<8;i++)
    {
        scale=leng_u[i]/(fabs(NACA10[0][0]-NACA10[15][0]));//0.08
        centr=x_finlet_u[i]+0.5*leng_u[i];
        for(j=0;j<NACAlen;j++)
        {
            xxx=(NACA10[j][0]-NACA10[15][0])*scale;
            yyy=(NACA10[j][1]*scale);
            zzz=0.0;

            Grids->Nodes[nodeCNT].Pos.X=(xxx*cos(teta_u[i]) - zzz*sin(teta_u[i])) + ((NACA10[15][0]*scale)+centr);
            Grids->Nodes[nodeCNT].Pos.Y=yyy;
            Grids->Nodes[nodeCNT].Pos.Z=(xxx*sin(teta_u[i]) + zzz*cos(teta_u[i])) + L*(0.075-0.076*tanh(6.3*(x_finlet_u[0]/L)-3.08)) + 0.0923633*(x_finlet_u[i]-x_finlet_u[0]);

            nodeCNT++;
        }
    }

    for(i=0;i<8;i++)
    {
        scale=leng_l[i]/(fabs(NACA10[0][0]-NACA10[15][0]));
        centr=x_finlet_l[i]+0.5*leng_l[i];
        for(j=0;j<NACAlen;j++)
        {
            xxx=(NACA10[j][0]-NACA10[15][0])*scale;
            yyy=(NACA10[j][1]*scale);
            zzz=0.0;

            Grids->Nodes[nodeCNT].Pos.X=(xxx*cos(-teta_l[i]) - zzz*sin(-teta_l[i])) + ((NACA10[15][0]*scale)+centr);
            Grids->Nodes[nodeCNT].Pos.Y=yyy;
            Grids->Nodes[nodeCNT].Pos.Z=(xxx*sin(-teta_l[i]) + zzz*cos(-teta_l[i])) - L*(0.075-0.076*tanh(6.3*(x_finlet_l[0]/L)-3.08)) - 0.6342763*(x_finlet_l[i]-x_finlet_l[0]);

            nodeCNT++;
        }
    }

    ///Cell Connectivity for dorsal/ventoral fins
    nodeCNT=(2*m*(n-2)+2) + (NACAlen*s);
    cellCNT=(2*m*(n-1)) + (NACAlen*(s-1));
    for(i=0;i<7;i++)
        for(j=0;j<NACAlen;j++)
        {
            Grids->Cells[cellCNT].Type=2;
            Grids->Cells[cellCNT].NodePerCell=4;
            Grids->Cells[cellCNT].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[cellCNT].NodeList[0]=nodeCNT+i*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[1]=(j==(NACAlen-1))?(nodeCNT+i*NACAlen):(nodeCNT+i*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[2]=(j==(NACAlen-1))?(nodeCNT+(i+1)*NACAlen):(nodeCNT+(i+1)*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[3]=nodeCNT+(i+1)*NACAlen+j;
            cellCNT++;
        }

    nodeCNT=(2*m*(n-2)+2) + (NACAlen*s) + NACAlen*8;
    for(i=0;i<7;i++)
        for(j=0;j<NACAlen;j++)
        {
            Grids->Cells[cellCNT].Type=2;
            Grids->Cells[cellCNT].NodePerCell=4;
            Grids->Cells[cellCNT].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[cellCNT].NodeList[0]=nodeCNT+i*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[1]=nodeCNT+(i+1)*NACAlen+j;
            Grids->Cells[cellCNT].NodeList[2]=(j==(NACAlen-1))?(nodeCNT+(i+1)*NACAlen):(nodeCNT+(i+1)*NACAlen+j+1);
            Grids->Cells[cellCNT].NodeList[3]=(j==(NACAlen-1))?(nodeCNT+i*NACAlen):(nodeCNT+i*NACAlen+j+1);
            cellCNT++;
        }






    FILE *fid;
    sprintf(path,"./mesh/%s.neu",FILENAME);
    fid=fopen(path,"w");
    fprintf(fid,"%26s\n","CONTROL INFO 2.4.6");
    fprintf(fid,"%s\n","** GAMBIT NEUTRAL FILE");
    fprintf(fid,"%s\n","FishGeometry");
    fprintf(fid,"%s\n","PROGRAM:                Gambit     VERSION:  2.4.6");
    fprintf(fid,"%s\n"," ");
    fprintf(fid,"%10s%10s%10s%10s%10s%10s\n","NUMNP","NELEM","NGRPS","NBSETS","NDFCD","NDFVL");
    fprintf(fid,"%10d%10d%10d%10d%10d%10d\n",Grids->NN,Grids->NC,0,0,3,3);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","NODAL COORDINATES 2.4.6");

    for(i=0;i<Grids->NN;i++)
        fprintf(fid,"%10d%20.11e%20.11e%20.11e\n",i+1,Grids->Nodes[i].Pos.X,Grids->Nodes[i].Pos.Y,Grids->Nodes[i].Pos.Z);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","ELEMENTS/CELLS 2.4.6");
    for(i=0;i<Grids->NC;i++)
    {
        fprintf(fid,"%8d%3d%3d ",i+1,Grids->Cells[i].Type,Grids->Cells[i].NodePerCell);
        for(j=0;j<Grids->Cells[i].NodePerCell;j++)
            fprintf(fid,"%8d",1+Grids->Cells[i].NodeList[j]);
        fprintf(fid,"\n");
    }
    fprintf(fid,"%s\n","ENDOFSECTION");
    fclose(fid);
}








void eel_Like_Fish_Generator(struct Grid *Grids, int n, int m, double L, char *FILENAME)
{
    int i,j,k;
//    int n=41; ///number of stations, take it odd for better result!!!!
//    int m=15; ///number of segments in each station
//    int s=15; ///number of stations for Cadeual fin(necessary to be odd,cause always a section in the middle)
   // double f=1.0; ///expansion ratio between stations
    double AR; ///Aspect ratio of Ovals
//    double L=1.0; ///characteristic length of fish body
    //double H=0.3*L; ///characteristic length of Cadeual fin
    double Dx[n-1];
    //double Dzf[s-1];
    double x[n];
    //double zf[s];
    double sum=0.0;
//    double npow=(double) (n-1.0)/2.0;
//    double nflr=(double) (n-2.0)/2.0;
    double r;
    double y[m][n],z[m][n];
   // double xNACte,xNACle;
   // double scale,centr;
    char path[80];

    double wh=0.04*L;
    double wt=0.01*L;
    double sb=0.04*L;
    double st=0.95*L;
    double a=0.55*L;
    double b=0.08*L;

    double W, H;

   // int NACAlen=(int) sizeof(NACA16)/sizeof(double)/2;
  //  double xNACA[s][NACAlen],yNACA[s][NACAlen],zNACA[s][NACAlen];

    x[0]=0.0;
//    Dx[0]=L/(double)n;//(0.5*L*(f-1.0))/(pow(f,npow)-1.0);

    for(i=0;i<(n-1);i++)
        Dx[i]=L/(double)n;

    for(i=1;i<n;i++)
    {
        sum+=Dx[i-1];
        x[i]=x[0]+sum;
    }

    for (i=1;i<(n-1);i++)
    {
        for (j=0;j<m;j++)
        {
            H=b*sqrt(1.0-pow((x[i]-a)/a,2.0));

            if( x[i]>=0.0 && x[i]<sb)
                W=sqrt(2.*wh*x[i]-x[i]*x[i]);

            else if( x[i]>=sb && x[i]<st)
                W=wt-(wt-wh)*pow((x[i]-st)/(sb-st),2.0);

            else
                W=wt*(L-x[i])/(L-st);

            r=sqrt(W*W+H*H);
            AR=H/W;
            y[j][i]=r*sin(j*pi/m)/AR;
            z[j][i]=r*cos(j*pi/m);
        }
    }

    Grids->NN=2*m*(n-2)+ 1 + (m+1);
    Grids->NC=2*m*(n-1);
    Grids->Nodes=(struct Node*) malloc(sizeof(struct Node)*Grids->NN);
    Grids->Cells=(struct Cell*) malloc(sizeof(struct Cell)*Grids->NC);



    Grids->Nodes[0].Pos.X=x[0];
    Grids->Nodes[0].Pos.Y=0.0;
    Grids->Nodes[0].Pos.Z=0.0;

    k=1;
    for (i=1;i<(n-1);i++)
    {
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=y[j][i];
            Grids->Nodes[k].Pos.Z=z[j][i];
            k++;
        }
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=-y[j][i];
            Grids->Nodes[k].Pos.Z=-z[j][i];
            k++;
        }
    }

    for(i=k; i<k+m+1; i++)
    {
        Grids->Nodes[i].Pos.X=x[n-1];
        Grids->Nodes[i].Pos.Y=0.0;
        Grids->Nodes[i].Pos.Z=b*sqrt(1.0-pow((x[n-1]-a)/a,2.0))*cos((i-k)*pi/m);
    }

    ///triangular cells of head part
    for(i=0;i<(2*m-1);i++)
    {
        Grids->Cells[i].Type=3;
        Grids->Cells[i].NodePerCell=3;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*3);
        Grids->Cells[i].NodeList[0]=0;
        Grids->Cells[i].NodeList[1]=i+1;
        Grids->Cells[i].NodeList[2]=i+2;
    }
    Grids->Cells[2*m-1].Type=3;
    Grids->Cells[2*m-1].NodePerCell=3;
    Grids->Cells[2*m-1].NodeList=(int *)malloc(sizeof(int)*3);
    Grids->Cells[2*m-1].NodeList[0]=0;
    Grids->Cells[2*m-1].NodeList[1]=2*m;
    Grids->Cells[2*m-1].NodeList[2]=1;

    ///quadrilateral cells of mid part
    for(i=1;i<(n-2);i++)
    {
        for(j=(2*i*m);j<(2*(i+1)*m-1);j++)
        {
            Grids->Cells[j].Type=2;
            Grids->Cells[j].NodePerCell=4;
            Grids->Cells[j].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[j].NodeList[0]=j-2*m+1;
            Grids->Cells[j].NodeList[1]=j+1;
            Grids->Cells[j].NodeList[2]=j+2;
            Grids->Cells[j].NodeList[3]=j-2*m+2;
        }
        Grids->Cells[2*(i+1)*m-1].Type=2;
        Grids->Cells[2*(i+1)*m-1].NodePerCell=4;
        Grids->Cells[2*(i+1)*m-1].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[2*(i+1)*m-1].NodeList[0]=2*m*i;
        Grids->Cells[2*(i+1)*m-1].NodeList[1]=2*m*(i+1);
        Grids->Cells[2*(i+1)*m-1].NodeList[2]=2*m*i+1;
        Grids->Cells[2*(i+1)*m-1].NodeList[3]=2*m*(i-1)+1;
    }

    ///quad cells of aft part
    int cnt0=0;
    for(i=2*m*(n-2);i<2*m*(n-2)+m;i++)
    {
        Grids->Cells[i].Type=2;
        Grids->Cells[i].NodePerCell=4;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[i].NodeList[0]=i-2*m+1;
        Grids->Cells[i].NodeList[1]=2*m*(n-2)+1+cnt0;
        Grids->Cells[i].NodeList[2]=2*m*(n-2)+1+cnt0+1;
        Grids->Cells[i].NodeList[3]=i-2*m+2;
        cnt0++;
    }

    cnt0=m;
    for(i=2*m*(n-2);i<(2*m*(n-2)+m)-1;i++)
    {
        Grids->Cells[i+m].Type=2;
        Grids->Cells[i+m].NodePerCell=4;
        Grids->Cells[i+m].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[i+m].NodeList[0]=i-m+1;
        Grids->Cells[i+m].NodeList[1]=2*m*(n-2)+cnt0 +1;
        Grids->Cells[i+m].NodeList[2]=2*m*(n-2)+cnt0;
        Grids->Cells[i+m].NodeList[3]=i-m+2;
        cnt0--;
    }

    Grids->Cells[i+m].Type=2;
    Grids->Cells[i+m].NodePerCell=4;
    Grids->Cells[i+m].NodeList=(int *)malloc(sizeof(int)*4);
    Grids->Cells[i+m].NodeList[0]=i-m+1;
    Grids->Cells[i+m].NodeList[1]=2*m*(n-2)+2;
    Grids->Cells[i+m].NodeList[2]=2*m*(n-2)+1;
    Grids->Cells[i+m].NodeList[3]=i-3*m+2;
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    FILE *fid;
    sprintf(path,"./mesh/%s.neu",FILENAME);
    fid=fopen(path,"w");
    fprintf(fid,"%26s\n","CONTROL INFO 2.4.6");
    fprintf(fid,"%s\n","** GAMBIT NEUTRAL FILE");
    fprintf(fid,"%s\n","FishGeometry");
    fprintf(fid,"%s\n","PROGRAM:                Gambit     VERSION:  2.4.6");
    fprintf(fid,"%s\n"," ");
    fprintf(fid,"%10s%10s%10s%10s%10s%10s\n","NUMNP","NELEM","NGRPS","NBSETS","NDFCD","NDFVL");
    fprintf(fid,"%10d%10d%10d%10d%10d%10d\n",Grids->NN,Grids->NC,0,0,3,3);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","NODAL COORDINATES 2.4.6");

    for(i=0;i<Grids->NN;i++)
        fprintf(fid,"%10d%20.11e%20.11e%20.11e\n",i+1,Grids->Nodes[i].Pos.X,Grids->Nodes[i].Pos.Y,Grids->Nodes[i].Pos.Z);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","ELEMENTS/CELLS 2.4.6");
    for(i=0;i<Grids->NC;i++)
    {
        fprintf(fid,"%8d%3d%3d ",i+1,Grids->Cells[i].Type,Grids->Cells[i].NodePerCell);
        for(j=0;j<Grids->Cells[i].NodePerCell;j++)
            fprintf(fid,"%8d",1+Grids->Cells[i].NodeList[j]);
        fprintf(fid,"\n");
    }
    fprintf(fid,"%s\n","ENDOFSECTION");
    fclose(fid);
}








void larval_Fish_Generator(struct Grid *Grids, int n, int m, double L, char *FILENAME)
{
    int i,j,k;
//    int n=41; ///number of stations, take it odd for better result!!!!
//    int m=15; ///number of segments in each station
//    int s=15; ///number of stations for Cadeual fin(necessary to be odd,cause always a section in the middle)
   // double f=1.0; ///expansion ratio between stations
    double AR; ///Aspect ratio of Ovals
//    double L=1.0; ///characteristic length of fish body
    //double H=0.3*L; ///characteristic length of Cadeual fin
    double Dx[n-1];
    //double Dzf[s-1];
    double x[n];
    //double zf[s];
    double sum=0.0;
//    double npow=(double) (n-1.0)/2.0;
//    double nflr=(double) (n-2.0)/2.0;
    double r;
    double y[m][n],z[m][n];
   // double xNACte,xNACle;
   // double scale,centr;
    char path[80];

    double wh=0.0635*L;
    double wt=0.0254*L;
    double sb=0.0862*L;
    double st=0.3448*L;
    double s1=0.284*L;
    double s2=0.844*L;
    double s3=0.957*L;
    double h1=0.072*L;
    double h2=0.041*L;
    double h3=0.071*L;

    double W, H;

   // int NACAlen=(int) sizeof(NACA16)/sizeof(double)/2;
  //  double xNACA[s][NACAlen],yNACA[s][NACAlen],zNACA[s][NACAlen];

    x[0]=0.0;
//    Dx[0]=L/(double)n;//(0.5*L*(f-1.0))/(pow(f,npow)-1.0);

    for(i=0;i<(n-1);i++)
        Dx[i]=L/(double)n;

    for(i=1;i<n;i++)
    {
        sum+=Dx[i-1];
        x[i]=x[0]+sum;
    }

    for (i=1;i<(n-1);i++)
    {
        for (j=0;j<m;j++)
        {
            if (x[i]>=0.0 && x[i]<sb){
                W=wh*sqrt(1.0-pow((sb-x[i])/sb,2.0));
                H=h1*sqrt(1.0-pow((x[i]-s1)/s1,2.0));
            }

            else if (x[i]>=sb && x[i]<s1){
                W=wh+(-2.*(wt-wh)-wt*(st-sb))*pow((x[i]-sb)/(st-sb),3.0)+(3.*(wt-wh)+wt*(st-sb))*pow((x[i]-sb)/(st-sb),2.0);
                H=h1*sqrt(1.0-pow((x[i]-s1)/s1,2.0));
            }

            else if (x[i]>=s1 && x[i]<st){
                W=wh+(-2.*(wt-wh)-wt*(st-sb))*pow((x[i]-sb)/(st-sb),3.0)+(3.*(wt-wh)+wt*(st-sb))*pow((x[i]-sb)/(st-sb),2.0);
                H=h1-2.*(h2-h1)*pow((x[i]-s1)/(s2-s1),3.0)+3.*(h2-h1)*pow((x[i]-s1)/(s2-s1),2.0);

            }

            else if (x[i]>=st && x[i]<s2){
                W=wt-wt*pow((x[i]-st)/(L-st),2.0);
                H=h1-2.*(h2-h1)*pow((x[i]-s1)/(s2-s1),3.0)+3.*(h2-h1)*pow((x[i]-s1)/(s2-s1),2.0);

            }

            else if (x[i]>=s2 && x[i]<s3){
                W=wt-wt*pow((x[i]-st)/(L-st),2.0);
                H=h2-2.*(h3-h2)*pow((x[i]-s2)/(s3-s2),3.0)+3.*(h3-h2)*pow((x[i]-s2)/(s3-s2),2.0);

            }

            else {
                W=wt-wt*pow((x[i]-st)/(L-st),2.0);
                H=h3*sqrt(1.0-pow((x[i]-s3)/s3,2.0));

            }

            r=sqrt(W*W+H*H);
            AR=H/W;
            y[j][i]=r*sin(j*pi/m)/AR;
            z[j][i]=r*cos(j*pi/m);
        }
    }

    Grids->NN=2*m*(n-2)+ 1 + (m+1);
    Grids->NC=2*m*(n-1);
    Grids->Nodes=(struct Node*) malloc(sizeof(struct Node)*Grids->NN);
    Grids->Cells=(struct Cell*) malloc(sizeof(struct Cell)*Grids->NC);



    Grids->Nodes[0].Pos.X=x[0];
    Grids->Nodes[0].Pos.Y=0.0;
    Grids->Nodes[0].Pos.Z=0.0;

    k=1;
    for (i=1;i<(n-1);i++)
    {
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=y[j][i];
            Grids->Nodes[k].Pos.Z=z[j][i];
            k++;
        }
        for (j=0;j<m;j++)
        {
            Grids->Nodes[k].Pos.X=x[i];
            Grids->Nodes[k].Pos.Y=-y[j][i];
            Grids->Nodes[k].Pos.Z=-z[j][i];
            k++;
        }
    }

    for(i=k; i<k+m+1; i++)
    {
        Grids->Nodes[i].Pos.X=x[n-1];
        Grids->Nodes[i].Pos.Y=0.0;
        Grids->Nodes[i].Pos.Z=h3*sqrt(1.0-pow((x[n-1]-s3)/s3,2.0))*cos((i-k)*pi/m);
    }

    ///triangular cells of head part
    for(i=0;i<(2*m-1);i++)
    {
        Grids->Cells[i].Type=3;
        Grids->Cells[i].NodePerCell=3;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*3);
        Grids->Cells[i].NodeList[0]=0;
        Grids->Cells[i].NodeList[1]=i+1;
        Grids->Cells[i].NodeList[2]=i+2;
    }
    Grids->Cells[2*m-1].Type=3;
    Grids->Cells[2*m-1].NodePerCell=3;
    Grids->Cells[2*m-1].NodeList=(int *)malloc(sizeof(int)*3);
    Grids->Cells[2*m-1].NodeList[0]=0;
    Grids->Cells[2*m-1].NodeList[1]=2*m;
    Grids->Cells[2*m-1].NodeList[2]=1;

    ///quadrilateral cells of mid part
    for(i=1;i<(n-2);i++)
    {
        for(j=(2*i*m);j<(2*(i+1)*m-1);j++)
        {
            Grids->Cells[j].Type=2;
            Grids->Cells[j].NodePerCell=4;
            Grids->Cells[j].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[j].NodeList[0]=j-2*m+1;
            Grids->Cells[j].NodeList[1]=j+1;
            Grids->Cells[j].NodeList[2]=j+2;
            Grids->Cells[j].NodeList[3]=j-2*m+2;
        }
        Grids->Cells[2*(i+1)*m-1].Type=2;
        Grids->Cells[2*(i+1)*m-1].NodePerCell=4;
        Grids->Cells[2*(i+1)*m-1].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[2*(i+1)*m-1].NodeList[0]=2*m*i;
        Grids->Cells[2*(i+1)*m-1].NodeList[1]=2*m*(i+1);
        Grids->Cells[2*(i+1)*m-1].NodeList[2]=2*m*i+1;
        Grids->Cells[2*(i+1)*m-1].NodeList[3]=2*m*(i-1)+1;
    }

    ///quad cells of aft part
    int cnt0=0;
    for(i=2*m*(n-2);i<2*m*(n-2)+m;i++)
    {
        Grids->Cells[i].Type=2;
        Grids->Cells[i].NodePerCell=4;
        Grids->Cells[i].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[i].NodeList[0]=i-2*m+1;
        Grids->Cells[i].NodeList[1]=2*m*(n-2)+1+cnt0;
        Grids->Cells[i].NodeList[2]=2*m*(n-2)+1+cnt0+1;
        Grids->Cells[i].NodeList[3]=i-2*m+2;
        cnt0++;
    }

    cnt0=m;
    for(i=2*m*(n-2);i<(2*m*(n-2)+m)-1;i++)
    {
        Grids->Cells[i+m].Type=2;
        Grids->Cells[i+m].NodePerCell=4;
        Grids->Cells[i+m].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[i+m].NodeList[0]=i-m+1;
        Grids->Cells[i+m].NodeList[1]=2*m*(n-2)+cnt0 +1;
        Grids->Cells[i+m].NodeList[2]=2*m*(n-2)+cnt0;
        Grids->Cells[i+m].NodeList[3]=i-m+2;
        cnt0--;
    }

    Grids->Cells[i+m].Type=2;
    Grids->Cells[i+m].NodePerCell=4;
    Grids->Cells[i+m].NodeList=(int *)malloc(sizeof(int)*4);
    Grids->Cells[i+m].NodeList[0]=i-m+1;
    Grids->Cells[i+m].NodeList[1]=2*m*(n-2)+2;
    Grids->Cells[i+m].NodeList[2]=2*m*(n-2)+1;
    Grids->Cells[i+m].NodeList[3]=i-3*m+2;
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    FILE *fid;
    sprintf(path,"./mesh/%s.neu",FILENAME);
    fid=fopen(path,"w");
    fprintf(fid,"%26s\n","CONTROL INFO 2.4.6");
    fprintf(fid,"%s\n","** GAMBIT NEUTRAL FILE");
    fprintf(fid,"%s\n","FishGeometry");
    fprintf(fid,"%s\n","PROGRAM:                Gambit     VERSION:  2.4.6");
    fprintf(fid,"%s\n"," ");
    fprintf(fid,"%10s%10s%10s%10s%10s%10s\n","NUMNP","NELEM","NGRPS","NBSETS","NDFCD","NDFVL");
    fprintf(fid,"%10d%10d%10d%10d%10d%10d\n",Grids->NN,Grids->NC,0,0,3,3);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","NODAL COORDINATES 2.4.6");

    for(i=0;i<Grids->NN;i++)
        fprintf(fid,"%10d%20.11e%20.11e%20.11e\n",i+1,Grids->Nodes[i].Pos.X,Grids->Nodes[i].Pos.Y,Grids->Nodes[i].Pos.Z);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","ELEMENTS/CELLS 2.4.6");
    for(i=0;i<Grids->NC;i++)
    {
        fprintf(fid,"%8d%3d%3d ",i+1,Grids->Cells[i].Type,Grids->Cells[i].NodePerCell);
        for(j=0;j<Grids->Cells[i].NodePerCell;j++)
            fprintf(fid,"%8d",1+Grids->Cells[i].NodeList[j]);
        fprintf(fid,"\n");
    }
    fprintf(fid,"%s\n","ENDOFSECTION");
    fclose(fid);
}







void FODGen(struct Grid *Grids, char *FILENAME, int sec, double radius, double MovX, double pitch_angle, double chord_mltply)
{
    int i,j,k;
    double teta[sec];
    char path[80];
    double tmpX, tmpY;

    int NACAlen=(int) sizeof(NACA6412)/sizeof(double)/2;
    double x[sec][NACAlen],y[sec][NACAlen],z[sec][NACAlen];

    Grids->NN=sec*NACAlen;
    Grids->NC=sec*NACAlen;
    Grids->Nodes=(struct Node*) malloc(sizeof(struct Node)*Grids->NN);
    Grids->Cells=(struct Cell*) malloc(sizeof(struct Cell)*Grids->NC);

    for(i=0; i<NACAlen; i++)
        NACA6412[i][1]+=chord_mltply*NACA6412[i][0]*(1.0-NACA6412[i][0]);

    for(i=0; i<NACAlen; i++)
        NACA6412[i][0]+=MovX;

    for(i=0; i<NACAlen; i++){
        tmpX=NACA6412[i][0]*cos(pitch_angle)-NACA6412[i][1]*sin(pitch_angle);
        tmpY=NACA6412[i][0]*sin(pitch_angle)+NACA6412[i][1]*cos(pitch_angle);
        NACA6412[i][0]=tmpX;
        NACA6412[i][1]=tmpY;
    }

    for(i=0; i<sec; i++)
        teta[i]=2.0*pi*i/sec;

    for(i=0; i<sec; i++)
        for(j=0; j<NACAlen; j++)
        {
            x[i][j]=NACA6412[j][0];
            y[i][j]=(NACA6412[j][1]+radius)*cos(teta[i]);
            z[i][j]=(NACA6412[j][1]+radius)*sin(teta[i]);
        }




    k=0;
    for(i=0; i<sec; i++)
        for(j=0; j<NACAlen; j++)
        {
            Grids->Nodes[k].Pos.X=x[i][j];
            Grids->Nodes[k].Pos.Y=y[i][j];
            Grids->Nodes[k].Pos.Z=z[i][j];

            k++;
        }


    k=0;
    for(i=0;i<sec-1;i++)
        for(j=0;j<NACAlen;j++)
        {
            Grids->Cells[k].Type=2;
            Grids->Cells[k].NodePerCell=4;
            Grids->Cells[k].NodeList=(int *)malloc(sizeof(int)*4);
            Grids->Cells[k].NodeList[0]=i*NACAlen+j;
            Grids->Cells[k].NodeList[1]=(j==(NACAlen-1))?(i*NACAlen):(i*NACAlen+j+1);
            Grids->Cells[k].NodeList[2]=(j==(NACAlen-1))?((i+1)*NACAlen):((i+1)*NACAlen+j+1);
            Grids->Cells[k].NodeList[3]=(i+1)*NACAlen+j;
            k++;
        }

    for(j=0;j<NACAlen;j++)
    {
        Grids->Cells[k].Type=2;
        Grids->Cells[k].NodePerCell=4;
        Grids->Cells[k].NodeList=(int *)malloc(sizeof(int)*4);
        Grids->Cells[k].NodeList[0]=i*NACAlen+j;
        Grids->Cells[k].NodeList[1]=(j==(NACAlen-1))?(i*NACAlen):(i*NACAlen+j+1);
        Grids->Cells[k].NodeList[2]=(j==(NACAlen-1))?(0):(j+1);
        Grids->Cells[k].NodeList[3]=j;
        k++;
    }

    FILE *fid;
    sprintf(path,"./mesh/%s.neu",FILENAME);
    fid=fopen(path,"w");
    fprintf(fid,"%26s\n","CONTROL INFO 2.4.6");
    fprintf(fid,"%s\n","** GAMBIT NEUTRAL FILE");
    fprintf(fid,"%s\n","FishGeometry");
    fprintf(fid,"%s\n","PROGRAM:                Gambit     VERSION:  2.4.6");
    fprintf(fid,"%s\n"," ");
    fprintf(fid,"%10s%10s%10s%10s%10s%10s\n","NUMNP","NELEM","NGRPS","NBSETS","NDFCD","NDFVL");
    fprintf(fid,"%10d%10d%10d%10d%10d%10d\n",Grids->NN,Grids->NC,0,0,3,3);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","NODAL COORDINATES 2.4.6");

    for(i=0;i<Grids->NN;i++)
        fprintf(fid,"%10d%20.11e%20.11e%20.11e\n",i+1,Grids->Nodes[i].Pos.X,Grids->Nodes[i].Pos.Y,Grids->Nodes[i].Pos.Z);
    fprintf(fid,"%s\n","ENDOFSECTION");
    fprintf(fid,"%26s\n","ELEMENTS/CELLS 2.4.6");
    for(i=0;i<Grids->NC;i++)
    {
        fprintf(fid,"%8d%3d%3d ",i+1,Grids->Cells[i].Type,Grids->Cells[i].NodePerCell);
        for(j=0;j<Grids->Cells[i].NodePerCell;j++)
            fprintf(fid,"%8d",1+Grids->Cells[i].NodeList[j]);
        fprintf(fid,"\n");
    }
    fprintf(fid,"%s\n","ENDOFSECTION");
    fclose(fid);



}
