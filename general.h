struct MyVector Sum(struct MyVector A ,struct MyVector B)
{
    struct MyVector res;
    res.X=A.X+B.X;
    res.Y=A.Y+B.Y;
    res.Z=A.Z+B.Z;
    return(res);
}

struct MyVector Sum3(struct MyVector A ,struct MyVector B ,struct MyVector C)
{
    struct MyVector res;
    res.X=A.X+B.X+C.X;
    res.Y=A.Y+B.Y+C.Y;
    res.Z=A.Z+B.Z+C.Z;
    return(res);
}

struct MyVector Minus(struct MyVector A,struct MyVector B)
{
    struct MyVector res;
    res.X=A.X-B.X;
    res.Y=A.Y-B.Y;
    res.Z=A.Z-B.Z;
    return(res);
}

struct MyVector MinusVec(struct MyVector A)
{
    struct MyVector res;
    res.X=-A.X;
    res.Y=-A.Y;
    res.Z=-A.Z;
    return(res);
}

double Dot(struct MyVector a ,struct MyVector b)
{
    double res;
    res=(a.X*b.X)+(a.Y*b.Y)+(a.Z*b.Z);
    return(res);
}

struct MyVector Cross(struct MyVector A ,struct MyVector B)
{
    struct MyVector res;
    res.X=A.Y*B.Z-A.Z*B.Y;
    res.Y=A.Z*B.X-A.X*B.Z;
    res.Z=A.X*B.Y-A.Y*B.X;
    return(res);
}

struct MyVector SVP(double a,struct MyVector B)
{
    struct MyVector res;
    res.X=a*B.X;
    res.Y=a*B.Y;
    res.Z=a*B.Z;
    return(res);
}

double sgn(double r)
{
    if (r>0.0)
        return (1.0);
    else if (r<0.0)
        return (-1.0);
    else
        return (0.0);
}

double VDotV3(double V1[4],double V2[4])
{
    return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2] + V1[3]*V2[3]);
}

double VDotV4(double V1[5],double V2[5])
{
    return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2] + V1[3]*V2[3] + V1[4]*V2[4]);
}

double Mag(struct MyVector A)
{
    double res;
    res=sqrt((A.X*A.X)+(A.Y*A.Y)+(A.Z*A.Z));
    return (res);
}

struct MyVector Normalized(struct MyVector Vec)
{
    struct MyVector res;
    res=SVP(1.0/Mag(Vec),Vec);
    return (res);
}

void InvMatrix(double InvA[3][3],double A[3][3])
{
    double Det;
    Det=(-A[0][0]*A[1][1]*A[2][2]+A[0][0]*A[2][1]*A[1][2]+A[1][0]*A[0][1]*A[2][2]-A[1][0]*A[2][1]*A[0][2]-A[2][0]*A[0][1]*A[1][2]+A[2][0]*A[1][1]*A[0][2]);
    InvA[0][0]=-(A[1][1]*A[2][2]-A[2][1]*A[1][2])/Det;
    InvA[0][1]= (A[0][1]*A[2][2]-A[2][1]*A[0][2])/Det;
    InvA[0][2]=-(A[0][1]*A[1][2]-A[1][1]*A[0][2])/Det;
    InvA[1][0]= (A[1][0]*A[2][2]-A[2][0]*A[1][2])/Det;
    InvA[1][1]=-(A[0][0]*A[2][2]-A[2][0]*A[0][2])/Det;
    InvA[1][2]= (A[0][0]*A[1][2]-A[1][0]*A[0][2])/Det;
    InvA[2][0]=-(A[1][0]*A[2][1]-A[2][0]*A[1][1])/Det;
    InvA[2][1]= (A[0][0]*A[2][1]-A[2][0]*A[0][1])/Det;
    InvA[2][2]=-(A[0][0]*A[1][1]-A[1][0]*A[0][1])/Det;
}

struct MyVector MatrixDotVector(double Matrix[3][3],struct MyVector V)
{
    struct MyVector U;
    U.X=V.X*Matrix[0][0]+V.Y*Matrix[0][1]+V.Z*Matrix[0][2];
    U.Y=V.X*Matrix[1][0]+V.Y*Matrix[1][1]+V.Z*Matrix[1][2];
    U.Z=V.X*Matrix[2][0]+V.Y*Matrix[2][1]+V.Z*Matrix[2][2];
    return U;
}

void CEA(double R[3][3],struct MyVector Or)
{
    double f,t;
    f=Or.X;t=Or.Y;
    R[0][0]=1.0;
    R[0][1]=sin(f)*tan(t);
    R[0][2]=cos(f)*tan(t);
    R[1][0]=0.0;
    R[1][1]=cos(f);
    R[1][2]=-sin(f);
    R[2][0]=0.0;
    R[2][1]=sin(f)/cos(t);
    R[2][2]=cos(f)/cos(t);
}

void DCM(double R[3][3],struct MyVector Or)
{
    double f,t,s;
    f=Or.X;t=Or.Y;s=Or.Z;
    R[0][0]=cos(t)*cos(s);
    R[0][1]=cos(t)*sin(s);
    R[0][2]=-sin(t);
    R[1][0]=sin(f)*sin(t)*cos(s)-cos(f)*sin(s);
    R[1][1]=sin(f)*sin(t)*sin(s)+cos(f)*cos(s);
    R[1][2]=sin(f)*cos(t);
    R[2][0]=cos(f)*sin(t)*cos(s)+sin(f)*sin(s);
    R[2][1]=cos(f)*sin(t)*sin(s)-sin(f)*cos(s);
    R[2][2]=cos(f)*cos(t);
}

struct MyVector getPointKinematicVelocity(struct Grid *G, struct MyVector pnt)
{
    double  L=1.5;
    double  w=2.*pi*0.1;
    double  fi=120.0*pi/180.0;
    double  ds0=0.3;
    double  ep=0.0;
    double  y0=0.0;

    double h=0.87;
    double p=20.*pi/180.;

    struct MyVector omg;
    struct MyVector V_pitch;
    struct MyVector RCG, VEC;
    struct MyVector result;
    double zl;

    omg.X=0.0;
    omg.Y=0.0;
    omg.Z=p*w*cos(w*Time+pi);

    zl=pnt.Z/L;

    if (zl<0.0)
        y0=(L*ds0)*(2.*pow(zl,2.) + (4./3.)*pow(zl,3.) + (1./3.)*pow(zl,4.))*(pow(fabs(zl),ep));
    else
        y0=(L*ds0)*(2.*pow(zl,2.) - (4./3.)*pow(zl,3.) + (1./3.)*pow(zl,4.))*(pow(fabs(zl),ep));

    RCG.X=0.0;
    RCG.Y=h*sin(w*Time+pi/2.);
    RCG.Z=0.0;

    VEC=Minus(pnt, RCG);

    V_pitch=Cross(omg,VEC);

//    result.X=V_pitch.X;
//    result.Y=V_pitch.Y + h*w*cos(w*Time+pi/2.) + y0*w*cos(w*Time+fi);
//    result.Z=V_pitch.Z;
    result.X=0.0;
    result.Y=0.0;
    result.Z=0.0;

    return result;




//
//    int i;
//    double w=2.0*6.28;
//    double a=0.5;
//    double tet0=25.*pi/180.;
//    double fi=90.*pi/180.;
//    struct MyVector result;
//    struct MyVector Rcg;
//    struct MyVector Rvec;
//    struct MyVector OMG;
//    struct MyVector V_pitch;
//    double rcg, omg;
//    double theta;
//    double u, MinX=1.0;
//
//    for(i=0; i<G->NN; i++)
//        if(G->Nodes[i].Pos.X<MinX)
//            MinX=G->Nodes[i].Pos.X;


//    u=pnt.X;//-MinX;
//    theta=atan2(pnt.Y, pnt.Z);
//    rcg=1.0+a*sin(w*Time); /// this phase angle, should be same with phase angle in the last line!!!
//
//    Rcg.X=0.0;
//    Rcg.Y=rcg*sin(theta);
//    Rcg.Z=rcg*cos(theta);
//
//    Rvec=Minus(pnt, Rcg);
//
//    omg=tet0*w*cos(w*Time+fi);
//    if (pnt.X<=Rcg.X)
//        OMG=Normalized(Cross(pnt, Rcg));
//    else
//        OMG=MinusVec(Normalized(Cross(pnt, Rcg)));
//
//    OMG=SVP(omg,OMG);
//
//    V_pitch=Cross(OMG, Rvec);

//    result.X= 0.0;//V_pitch.X;
//    result.Y= 0.1*u*(1.0-u)*w*cos(w*Time-fi);//*sin(theta) + V_pitch.Y + a*w*cos(w*Time)*sin(theta);
//    result.Z= 0.0;//0.2*u*(1.0-u)*w*cos(w*Time-fi)*cos(theta);// + V_pitch.Z + a*w*cos(w*Time)*cos(theta);
//
//    return result;



//    double R[3][3],InvR[3][3];
//    struct MyVector pnt_loc;
//    float c1=0.00236;
//    float c2=-0.114;
//    float  k=2.0*pi/1.675;
//    float  w=13.0;
//    float  a=18.76*pi/180.;
//    float  fi=85.0*pi/180.;
//    struct MyVector result, theta_dot, v_rel, v_lin, R_peduncle, rVec;
//
////    pnt_loc=Minus(pnt,G->DynaVariables.MassCenter);
//
//    R_peduncle=G->Nodes[2101].Pos;
//    rVec=Minus(pnt,R_peduncle);
//    if (pnt.X<R_peduncle.X)
//    {
//        result.X=0.0;
//        result.Y=-w*(c1*pnt.X+c2*pnt.X*pnt.X)*cos(k*pnt.X-w*Time);
//        result.Z=0.0;//
//    }
//    else
//    {
//        theta_dot.X=0.0;theta_dot.Y=0.0;theta_dot.Z=-a*w*cos(k*R_peduncle.X-w*Time-fi);
//        v_rel=Cross(theta_dot, rVec);
//        v_lin.X=0.0;v_lin.Y=-w*(c1*R_peduncle.X+c2*R_peduncle.X*R_peduncle.X)*cos(k*R_peduncle.X-w*Time);v_lin.Z=0.0;
//        result.X=0.0;//v_rel.X + v_lin.X;
//        result.Y=v_rel.Y + v_lin.Y;
//        result.Z=0.0;//v_rel.Z + v_lin.Z;
//    }
//
//
//
//
////    rVec=Minus(pnt,G->DynaVariables.MassCenter);
////    result=Sum(G->DynaVariables.LinVel , Cross(G->DynaVariables.AngVel,rVec));  /// BASEC ON DYNAMICS, MERRIAM
//
////    rVec=Minus(pnt,G->DynaVariables.MassCenter);
////    result=Sum(Sum(G->DynaVariables.LinVel , Cross(omega,Minus(G->DynaVariables.MassCenter,ref_p))),
////               Sum(Cross(omega,rVec) , Cross(G->DynaVariables.AngVel,rVec)));  /// BASEC ON DYNAMICS, MERRIAM
//
//
//    return result;
};

void CopyVectors(double **Src_Vec, double **Des_Vec, int noElem)
{
    int i;
    for(i=0;i<noElem;i++)
        (*Des_Vec)[i]=(*Src_Vec)[i];

}

void GlobalGradient(struct Grid *Grids,double **fi,double **u,double **v,double **w)
{
    int i,j,NgbCell;
    double Area,PhiF,coeff,GradDotN;
    struct MyVector Grad, NormalArea, OldGrad, P2F, N2F, EdgeNormal, EdgeCenter;
    struct MyVector n_Vec, rVec, vBT, vREF;

    for(j=0;j<Grids->NC;j++)
    {
        Grad.X=Grad.Y=Grad.Z=0.0;
        OldGrad.X=(*u)[j];OldGrad.Y=(*v)[j];OldGrad.Z=(*w)[j];

        Area=Mag(Grids->Cells[j].Area);
        n_Vec=Normalized(Grids->Cells[j].Area);
        for (i=0;i<(Grids->Cells[j].NodePerCell);i++)
        {
            NgbCell=Grids->Cells[j].Ngb[i];
            if (NgbCell==-1)
            {
                P2F=Minus(Grids->Cells[j].EdgeCenter[i],Grids->Cells[j].CellCenter);
                PhiF=(*fi)[j]+Dot(OldGrad,P2F);
            }
            else
            {
                P2F=Minus(Grids->Cells[j].EdgeCenter[i],Grids->Cells[j].CellCenter);
                N2F=Minus(Grids->Cells[j].EdgeCenter[i],Grids->Cells[NgbCell].CellCenter);
                coeff=Mag(N2F)/(Mag(P2F)+Mag(N2F));
                PhiF=coeff * (*fi)[j] + (1.0-coeff) * (*fi)[NgbCell];
            }

            EdgeNormal=Grids->Cells[j].EdgeNormal[i];
            Grad.X+=PhiF*EdgeNormal.X/Area;
            Grad.Y+=PhiF*EdgeNormal.Y/Area;
            Grad.Z+=PhiF*EdgeNormal.Z/Area;
        }

//        rVec=Minus(Grids->Cells[j].CellCenter,Grids->DynaVariables.MassCenter);
//        vBT=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        vBT=getPointKinematicVelocity(Grids, Grids->Cells[j].CellCenter);
        vREF=Minus(vBT,inf_Vel);  /// not correct, should be summewd!!!
        Grad=Minus(Grad,vREF);
        Grad=Minus(Grad,SVP(Dot(Grad,n_Vec),n_Vec));

        (*u)[j] =Grad.X;
        (*v)[j] =Grad.Y;
        (*w)[j] =Grad.Z;
    }
}

void LSEGradient_Nodal(struct Grid *Grids,double **phi,double **phi_nodal,double **u,double **v,double **w)
{
    double x2bar,xybar,xzbar,y2bar,yzbar,z2bar;
    double xfeebar,yfeebar,zfeebar;
    double determinan;
    int j,i,npc;
    struct MyVector Center,Area,u_Vec,n_Vec,o_Vec;

    struct MyVector N0,N1,N2,N3,Grad;
    double MMM[3][3];

    double *fi;
    struct MyVector *r;

    for(j=0;j<Grids->NC;j++)
    {
        Grad.X=0.0;Grad.Y=0.0;Grad.Z=0.0;
        x2bar=0.0,y2bar=0.0,xybar=0.0;
        xfeebar=0.0,yfeebar=0.0;
        determinan=0.0;

        Center=Grids->Cells[j].CellCenter;
        Area=Grids->Cells[j].Area;
        u_Vec=Normalized(Minus(Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos,Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos));
        n_Vec=Normalized(Area);
        o_Vec=Cross(n_Vec,u_Vec);
        MMM[0][0]=u_Vec.X;  MMM[0][1]=u_Vec.Y;  MMM[0][2]=u_Vec.Z;
        MMM[1][0]=o_Vec.X;  MMM[1][1]=o_Vec.Y;  MMM[1][2]=o_Vec.Z;
        MMM[2][0]=n_Vec.X;  MMM[2][1]=n_Vec.Y;  MMM[2][2]=n_Vec.Z;

        npc=Grids->Cells[j].NodePerCell;
        fi=(double *) calloc(npc+1,sizeof(double));
        r=(struct MyVector *) calloc(npc+1,sizeof(struct MyVector));

        r[0].X=0.0;r[0].Y=0.0;r[0].Z=0.0;
        fi[0]=(*phi)[j];

        if (npc==3)
        {
            N0=Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos;
            N1=Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos;
            N2=Grids->Nodes[Grids->Cells[j].NodeList[2]].Pos;
            r[1]=MatrixDotVector(MMM,Minus(N0,Center));
            r[2]=MatrixDotVector(MMM,Minus(N1,Center));
            r[3]=MatrixDotVector(MMM,Minus(N2,Center));
            fi[1]=(*phi_nodal)[Grids->Cells[j].NodeList[0]];
            fi[2]=(*phi_nodal)[Grids->Cells[j].NodeList[1]];
            fi[3]=(*phi_nodal)[Grids->Cells[j].NodeList[2]];

        }
        else
        {
            N0=Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos;
            N1=Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos;
            N2=Grids->Nodes[Grids->Cells[j].NodeList[2]].Pos;
            N3=Grids->Nodes[Grids->Cells[j].NodeList[3]].Pos;
            r[1]=MatrixDotVector(MMM,Minus(N0,Center));
            r[2]=MatrixDotVector(MMM,Minus(N1,Center));
            r[3]=MatrixDotVector(MMM,Minus(N2,Center));
            r[4]=MatrixDotVector(MMM,Minus(N3,Center));
            fi[1]=(*phi_nodal)[Grids->Cells[j].NodeList[0]];
            fi[2]=(*phi_nodal)[Grids->Cells[j].NodeList[1]];
            fi[3]=(*phi_nodal)[Grids->Cells[j].NodeList[2]];
            fi[4]=(*phi_nodal)[Grids->Cells[j].NodeList[3]];
        }

        for (i=0;i<npc;i++)
        {
            x2bar+=(r[i+1].X-r[0].X)*(r[i+1].X-r[0].X);
            y2bar+=(r[i+1].Y-r[0].Y)*(r[i+1].Y-r[0].Y);

            xybar+=(r[i+1].X-r[0].X)*(r[i+1].Y-r[0].Y);

            xfeebar+=(r[i+1].X-r[0].X)*(fi[i+1]-fi[0]);
            yfeebar+=(r[i+1].Y-r[0].Y)*(fi[i+1]-fi[0]);
        }
        determinan=(x2bar*y2bar)-(xybar*xybar);
        Grad.X=(xfeebar*y2bar - yfeebar*xybar)/determinan;
        Grad.Y=(yfeebar*x2bar - xfeebar*xybar)/determinan;
        Grad.Z=0.0;

        MMM[0][0]=u_Vec.X;  MMM[0][1]=o_Vec.X;  MMM[0][2]=n_Vec.X;
        MMM[1][0]=u_Vec.Y;  MMM[1][1]=o_Vec.Y;  MMM[1][2]=n_Vec.Y;
        MMM[2][0]=u_Vec.Z;  MMM[2][1]=o_Vec.Z;  MMM[2][2]=n_Vec.Z;

        Grad=MatrixDotVector(MMM,Grad);

        (*u)[j]=Grad.X;
        (*v)[j]=Grad.Y;
        (*w)[j]=Grad.Z;
    }
}

void LSEGradient_Ngb(struct Grid *Grids,double **phi,double **u,double **v,double **w, int step)
{
    double x2bar,xybar,xzbar,y2bar,yzbar,z2bar;
    double xfeebar,yfeebar,zfeebar;
    double determinan;
    int j,i,npc;
    struct MyVector Center,Area,u_Vec,n_Vec,o_Vec;
    struct MyVector rVec,vBT, vREF;

    struct MyVector N0,N1,N2,N3,Grad, old_Grad;
    double MMM[3][3];

    double *fi;
    struct MyVector *r;


    for(j=0;j<Grids->NC;j++)
    {
        old_Grad.X=(*u)[j];
        old_Grad.Y=(*v)[j];
        old_Grad.Z=(*w)[j];


        Grad.X=0.0;Grad.Y=0.0;Grad.Z=0.0;
        x2bar=0.0,y2bar=0.0,xybar=0.0;
        xfeebar=0.0,yfeebar=0.0;
        determinan=0.0;

        Center=Grids->Cells[j].CellCenter;
        Area=Grids->Cells[j].Area;
        u_Vec=Normalized(Minus(Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos,Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos));
        n_Vec=Normalized(Area);
        o_Vec=Cross(n_Vec,u_Vec);
        MMM[0][0]=u_Vec.X;  MMM[0][1]=u_Vec.Y;  MMM[0][2]=u_Vec.Z;
        MMM[1][0]=o_Vec.X;  MMM[1][1]=o_Vec.Y;  MMM[1][2]=o_Vec.Z;
        MMM[2][0]=n_Vec.X;  MMM[2][1]=n_Vec.Y;  MMM[2][2]=n_Vec.Z;

        npc=Grids->Cells[j].NodePerCell;
        fi=(double *) calloc(npc+1,sizeof(double));
        r=(struct MyVector *) calloc(npc+1,sizeof(struct MyVector));

        r[0].X=0.0;r[0].Y=0.0;r[0].Z=0.0;
        fi[0]=(*phi)[j];

        if (npc==3)
        {
            N0=(Grids->Cells[j].Ngb[0]!=-1)?(Grids->Cells[Grids->Cells[j].Ngb[0]].CellCenter):(Grids->Cells[j].EdgeCenter[0]);
            N1=(Grids->Cells[j].Ngb[1]!=-1)?(Grids->Cells[Grids->Cells[j].Ngb[1]].CellCenter):(Grids->Cells[j].EdgeCenter[1]);
            N2=(Grids->Cells[j].Ngb[2]!=-1)?(Grids->Cells[Grids->Cells[j].Ngb[2]].CellCenter):(Grids->Cells[j].EdgeCenter[2]);
            r[1]=MatrixDotVector(MMM,Minus(N0,Center));
            r[2]=MatrixDotVector(MMM,Minus(N1,Center));
            r[3]=MatrixDotVector(MMM,Minus(N2,Center));
            fi[1]=(Grids->Cells[j].Ngb[0]!=-1)?((*phi)[Grids->Cells[j].Ngb[0]]):(*phi)[j] ;
            fi[2]=(Grids->Cells[j].Ngb[1]!=-1)?((*phi)[Grids->Cells[j].Ngb[1]]):(*phi)[j] ;
            fi[3]=(Grids->Cells[j].Ngb[2]!=-1)?((*phi)[Grids->Cells[j].Ngb[2]]):(*phi)[j] ;
        }
        else
        {
            N0=(Grids->Cells[j].Ngb[0]!=-1)?(Grids->Cells[Grids->Cells[j].Ngb[0]].CellCenter):(Grids->Cells[j].EdgeCenter[0]);
            N1=(Grids->Cells[j].Ngb[1]!=-1)?(Grids->Cells[Grids->Cells[j].Ngb[1]].CellCenter):(Grids->Cells[j].EdgeCenter[1]);
            N2=(Grids->Cells[j].Ngb[2]!=-1)?(Grids->Cells[Grids->Cells[j].Ngb[2]].CellCenter):(Grids->Cells[j].EdgeCenter[2]);
            N3=(Grids->Cells[j].Ngb[3]!=-1)?(Grids->Cells[Grids->Cells[j].Ngb[3]].CellCenter):(Grids->Cells[j].EdgeCenter[3]);
            r[1]=MatrixDotVector(MMM,Minus(N0,Center));
            r[2]=MatrixDotVector(MMM,Minus(N1,Center));
            r[3]=MatrixDotVector(MMM,Minus(N2,Center));
            r[4]=MatrixDotVector(MMM,Minus(N3,Center));
            fi[1]=(Grids->Cells[j].Ngb[0]!=-1)?((*phi)[Grids->Cells[j].Ngb[0]]):(*phi)[j] ;
            fi[2]=(Grids->Cells[j].Ngb[1]!=-1)?((*phi)[Grids->Cells[j].Ngb[1]]):(*phi)[j] ;
            fi[3]=(Grids->Cells[j].Ngb[2]!=-1)?((*phi)[Grids->Cells[j].Ngb[2]]):(*phi)[j] ;
            fi[4]=(Grids->Cells[j].Ngb[3]!=-1)?((*phi)[Grids->Cells[j].Ngb[3]]):(*phi)[j] ;
        }

        for (i=0;i<npc;i++)
        {
            x2bar+=(r[i+1].X-r[0].X)*(r[i+1].X-r[0].X);
            y2bar+=(r[i+1].Y-r[0].Y)*(r[i+1].Y-r[0].Y);

            xybar+=(r[i+1].X-r[0].X)*(r[i+1].Y-r[0].Y);

            xfeebar+=(r[i+1].X-r[0].X)*(fi[i+1]-fi[0]);
            yfeebar+=(r[i+1].Y-r[0].Y)*(fi[i+1]-fi[0]);
        }
        determinan=(x2bar*y2bar)-(xybar*xybar);
        Grad.X=(xfeebar*y2bar - yfeebar*xybar)/determinan;
        Grad.Y=(yfeebar*x2bar - xfeebar*xybar)/determinan;
        Grad.Z=0.0;

        MMM[0][0]=u_Vec.X;  MMM[0][1]=o_Vec.X;  MMM[0][2]=n_Vec.X;
        MMM[1][0]=u_Vec.Y;  MMM[1][1]=o_Vec.Y;  MMM[1][2]=n_Vec.Y;
        MMM[2][0]=u_Vec.Z;  MMM[2][1]=o_Vec.Z;  MMM[2][2]=n_Vec.Z;
        Grad=MatrixDotVector(MMM,Grad);

        (*u)[j]=Grad.X;
        (*v)[j]=Grad.Y;
        (*w)[j]=Grad.Z;
    }

}





void gradient_vortexje(struct Grid *Grids,double **phi,double **u,double **v,double **w)
{
    struct MyVector center;
    struct MyVector area;
    struct MyVector u_vect;
    struct MyVector n_vect;
    struct MyVector o_vect;
    struct MyVector grad;
    struct MyVector dummy;
    int i, j, k, l, cnt, no_Ngb;
    double MMM[3][3];


    double wmax,wmin,*w_0,*x_0,*b_0;
	double **a_0,**u_0,**v_0;

    for(j=0;j<Grids->NC;j++)
    {
        grad.X=0.0;grad.Y=0.0;grad.Z=0.0;
        center=Grids->Cells[j].CellCenter;
        area=Grids->Cells[j].Area;
        u_vect=Normalized(Minus(Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos,Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos));
        n_vect=Normalized(area);
        o_vect=Cross(n_vect,u_vect);
        MMM[0][0]=u_vect.X;  MMM[0][1]=u_vect.Y;  MMM[0][2]=u_vect.Z;
        MMM[1][0]=o_vect.X;  MMM[1][1]=o_vect.Y;  MMM[1][2]=o_vect.Z;
        MMM[2][0]=n_vect.X;  MMM[2][1]=n_vect.Y;  MMM[2][2]=n_vect.Z;


        no_Ngb = 0;
        for (k = 0; k < Grids->Cells[j].NodePerCell; k++)
            if (Grids->Cells[j].Ngb[k] >= 0)
                no_Ngb++;


        w_0=dvector(1,3);
        x_0=dvector(1,3);
        b_0=dvector(1,1+no_Ngb);
        a_0=dmatrix(1,1+no_Ngb,1,3);
        u_0=dmatrix(1,1+no_Ngb,1,3);
        v_0=dmatrix(1,3,1,3);


        b_0[1] = (*phi)[j];
        cnt = 0;
        for (k = 0; k < Grids->Cells[j].NodePerCell; k++)
            if (Grids->Cells[j].Ngb[k] >= 0)
            {
                cnt++;
                b_0[1+cnt] = (*phi)[Grids->Cells[j].Ngb[k]];
            }

        a_0[1][1] = 0.0;
        a_0[1][2] = 0.0;
        a_0[1][3] = 1.0;
        cnt = 0;
        for (k = 0; k < Grids->Cells[j].NodePerCell; k++)
        {
            if (Grids->Cells[j].Ngb[k] >= 0)
            {
                cnt++;
                dummy = MatrixDotVector(MMM, Minus(Grids->Cells[Grids->Cells[j].Ngb[k]].CellCenter, center));
                a_0[1+cnt][1] = dummy.X;
                a_0[1+cnt][2] = dummy.Y;
                a_0[1+cnt][3] = 1.0;
            }
        }
        /* copy a into u */
        for (k=1;k<=(1+no_Ngb);k++)
            for (l=1;l<=3;l++)
                u_0[k][l]=a_0[k][l];

        /* decompose matrix a */
        svdcmp(u_0,1+no_Ngb,3,w_0,v_0);


        /* find maximum singular value and zero the "small" singular values */
        wmax=0.0;
        for (k=1;k<=3;k++)
            if (w_0[k] > wmax)
                wmax=w_0[k];
        wmin=wmax*(1.0e-6);
        for (k=1;k<=3;k++)
            if (w_0[k] < wmin)
                w_0[k]=0.0;

        svbksb(u_0,w_0,v_0,1+no_Ngb,3,b_0,x_0);


        grad.X=x_0[1];
        grad.Y=x_0[2];
        grad.Z=x_0[3];

        MMM[0][0]=u_vect.X;  MMM[1][0]=u_vect.Y;  MMM[2][0]=u_vect.Z;
        MMM[0][1]=o_vect.X;  MMM[1][1]=o_vect.Y;  MMM[2][1]=o_vect.Z;
        MMM[0][2]=n_vect.X;  MMM[1][2]=n_vect.Y;  MMM[2][2]=n_vect.Z;
        grad = MatrixDotVector(MMM, grad);

		grad=Minus(grad, SVP(Dot(grad,n_vect),n_vect));
        (*u)[j]=-grad.X;
        (*v)[j]=-grad.Y;
        (*w)[j]=-grad.Z;


        free_dmatrix(v_0,1,3,1,3);
        free_dmatrix(u_0,1,1+no_Ngb,1,3);
        free_dmatrix(a_0,1,1+no_Ngb,1,3);
        free_dvector(b_0,1,1+no_Ngb);
        free_dvector(x_0,1,3);
        free_dvector(w_0,1,3);

    }
}






void SurfaceGradient_QUAD(struct Grid *Grids,double **fi,double **u,double **v,double **w)
{
    int j,Ngb0,Ngb1,Ngb2,Ngb3;
    double fi_n0,fi_n1,fi_n2,fi_n3,rphi_rk,rphi_re, sigma;
    struct MyVector N0,N1,N2,N3,drdk,drde,e_k,e_e;
    struct MyVector PN0,PN1,PN2,PN3,Grad;
    struct MyVector n_Vec,rVec,vBT,vREF;

    for(j=0;j<Grids->NC;j++)
    {
        Grad.X=0.0;
        Grad.Y=0.0;
        Grad.Z=0.0;

        N0=Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos;
        N1=Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos;
        N2=Grids->Nodes[Grids->Cells[j].NodeList[2]].Pos;
        N3=Grids->Nodes[Grids->Cells[j].NodeList[3]].Pos;
        n_Vec=Normalized(Grids->Cells[j].Area);
        drdk=Sum(SVP(0.5,Sum(N1,N2)),SVP(-0.5,Sum(N0,N3)));
        drde=Sum(SVP(0.5,Sum(N2,N3)),SVP(-0.5,Sum(N0,N1)));
        e_k=Normalized(drdk);
        e_e=Normalized(drde);
        Ngb0=Grids->Cells[j].Ngb[0];
        Ngb1=Grids->Cells[j].Ngb[1];
        Ngb2=Grids->Cells[j].Ngb[2];
        Ngb3=Grids->Cells[j].Ngb[3];

        if ((Ngb0!=-1)&&(Ngb2!=-1))
        {
            fi_n0=(*fi)[Ngb0];
            fi_n2=(*fi)[Ngb2];
            PN0=Minus(Grids->Cells[Ngb0].CellCenter,Grids->Cells[j].CellCenter);
            PN2=Minus(Grids->Cells[Ngb2].CellCenter,Grids->Cells[j].CellCenter);
            rphi_re=((((*fi)[j]-fi_n0)/Mag(PN0))*Mag(PN2)+((fi_n2-(*fi)[j])/Mag(PN2))*Mag(PN0))/(Mag(PN0)+Mag(PN2));
        }
        if ((Ngb0!=-1)&&(Ngb2==-1))
        {
            fi_n0=(*fi)[Ngb0];
            PN0=Minus(Grids->Cells[Ngb0].CellCenter,Grids->Cells[j].CellCenter);
            rphi_re=((*fi)[j]-fi_n0)/Mag(PN0);
        }
        if ((Ngb0==-1)&&(Ngb2!=-1))
        {
            fi_n2=(*fi)[Ngb2];
            PN2=Minus(Grids->Cells[Ngb2].CellCenter,Grids->Cells[j].CellCenter);
            rphi_re=(fi_n2-(*fi)[j])/Mag(PN2);
        }
        if ((Ngb0==-1)&&(Ngb2==-1))
        {
            rphi_re=0.0;
        }


        if ((Ngb1!=-1)&&(Ngb3!=-1))
        {
            fi_n1=(*fi)[Ngb1];
            fi_n3=(*fi)[Ngb3];
            PN1=Minus(Grids->Cells[Ngb1].CellCenter,Grids->Cells[j].CellCenter);
            PN3=Minus(Grids->Cells[Ngb3].CellCenter,Grids->Cells[j].CellCenter);
            rphi_rk=((((*fi)[j]-fi_n3)/Mag(PN3))*Mag(PN1)+((fi_n1-(*fi)[j])/Mag(PN1))*Mag(PN3))/(Mag(PN1)+Mag(PN3));
        }
        if ((Ngb1!=-1)&&(Ngb3==-1))
        {
            fi_n1=(*fi)[Ngb1];
            PN1=Minus(Grids->Cells[Ngb1].CellCenter,Grids->Cells[j].CellCenter);
            rphi_rk=(fi_n1-(*fi)[j])/Mag(PN1);
        }
        if ((Ngb1==-1)&&(Ngb3!=-1))
        {
            fi_n3=(*fi)[Ngb3];
            PN3=Minus(Grids->Cells[Ngb3].CellCenter,Grids->Cells[j].CellCenter);
            rphi_rk=((*fi)[j]-fi_n3)/Mag(PN3);
        }
        if ((Ngb1==-1)&&(Ngb3==-1))
        {
            rphi_rk=0.0;
        }

        Grad=Sum(SVP(rphi_rk,e_k),SVP(rphi_re,e_e));

//        rVec=Minus(Grids->Cells[j].CellCenter,Grids->DynaVariables.MassCenter);
//        vBT=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
                vBT=getPointKinematicVelocity(Grids, Grids->Cells[j].CellCenter);
        vREF=Minus(vBT,inf_Vel);

        sigma=Dot(n_Vec, vREF);

        Grad=Sum(Grad,SVP(sigma,n_Vec));

        (*u)[j]=Grad.X;
        (*v)[j]=Grad.Y;
        (*w)[j]=Grad.Z;
    }
}

void TangentialGradient_Nodal(struct Grid *Grids,double **fi,double**fi_nodal,double **u,double **v,double **w)
{
    int i,j,npc;
    double fi_n0,fi_n1,fi_n2,fi_n3;
    double a,b,c,d,e,f;
    double r1,r2;
    struct MyVector PC,N0,N1,N2,N3;
    struct MyVector t1,t2,n,Grad;
    struct MyVector rVec,vBT,vREF;

    for(i=0;i<Grids->NC;i++)
    {
        Grad.X=0.0;Grad.Y=0.0;Grad.Z=0.0;
        r1=0.0;r2=0.0;
        npc=Grids->Cells[i].NodePerCell;
        if (npc==4)
        {
            fi_n0=(*fi_nodal)[Grids->Cells[i].NodeList[0]];
            fi_n1=(*fi_nodal)[Grids->Cells[i].NodeList[1]];
            fi_n2=(*fi_nodal)[Grids->Cells[i].NodeList[2]];
            fi_n3=(*fi_nodal)[Grids->Cells[i].NodeList[3]];
            N0=Grids->Nodes[Grids->Cells[i].NodeList[0]].Pos;
            N1=Grids->Nodes[Grids->Cells[i].NodeList[1]].Pos;
            N2=Grids->Nodes[Grids->Cells[i].NodeList[2]].Pos;
            N3=Grids->Nodes[Grids->Cells[i].NodeList[3]].Pos;
            PC=Grids->Cells[i].CellCenter;
            t1=Normalized(Minus(Grids->Cells[i].EdgeCenter[0],Grids->Cells[i].EdgeCenter[2]));
            n=Normalized(Grids->Cells[i].Area);
            t2=Cross(n,t1);

            a=Dot(t1,Minus(N0,PC));b=Dot(t2,Minus(N0,PC));c=(fi_n0-(*fi)[i]);///Mag(Minus(N0,PC));
            d=Dot(t1,Minus(N1,PC));e=Dot(t2,Minus(N1,PC));f=(fi_n1-(*fi)[i]);///Mag(Minus(N1,PC));
            r1+=((e*c-b*f))/((a*e-b*d));r2+=((a*f-c*d))/((a*e-b*d));

//            a=Dot(t1,Minus(N0,PC));b=Dot(t2,Minus(N0,PC));c=(fi_n0-(*fi)[i]);///Mag(Minus(N0,PC));
//            d=Dot(t1,Minus(N2,PC));e=Dot(t2,Minus(N2,PC));f=(fi_n2-(*fi)[i]);///Mag(Minus(N2,PC));
//            r1+=((e*c-b*f))/((a*e-b*d));r2+=((a*f-c*d))/((a*e-b*d));

            a=Dot(t1,Minus(N3,PC));b=Dot(t2,Minus(N3,PC));c=(fi_n3-(*fi)[i]);///Mag(Minus(N3,PC));
            d=Dot(t1,Minus(N0,PC));e=Dot(t2,Minus(N0,PC));f=(fi_n0-(*fi)[i]);///Mag(Minus(N0,PC));
            r1+=((e*c-b*f))/((a*e-b*d));r2+=((a*f-c*d))/((a*e-b*d));

            a=Dot(t1,Minus(N1,PC));b=Dot(t2,Minus(N1,PC));c=(fi_n1-(*fi)[i]);///Mag(Minus(N1,PC));
            d=Dot(t1,Minus(N2,PC));e=Dot(t2,Minus(N2,PC));f=(fi_n2-(*fi)[i]);///Mag(Minus(N2,PC));
            r1+=((e*c-b*f))/((a*e-b*d));r2+=((a*f-c*d))/((a*e-b*d));

//            a=Dot(t1,Minus(N1,PC));b=Dot(t2,Minus(N1,PC));c=(fi_n1-(*fi)[i]);///Mag(Minus(N1,PC));
//            d=Dot(t1,Minus(N3,PC));e=Dot(t2,Minus(N3,PC));f=(fi_n3-(*fi)[i]);///Mag(Minus(N3,PC));
//            r1+=((e*c-b*f))/((a*e-b*d));r2+=((a*f-c*d))/((a*e-b*d));

            a=Dot(t1,Minus(N2,PC));b=Dot(t2,Minus(N2,PC));c=(fi_n2-(*fi)[i]);///Mag(Minus(N2,PC));
            d=Dot(t1,Minus(N3,PC));e=Dot(t2,Minus(N3,PC));f=(fi_n3-(*fi)[i]);///Mag(Minus(N3,PC));
            r1+=((e*c-b*f))/((a*e-b*d));r2+=((a*f-c*d))/((a*e-b*d));

            r1=r1/4.0;r2=r2/4.0;
            Grad=Sum(SVP(r1,t1),SVP(r2,t2));
        }
        else
        {
            fi_n0=(*fi_nodal)[Grids->Cells[i].NodeList[0]];
            fi_n1=(*fi_nodal)[Grids->Cells[i].NodeList[1]];
            fi_n2=(*fi_nodal)[Grids->Cells[i].NodeList[2]];
            N0=Grids->Nodes[Grids->Cells[i].NodeList[0]].Pos;
            N1=Grids->Nodes[Grids->Cells[i].NodeList[1]].Pos;
            N2=Grids->Nodes[Grids->Cells[i].NodeList[2]].Pos;
            PC=Grids->Cells[i].CellCenter;
            t1=Normalized(Minus(Grids->Cells[i].EdgeCenter[0],PC));
            n=Normalized(Grids->Cells[i].Area);
            t2=Cross(n,t1);

            a=Dot(t1,Minus(N0,PC));b=Dot(t2,Minus(N0,PC));c=(fi_n0-(*fi)[i]);///Mag(Minus(N0,PC));
            d=Dot(t1,Minus(N1,PC));e=Dot(t2,Minus(N1,PC));f=(fi_n1-(*fi)[i]);///Mag(Minus(N1,PC));
            r1+=(e*c-b*f)/(a*e-b*d);r2+=(a*f-c*d)/(a*e-b*d);

            a=Dot(t1,Minus(N1,PC));b=Dot(t2,Minus(N1,PC));c=(fi_n1-(*fi)[i]);///Mag(Minus(N1,PC));
            d=Dot(t1,Minus(N2,PC));e=Dot(t2,Minus(N2,PC));f=(fi_n2-(*fi)[i]);///Mag(Minus(N2,PC));
            r1+=(e*c-b*f)/(a*e-b*d);r2+=(a*f-c*d)/(a*e-b*d);

            a=Dot(t1,Minus(N2,PC));b=Dot(t2,Minus(N2,PC));c=(fi_n2-(*fi)[i]);///Mag(Minus(N2,PC));
            d=Dot(t1,Minus(N0,PC));e=Dot(t2,Minus(N0,PC));f=(fi_n0-(*fi)[i]);///Mag(Minus(N0,PC));
            r1+=(e*c-b*f)/(a*e-b*d);r2+=(a*f-c*d)/(a*e-b*d);

            r1=r1/3.0;r2=r2/3.0;
            Grad=Sum(SVP(r1,t1),SVP(r2,t2));
        }

//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        vBT=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        vBT=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        vREF=Minus(inf_Vel,vBT);   /// not correct, should be summed!!!

        Grad=Minus(vREF,Grad);
        Grad=Minus(Grad,SVP(Dot(Grad,n),n));

        (*u)[i]=Grad.X;
        (*v)[i]=Grad.Y;
        (*w)[i]=Grad.Z;
    }
}

//void GetMassCenter (struct MyVector *MassCen,struct Grid *Grids)
//{
//    int i;
//    double AMag,TotalAMag=0.0,Xc=0.0,Yc=0.0,Zc=0.0;
//    struct MyVector CellCen;
//
//    for (i=0;i<Grids->NC;i++)
//    {
//        AMag=Mag(Grids->Cells[i].Area);
//        TotalAMag+=AMag;
//        CellCen=Grids->Cells[i].CellCenter;
//        Xc+=(CellCen.X)*AMag;
//        Yc+=(CellCen.Y)*AMag;
//        Zc+=(CellCen.Z)*AMag;
//    }
////    (*MassCen).X=Xc/TotalAMag;
////    (*MassCen).Y=Yc/TotalAMag;
////    (*MassCen).Z=Zc/TotalAMag;
//    (*MassCen).X=Grids->DynaVariables.MassCenter.X;
//    (*MassCen).Y=Grids->DynaVariables.MassCenter.Y;
//    (*MassCen).Z=Grids->DynaVariables.MassCenter.Z;
//}

//void CopyToLaspackMatrix(QMatrix *lmc,struct MatrixCoefficient *g,struct Grid *Grids)
//{
//    int i,j;
////    int n=(Grids->NC)+(Grids->NTE)*(Grids->NST);
//    int n=(Grids->NC);
//
//    for (i=0;i<n;i++)
//    {
////        Q_SetLen(lmc,i+1,n);
//        for(j=0;j<n;j++)
//            Q_SetEntry(lmc,i+1,j,j+1,g->Elem[i][j]);
//    }
//}
//
//void CopyToLaspackVector(Vector *GlobU,double *U,struct Grid *Grids)
//{
//    int PCell;
////    int n=(Grids->NC)+(Grids->NTE)*(Grids->NST);
//    int n=(Grids->NC);
//    for (PCell=0;PCell<n;PCell++)
//        V_SetCmp(GlobU,PCell+1,U[PCell]);
//}
//
//void CopyFromLaspackVector(double *U,Vector *GlobU,struct Grid *Grids)
//{
//    int PCell;
////    int n=(Grids->NC)+(Grids->NTE)*(Grids->NST);
//    int n=(Grids->NC);
//    for (PCell=0;PCell<n;PCell++)
//        U[PCell]=V_GetCmp(GlobU,PCell+1);
//}

