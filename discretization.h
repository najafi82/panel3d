//void InitializeStaticArrays
//               (double ***g,double ***phi,double ***rhs,
//                double ***p,double ***bu,double ***bv,double ***bw, ///press, Body Purt. Velocity
//                int nc, int nte, int nst)
//{
//    int i;
//
//    *g=AllocMatrix(1,nc+(nte*nst),1,nc+(nte*nst));
//    *phi=AllocMatrix(1,nc+(nte*nst),1,1);
//    *rhs=AllocMatrix(1,nc+(nte*nst),1,1);
//
//    *p=AllocMatrix(1,nc,1,1);
//    *bu=AllocMatrix(1,nc,1,1);
//    *bv=AllocMatrix(1,nc,1,1);
//    *bw=AllocMatrix(1,nc,1,1);
//
//    for(i=1;i<=(nc+(nte*nst));i++)
//        (*phi)[i][1]=0.0;
//}

void AllocateStaticArrays
               (struct MatrixCoefficient *a,struct MatrixCoefficient *b,
                struct MatrixCoefficient *g,double **phi,double **phi_old,double **phibar,double **phibar_old,double **phi_nodal,double **rhs,double **src,
                double **p,double **bu,double **bv,double **bw, ///pressure and Body Purt. Velocity
                int nc, int nn, int nte)
{
    int i;

//    (*g).Elem =(double**) calloc((nc+(nte*nst)),sizeof(double*));
    (*a).Elem =(double**) calloc(nc,sizeof(double*));
    (*b).Elem =(double**) calloc(nc,sizeof(double*));
    for(i=0;i<nc;i++)
    {
        (*a).Elem[i]=(double*) calloc(nc,sizeof(double));
        (*b).Elem[i]=(double*) calloc(nc,sizeof(double));
    }

//    for(i=0;i<(nc+(nte*nst));i++)
//        (*g).Elem[i]=(double*) calloc((nc+(nte*nst)),sizeof(double));

    *p=(double*) calloc(nc,sizeof(double));
    *bu=(double*) calloc(nc,sizeof(double));
    *bv=(double*) calloc(nc,sizeof(double));
    *bw=(double*) calloc(nc,sizeof(double));

    /// phi,rhs va Laspack ghablan bejaye nc, nc+(nte*nst) alloc shode boodand!!!
    *phi=(double*) calloc(nc,sizeof(double));
    *phi_old=(double*) calloc(nc,sizeof(double));
    *phibar=(double*) calloc(nc,sizeof(double));
    *phibar_old=(double*) calloc(nc,sizeof(double));
    *phi_nodal=(double*) calloc(nn,sizeof(double));
    *rhs=(double*) calloc(nc,sizeof(double));
    *src=(double*) calloc(nc,sizeof(double));

//    Q_Constr(lmc,"LaspackMC",nc,False,Rowws,Normal,True);
//    V_Constr(st,"ST",nc,Normal,True);
//    V_Constr(fee,"Fee",nc,Normal,True);

//    for (i=0;i<nc;i++)
//        Q_SetLen(lmc,i+1,nc);

}

void AllocateDynamicArrays
               (double ***m,double ***iu,double ***iv,double ***iw,///Miu and Induced Velocity of Vortex Nodes
                int nc, int nte, int nst)
{
    int i;

    *m =(double**) calloc( nst   *MAX_SIM_STEP, sizeof(double*));
    *iu=(double**) calloc((nst+1)*MAX_SIM_STEP, sizeof(double*));
    *iv=(double**) calloc((nst+1)*MAX_SIM_STEP, sizeof(double*));
    *iw=(double**) calloc((nst+1)*MAX_SIM_STEP, sizeof(double*));

    for (i=0;i<nst*MAX_SIM_STEP;i++)
         (*m)[i]=(double*) calloc(nte, sizeof(double));

    for (i=0;i<(nst+1)*MAX_SIM_STEP;i++)
    {
        (*iu)[i]=(double*) calloc(nte, sizeof(double));
        (*iv)[i]=(double*) calloc(nte, sizeof(double));
        (*iw)[i]=(double*) calloc(nte, sizeof(double));
    }
}


//void Free(double ***g,double ***phi,double ***u,double ***v,double ***w,double ***p,
//          double ***rhs,int NoRow)
//{
//    FreeMatrix(*g,1,NoRow,1,NoRow);
//    FreeMatrix(*rhs,1,NoRow,1,1);
//    FreeMatrix(*phi,1,NoRow,1,1);
////    FreeMatrix(*U,1,NoRow,1,1);
////    FreeMatrix(*V,1,NoRow,1,1);
////    FreeMatrix(*W,1,NoRow,1,1);
////    FreeMatrix(*P,1,NoRow,1,1);
//}

/// //////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////
void PanelMapping(struct MyVector *q0,struct MyVector *q1,
                  struct MyVector *q2,struct MyVector *q3,
                  struct MyVector *SPanel)
{
    *q0=SVP(0.25,Sum(SPanel[0],Sum(SPanel[3],Sum(SPanel[1],SPanel[2]))));
    *q1=SVP(0.25,Minus(Sum(SPanel[0],SPanel[3]),Sum(SPanel[1],SPanel[2])));
    *q2=SVP(0.25,Minus(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
    *q3=SVP(0.25,Minus(Sum(SPanel[0],SPanel[2]),Sum(SPanel[3],SPanel[1])));
}

struct MyVector Q(double xi, double et, struct MyVector *SPanel)
{
    struct MyVector q0,q1,q2,q3;

    PanelMapping(&q0,&q1,&q2,&q3,SPanel);
    struct MyVector V0=q0;
    struct MyVector V1=SVP(xi,q1);
    struct MyVector V2=SVP(et,q2);
    struct MyVector V3=SVP(et*xi,q3);

    return Sum(V0,Sum(V1,Sum(V2,V3)));
};

void calcSourceDipole_Quad_Std(double *Src, double *Dpl, struct MyVector FPoint, struct MyVector *SPanel)
{
    struct MyVector q1,q2,q3;

    struct MyVector Q_00,Q_11,Q_m11,Q_1m1,Q_m1m1;
    struct MyVector R_00,R_11,R_m11,R_1m1,R_m1m1;
    struct MyVector S1_00,S1_11,S1_m11,S1_1m1,S1_m1m1;
    struct MyVector S2_00,S2_11,S2_m11,S2_1m1,S2_m1m1;
    struct MyVector S3;

    double ID_11,ID_m11,ID_1m1,ID_m1m1;
    double IS_11,IS_m11,IS_1m1,IS_m1m1;
    double ep_11,ep_m11,ep_1m1,ep_m1m1;

    double Phi_D;
    double IDSUM,EPSMPLY;

    q1=SVP(0.25,Minus(Sum(SPanel[0],SPanel[3]),Sum(SPanel[1],SPanel[2])));
    q2=SVP(0.25,Minus(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
    q3=SVP(0.25,Minus(Sum(SPanel[0],SPanel[2]),Sum(SPanel[3],SPanel[1])));

    Q_00=Q(0.,0.,SPanel);
    Q_11=Q(1.,1.,SPanel);
    Q_m11=Q(-1.,1.,SPanel);
    Q_1m1=Q(1.,-1,SPanel);
    Q_m1m1=Q(-1.,-1.,SPanel);

    R_00=Minus(Q_00,FPoint);
    R_11=Minus(Q_11,FPoint);
    R_m11=Minus(Q_m11,FPoint);
    R_1m1=Minus(Q_1m1,FPoint);
    R_m1m1=Minus(Q_m1m1,FPoint);

    S1_00=q1;
    S1_11=Sum(q1,q3);
    S1_m11=Sum(q1,q3);
    S1_1m1=Minus(q1,q3);
    S1_m1m1=Minus(q1,q3);

    S2_00=q2;
    S2_11=Sum(q2,q3);
    S2_m11=Minus(q2,q3);
    S2_1m1=Sum(q2,q3);
    S2_m1m1=Minus(q2,q3);

    S3=Normalized(Cross(S1_00,S2_00));

    ID_11=(0.25/pi)*atan((Dot(Cross(R_11,S1_11),Cross(R_11,S2_11)))/((Mag(R_11)*fabs(Dot(R_11,Cross(S1_11,S2_11))))));
    ID_m11=(0.25/pi)*atan((Dot(Cross(R_m11,S1_m11),Cross(R_m11,S2_m11)))/((Mag(R_m11)*fabs(Dot(R_m11,Cross(S1_m11,S2_m11))))));
    ID_1m1=(0.25/pi)*atan((Dot(Cross(R_1m1,S1_1m1),Cross(R_1m1,S2_1m1)))/((Mag(R_1m1)*fabs(Dot(R_1m1,Cross(S1_1m1,S2_1m1))))));
    ID_m1m1=(0.25/pi)*atan((Dot(Cross(R_m1m1,S1_m1m1),Cross(R_m1m1,S2_m1m1)))/((Mag(R_m1m1)*fabs(Dot(R_m1m1,Cross(S1_m1m1,S2_m1m1))))));

    IS_11=(-0.25/pi)* ((Dot(Cross(R_11,S2_11),S3)/Mag(S2_11))*asinh((Dot(R_11,S2_11))/(Mag(Cross(R_11,S2_11))))
                     -(Dot(Cross(R_11,S1_11),S3)/Mag(S1_11))*asinh((Dot(R_11,S1_11))/(Mag(Cross(R_11,S1_11)))));
    IS_m11=(-0.25/pi)* ((Dot(Cross(R_m11,S2_m11),S3)/Mag(S2_m11))*asinh((Dot(R_m11,S2_m11))/(Mag(Cross(R_m11,S2_m11))))
                     -(Dot(Cross(R_m11,S1_m11),S3)/Mag(S1_m11))*asinh((Dot(R_m11,S1_m11))/(Mag(Cross(R_m11,S1_m11)))));
    IS_1m1=(-0.25/pi)* ((Dot(Cross(R_1m1,S2_1m1),S3)/Mag(S2_1m1))*asinh((Dot(R_1m1,S2_1m1))/(Mag(Cross(R_1m1,S2_1m1))))
                     -(Dot(Cross(R_1m1,S1_1m1),S3)/Mag(S1_1m1))*asinh((Dot(R_1m1,S1_1m1))/(Mag(Cross(R_1m1,S1_1m1)))));
    IS_m1m1=(-0.25/pi)* ((Dot(Cross(R_m1m1,S2_m1m1),S3)/Mag(S2_m1m1))*asinh((Dot(R_m1m1,S2_m1m1))/(Mag(Cross(R_m1m1,S2_m1m1))))
                     -(Dot(Cross(R_m1m1,S1_m1m1),S3)/Mag(S1_m1m1))*asinh((Dot(R_m1m1,S1_m1m1))/(Mag(Cross(R_m1m1,S1_m1m1)))));

    ep_11=(Dot(R_11,Cross(S1_11,S2_11)))/(fabs(Dot(R_11,Cross(S1_11,S2_11))));
    ep_m11=(Dot(R_m11,Cross(S1_m11,S2_m11)))/(fabs(Dot(R_m11,Cross(S1_m11,S2_m11))));
    ep_1m1=(Dot(R_1m1,Cross(S1_1m1,S2_1m1)))/(fabs(Dot(R_1m1,Cross(S1_1m1,S2_1m1))));
    ep_m1m1=(Dot(R_m1m1,Cross(S1_m1m1,S2_m1m1)))/(fabs(Dot(R_m1m1,Cross(S1_m1m1,S2_m1m1))));

    Phi_D = ep_11*ID_11 - ep_1m1*ID_1m1 - ep_m11*ID_m11 + ep_m1m1*ID_m1m1 ;

    IDSUM  =ID_11 + ID_1m1 + ID_m11 + ID_m1m1;
    EPSMPLY=ep_11 * ep_1m1 * ep_m11 * ep_m1m1;

    if ((IDSUM<1.0+1.0e-7) && (IDSUM>1.0-1.0e-7))
        (*Dpl)=0.0;
    else if (EPSMPLY<0.0)
        (*Dpl)=(Phi_D - 0.5*sgn(Phi_D));
    else
        (*Dpl)=Phi_D;

    (*Src)= IS_11 - IS_1m1 - IS_m11 + IS_m1m1 - (*Dpl)*Dot(R_00,S3);
}

void calcSourceDipole_Quad_APAME(double *Src, double *Dpl, struct MyVector FPoint, struct MyVector *SPanel)
{
    double a1,a2,a3,a4,b1,b2,b3,b4,e1,e2,e3,e4,r1,r2,r3,r4;
    double d1,d2,d3,d4,h1,h2,h3,h4,F,G;
    double cpx1,cpx2,cpx3,cpx4,cpy1,cpy2,cpy3,cpy4;
    double x21,x32,x43,x14,y21,y32,y43,y14;
    (*Src)=0.0;
    (*Dpl)=0.0;

    double MMM[3][3];
    struct MyVector N0,N1,N2,N3;
    struct MyVector Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
    struct MyVector Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));
    struct MyVector u_Vec=Normalized(Minus(SPanel[1],SPanel[0]));
    struct MyVector n_Vec=Normalized(Area);
    struct MyVector o_Vec=Cross(n_Vec,u_Vec);

    MMM[0][0]=u_Vec.X;
    MMM[0][1]=u_Vec.Y;
    MMM[0][2]=u_Vec.Z;
    MMM[1][0]=o_Vec.X;
    MMM[1][1]=o_Vec.Y;
    MMM[1][2]=o_Vec.Z;
    MMM[2][0]=n_Vec.X;
    MMM[2][1]=n_Vec.Y;
    MMM[2][2]=n_Vec.Z;

    N0=MatrixDotVector(MMM,Minus(SPanel[0],Center));
    N1=MatrixDotVector(MMM,Minus(SPanel[1],Center));
    N2=MatrixDotVector(MMM,Minus(SPanel[2],Center));
    N3=MatrixDotVector(MMM,Minus(SPanel[3],Center));

    d1=Mag(Minus(N1,N0));
    d2=Mag(Minus(N2,N1));
    d3=Mag(Minus(N3,N2));
    d4=Mag(Minus(N0,N3));
    double S=Mag(Area);
    double dist_x = FPoint.X-Center.X;
    double dist_y = FPoint.Y-Center.Y;
    double dist_z = FPoint.Z-Center.Z;
    double cpx = dist_x*MMM[0][0] + dist_y*MMM[0][1] + dist_z*MMM[0][2];
    double cpy = dist_x*MMM[1][0] + dist_y*MMM[1][1] + dist_z*MMM[1][2];
    double cpz = dist_x*MMM[2][0] + dist_y*MMM[2][1] + dist_z*MMM[2][2];
    double rad = pow(cpx,2.0)+pow(cpy,2.0)+pow(cpz,2.0);
    double fourpi=1.0/(4.0*pi);

//    if (rad > FF_sqr)
//    {
//        (*Dpl) = fourpi*S*cpz*pow(rad,(-1.5));
//        (*Src) = fourpi*S/sqrt(rad);
//    }
//    else
//    {
        cpx1 = cpx - N0.X;
        cpx2 = cpx - N1.X;
        cpx3 = cpx - N2.X;
        cpx4 = cpx - N3.X;
        cpy1 = cpy - N0.Y;
        cpy2 = cpy - N1.Y;
        cpy3 = cpy - N2.Y;
        cpy4 = cpy - N3.Y;
        e1 = pow(cpx1,2.0)+pow(cpz,2.0);
        e2 = pow(cpx2,2.0)+pow(cpz,2.0);
        e3 = pow(cpx3,2.0)+pow(cpz,2.0);
        e4 = pow(cpx4,2.0)+pow(cpz,2.0);
        r1 = sqrt(e1 + pow(cpy1,2.0));
        r2 = sqrt(e2 + pow(cpy2,2.0));
        r3 = sqrt(e3 + pow(cpy3,2.0));
        r4 = sqrt(e4 + pow(cpy4,2.0));
        x21 = N1.X-N0.X;
        x32 = N2.X-N1.X;
        x43 = N3.X-N2.X;
        x14 = N0.X-N3.X;
        y21 = N1.Y-N0.Y;
        y32 = N2.Y-N1.Y;
        y43 = N3.Y-N2.Y;
        y14 = N0.Y-N3.Y;

        if (fabs(cpz) < 1.e-7)
        {
            if (d1 < 1.e-7)
                b1 = 0.0;
            else
                b1 = (cpx1*y21-cpy1*x21)/d1*log((r1+r2+d1)/(r1+r2-d1));


            if (d2 < 1.e-7)
                b2 = 0.0;
            else
                b2 = (cpx2*y32-cpy2*x32)/d2*log((r2+r3+d2)/(r2+r3-d2));


            if (d3 < 1.e-7)
                b3 = 0.0;
            else
                b3 = (cpx3*y43-cpy3*x43)/d3*log((r3+r4+d3)/(r3+r4-d3));


            if (d4 < 1.e-7)
                b4 = 0.0;
            else
                b4 = (cpx4*y14-cpy4*x14)/d4*log((r4+r1+d4)/(r4+r1-d4));

            (*Dpl) = 0.0;
            (*Src) = -(b1+b2+b3+b4)*fourpi;
        }
        else
        {
            h1 = cpx1*cpy1;
            h2 = cpx2*cpy2;
            h3 = cpx3*cpy3;
            h4 = cpx4*cpy4;

            if (d1 < 1.e-7)
            {
                a1 = 0.0;
                b1 = 0.0;
            }
            else
            {
                F = y21*e1 - x21*h1;
                G = y21*e2 - x21*h2;
                a1 = atan2(cpz*x21*(F*r2-G*r1), pow(cpz,2.0)*pow(x21,2.0)*r1*r2+F*G);
                b1 = (cpx1*y21-cpy1*x21)/d1*log((r1+r2+d1)/(r1+r2-d1));
            }

            if (d2 < 1.e-7)
            {
                a2 = 0.0;
                b2 = 0.0;
            }
            else
            {
                F = y32*e2 - x32*h2;
                G = y32*e3 - x32*h3;
                a2 = atan2(cpz*x32*(F*r3-G*r2), pow(cpz,2.0)*pow(x32,2.0)*r2*r3+F*G);
                b2 = (cpx2*y32-cpy2*x32)/d2*log((r2+r3+d2)/(r2+r3-d2));
            }
            if (d3 < 1.e-7)
            {
                a3 = 0.0;
                b3 = 0.0;
            }
            else
            {
                F = y43*e3 - x43*h3;
                G = y43*e4 - x43*h4;
                a3 = atan2(cpz*x43*(F*r4-G*r3), pow(cpz,2.0)*pow(x43,2.0)*r3*r4+F*G);
                b3 = (cpx3*y43-cpy3*x43)/d3*log((r3+r4+d3)/(r3+r4-d3));
            }
            if (d4 < 1.e-7)
            {
                a4 = 0.0;
                b4 = 0.0;
            }
            else
            {
                F = y14*e4 - x14*h4;
                G = y14*e1 - x14*h1;
                a4 = atan2(cpz*x14*(F*r1-G*r4), pow(cpz,2.0)*pow(x14,2.0)*r4*r1+F*G);
                b4 = (cpx4*y14-cpy4*x14)/d4*log((r4+r1+d4)/(r4+r1-d4));
            }
            (*Dpl) = -(a1+a2+a3+a4)*fourpi;
            (*Src) = -(b1+b2+b3+b4)*fourpi-cpz*(*Dpl);
        }
//    }
}

void calcSourceDipole_Quad_Vortexje(double *Src, double *Dpl, struct MyVector FPoint, struct MyVector *SPanel)
{
    int i;
    int prev_idx;
    double d,z,m,e1,e2,r1,r2,h1,h2,u,v,delta_theta;
    struct MyVector node_a;
    struct MyVector node_b;
    (*Src)=0.0;
    (*Dpl)=0.0;

    double MMM[3][3];
    struct MyVector F;
    struct MyVector N[4];

    struct MyVector Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
    struct MyVector Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));

    struct MyVector u_Vec=Normalized(Minus(SPanel[1],SPanel[0]));
    struct MyVector n_Vec=Normalized(Area);
    struct MyVector o_Vec=Cross(n_Vec,u_Vec);

    MMM[0][0]=u_Vec.X;
    MMM[0][1]=u_Vec.Y;
    MMM[0][2]=u_Vec.Z;
    MMM[1][0]=o_Vec.X;
    MMM[1][1]=o_Vec.Y;
    MMM[1][2]=o_Vec.Z;
    MMM[2][0]=n_Vec.X;
    MMM[2][1]=n_Vec.Y;
    MMM[2][2]=n_Vec.Z;

    N[0]=MatrixDotVector(MMM,Minus(SPanel[0],Center));
    N[1]=MatrixDotVector(MMM,Minus(SPanel[1],Center));
    N[2]=MatrixDotVector(MMM,Minus(SPanel[2],Center));
    N[3]=MatrixDotVector(MMM,Minus(SPanel[3],Center));

    F=MatrixDotVector(MMM,Minus(FPoint,Center));
    z = F.Z;

    if (fabs(z)<1.e-7)
    {
        for ( i = 0; i < 4; i++)
        {
            if (i == 0)
                prev_idx = 3;
            else
                prev_idx = i - 1;

            node_a = N[prev_idx];
            node_b = N[i];

            d = sqrt(pow(node_b.X - node_a.X, 2.0) + pow(node_b.Y - node_a.Y, 2.0));

            if (d < 1.0e-7)
            {
                (*Dpl)+=0.0;
                (*Src)+=0.0;
            }
            else
            {
                e1 = pow(F.X - node_a.X, 2.0) + pow(z, 2.0);
                e2 = pow(F.X - node_b.X, 2.0) + pow(z, 2.0);
                r1 = sqrt(e1 + pow(F.Y - node_a.Y, 2.0));
                r2 = sqrt(e2 + pow(F.Y - node_b.Y, 2.0));
                (*Dpl)+=0.0;
                (*Src)+=((F.X - node_a.X)*(node_b.Y-node_a.Y)-(F.Y-node_a.Y)*(node_b.X-node_a.X))/d * log((r1 + r2 + d) / (r1 + r2 - d));
            }
        }
    }
    else
    {
        for ( i = 0; i < 4; i++)
        {
            if (i == 0)
                prev_idx = 3;
            else
                prev_idx = i - 1;

            node_a = N[prev_idx];
            node_b = N[i];

            d = sqrt(pow(node_b.X - node_a.X, 2.0) + pow(node_b.Y - node_a.Y, 2.0));

            if (d < 1.0e-7)
            {
                (*Src)+=0.0;
                (*Dpl)+=0.0;
            }
            else
            {
                m = (node_b.Y - node_a.Y) / (node_b.X - node_a.X);
                e1 = pow(F.X - node_a.X, 2.0) + pow(z, 2.0);
                e2 = pow(F.X - node_b.X, 2.0) + pow(z, 2.0);
                r1 = sqrt(e1 + pow(F.Y - node_a.Y, 2.0));
                r2 = sqrt(e2 + pow(F.Y - node_b.Y, 2.0));
                h1 = (F.X - node_a.X) * (F.Y - node_a.Y);
                h2 = (F.X - node_b.X) * (F.Y - node_b.Y);
                u = (m * e1 - h1) / (z * r1);
                v = (m * e2 - h2) / (z * r2);

                if (u==v)
                    delta_theta = 0.0;
                else
                    delta_theta = atan2(u - v, 1.0 + u * v);

                (*Src)+=((F.X - node_a.X)*(node_b.Y-node_a.Y)-(F.Y-node_a.Y)*(node_b.X-node_a.X))/d * log((r1 + r2 + d) / (r1 + r2 - d)) + fabs(z) * delta_theta;
                (*Dpl)+=-delta_theta;
            }
        }
    }

    (*Src) *= 1.0/(-4.0*pi);
    (*Dpl) *= 1.0/(4.0*pi);
}

void calcSourceDipole_Quad_Vortexje_Original(double *Src, double *Dpl, struct MyVector FPoint, struct MyVector *SPanel)
{
    int i;
    int prev_idx;
    double d,z,m,e1,e2,r1,r2,h1,h2,u,v,delta_theta;
    struct MyVector node_a;
    struct MyVector node_b;
    (*Src)=0.0;
    (*Dpl)=0.0;

    double MMM[3][3];
    struct MyVector F;
    struct MyVector N[4];

    struct MyVector Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
    struct MyVector Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));

    struct MyVector u_Vec=Normalized(Minus(SPanel[1],SPanel[0]));
    struct MyVector n_Vec=Normalized(Area);
    struct MyVector o_Vec=Cross(n_Vec,u_Vec);

    MMM[0][0]=u_Vec.X;
    MMM[0][1]=u_Vec.Y;
    MMM[0][2]=u_Vec.Z;
    MMM[1][0]=o_Vec.X;
    MMM[1][1]=o_Vec.Y;
    MMM[1][2]=o_Vec.Z;
    MMM[2][0]=n_Vec.X;
    MMM[2][1]=n_Vec.Y;
    MMM[2][2]=n_Vec.Z;

    N[0]=MatrixDotVector(MMM,Minus(SPanel[0],Center));
    N[1]=MatrixDotVector(MMM,Minus(SPanel[1],Center));
    N[2]=MatrixDotVector(MMM,Minus(SPanel[2],Center));
    N[3]=MatrixDotVector(MMM,Minus(SPanel[3],Center));

    F=MatrixDotVector(MMM,Minus(FPoint,Center));
    z = F.Z;

    for ( i = 0; i < 4; i++)
    {
        if (i == 0)
            prev_idx = 3;
        else
            prev_idx = i - 1;

        node_a = N[prev_idx];
        node_b = N[i];

        d = sqrt(pow(node_b.X - node_a.X, 2.0) + pow(node_b.Y - node_a.Y, 2.0));

        if (d < 1.0e-7)
        {
            (*Src)+=0.0;
            (*Dpl)+=0.0;
        }
        else
        {
            m = (node_b.Y - node_a.Y) / (node_b.X - node_a.X);
            e1 = pow(F.X - node_a.X, 2.0) + pow(z, 2.0);
            e2 = pow(F.X - node_b.X, 2.0) + pow(z, 2.0);
            r1 = sqrt(e1 + pow(F.Y - node_a.Y, 2.0));
            r2 = sqrt(e2 + pow(F.Y - node_b.Y, 2.0));
            h1 = (F.X - node_a.X) * (F.Y - node_a.Y);
            h2 = (F.X - node_b.X) * (F.Y - node_b.Y);
            u = (m * e1 - h1) / (z * r1);
            v = (m * e2 - h2) / (z * r2);

            if (fabs(u-v)<1.0e-7) ///if(u==v)
                delta_theta = 0.0;
            else
                delta_theta = atan2(u - v, 1.0 + u * v);

            (*Src)+=((F.X - node_a.X)*(node_b.Y-node_a.Y)-(F.Y-node_a.Y)*(node_b.X-node_a.X))/d * log((r1 + r2 + d) / (r1 + r2 - d)) - fabs(z) * delta_theta;
            (*Dpl)+=delta_theta;
        }
    }

    (*Src) *= -1.0/(4.0*pi);
    (*Dpl) *= 1.0/(4.0*pi);
}

void calcSourceDipole_Delhommau(double *Src, double *Dpl, struct MyVector FPoint, struct MyVector *SPanel, int m)
{
    int i,ip1;
    struct MyVector Area,NormalArea,Center,DummyVector;
    double Z;
    double d[m];
    double Y[m];
    double R[m];
    double Nl[m];
    double Dl[m];
    double Nt[m];
    double Dt[m];
    double J[m];
    double K[m];

    (*Src)=(*Dpl)=0.0;

    if (m==4)
    {
        Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
        Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));
    }
    else
    {
        Area=SVP(0.5,Cross(Minus(SPanel[1],SPanel[0]),Minus(SPanel[2],SPanel[0])));
        Center=SVP(1.0/3.0,Sum3(SPanel[0],SPanel[1],SPanel[2]));
    }

    NormalArea=Normalized(Area);
    Z=Dot(NormalArea,Minus(FPoint,Center));
//    Z=Dot(NormalArea,Minus(FPoint,Sum(Center,SVP(-1.0e-6,NormalArea))));

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;
        DummyVector=Cross(NormalArea,Minus(SPanel[ip1],SPanel[i]));
        d[i]=Mag(Minus(SPanel[ip1],SPanel[i]));
        Y[i]=Dot(Minus(FPoint,SPanel[i]),DummyVector)/d[i];
        R[i]=Mag(Minus(FPoint,SPanel[i]));
    }

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;
        Nl[i]=R[ip1]+R[i]+d[i];
        Dl[i]=R[ip1]+R[i]-d[i];

        Nt[i]=2.0*Y[i]*d[i];
        Dt[i]=pow(R[ip1]+R[i],2.0)-pow(d[i],2.0)+2.0*fabs(Z)*(R[ip1]+R[i]);
    }

    for (i=0;i<m;i++)
    {
        J[i]=2.0*fabs(Z)*atan(Nt[i]/Dt[i])-Y[i]*log((Nl[i])/(Dl[i]));
        K[i]=2.0*sgn(Z)*atan(Nt[i]/Dt[i]);
    }

    for (i=0;i<m;i++)
    {
        (*Src)+=J[i]/(-4.0*pi);
        (*Dpl)+=K[i]/(4.0*pi);
    }
}

void calcSource_Delhommau(double *Src, struct MyVector FPoint, struct MyVector *SPanel, int m)
{
    int i,ip1;
    struct MyVector Area,NormalArea,Center,DummyVector;
    double Z;
    double d[m];
    double Y[m];
    double R[m];
    double Nl[m];
    double Dl[m];
    double Nt[m];
    double Dt[m];
    double J[m];
    double K[m];

    (*Src)=0.0;

    if (m==4)
    {
        Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
        Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));
    }
    else
    {
        Area=SVP(0.5,Cross(Minus(SPanel[1],SPanel[0]),Minus(SPanel[2],SPanel[0])));
        Center=SVP(1.0/3.0,Sum3(SPanel[0],SPanel[1],SPanel[2]));
    }

    NormalArea=SVP(1.0/Mag(Area),Area);
    Z=Dot(NormalArea,Minus(FPoint,Center));
//    Z=Dot(NormalArea,Minus(FPoint,Sum(Center,SVP(-1.0e-6,NormalArea))));

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;
        DummyVector=Cross(NormalArea,Minus(SPanel[ip1],SPanel[i]));
        d[i]=Mag(Minus(SPanel[ip1],SPanel[i]));
        Y[i]=Dot(Minus(FPoint,SPanel[i]),DummyVector)/d[i];
        R[i]=Mag(Minus(FPoint,SPanel[i]));
    }

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;
        Nl[i]=R[ip1]+R[i]+d[i];
        Dl[i]=R[ip1]+R[i]-d[i];

        Nt[i]=2.0*Y[i]*d[i];
        Dt[i]=pow(R[ip1]+R[i],2.0)-pow(d[i],2.0)+2.0*fabs(Z)*(R[ip1]+R[i]);
    }

    for (i=0;i<m;i++)
        J[i]=2.0*fabs(Z)*atan(Nt[i]/Dt[i])-Y[i]*log((Nl[i])/(Dl[i]));

    for (i=0;i<m;i++)
        (*Src)+=J[i]/(-4.0*pi);

}

void calcDipole_Delhommau(double *Dpl, struct MyVector FPoint, struct MyVector *SPanel, int m)
{
    int i,ip1;
    struct MyVector Area,NormalArea,Center,DummyVector;
    double Z;
    double d[m];
    double Y[m];
    double R[m];
    double Nl[m];
    double Dl[m];
    double Nt[m];
    double Dt[m];
    double J[m];
    double K[m];

    (*Dpl)=0.0;

    if (m==4)
    {
        Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
        Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));
    }
    else
    {
        Area=SVP(0.5,Cross(Minus(SPanel[1],SPanel[0]),Minus(SPanel[2],SPanel[0])));
        Center=SVP(1.0/3.0,Sum3(SPanel[0],SPanel[1],SPanel[2]));
    }

    NormalArea=SVP(1.0/Mag(Area),Area);
    Z=Dot(NormalArea,Minus(FPoint,Center));
//    Z=Dot(NormalArea,Minus(FPoint,Sum(Center,SVP(-1.0e-6,NormalArea))));

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;
        DummyVector=Cross(NormalArea,Minus(SPanel[ip1],SPanel[i]));
        d[i]=Mag(Minus(SPanel[ip1],SPanel[i]));
        Y[i]=Dot(Minus(FPoint,SPanel[i]),DummyVector)/d[i];
        R[i]=Mag(Minus(FPoint,SPanel[i]));
    }

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;

        Nt[i]=2.0*Y[i]*d[i];
        Dt[i]=pow(R[ip1]+R[i],2.0)-pow(d[i],2.0)+2.0*fabs(Z)*(R[ip1]+R[i]);
    }

    for (i=0;i<m;i++)
        K[i]=2.0*sgn(Z)*atan(Nt[i]/Dt[i]);


    for (i=0;i<m;i++)
        (*Dpl)+=K[i]/(4.0*pi);
}

//struct MyVector calcInducedVelocity_LineIntegral(struct MyVector FPoint, struct MyVector *SLine)
//{
//    struct MyVector a1,b1,l,dummyVec;
//    double d,s0,s1,s2,f,g,molifier;
//    double maga1, magb1;
//
//    dummyVec.X=0.0;
//    dummyVec.Y=0.0;
//    dummyVec.Z=0.0;
//
//    a1=Minus(SLine[0],FPoint);
//    b1=Minus(SLine[1],FPoint);
//
//    maga1=Mag(a1);
//    magb1=Mag(b1);
//
//    d=Mag(Minus(b1,a1));
//    l=SVP(1.0/d, Minus(b1,a1));
//    s0=pow(maga1,2.0);
//    s1=Dot(l,a1);
//    f=sqrt(d*d + 2.0*d*s1 + s0);
//    g=fabs(s0 - s1*s1);
//
//    if (sqrt(g)>1.e-7)
//    {
//        molifier=1.0-exp(-0.02*1.e5*g);
//        s2=(maga1*(d+s1)-s1*f)/(f*g*maga1);
//        dummyVec= SVP(molifier*s2/(-4.0*pi),Cross(a1,l));
//    }
//
//    return dummyVec;
//}

struct MyVector calcInducedVelocity_LineIntegral_Intelligent_Filter(struct MyVector FPoint, struct MyVector *SLine, double RF)
{

    struct MyVector a1,b1,l,dummyVec;
    double d,s0,s1,s2,f,g,molifier;
    double maga1, magb1;

    dummyVec.X=0.0;
    dummyVec.Y=0.0;
    dummyVec.Z=0.0;

    a1=Minus(SLine[0],FPoint);
    b1=Minus(SLine[1],FPoint);

    maga1=Mag(a1);
    magb1=Mag(b1);

    d=Mag(Minus(b1,a1));
    l=SVP(1.0/d, Minus(b1,a1));
    s0=pow(maga1,2.0);
    s1=Dot(l,a1);
    f=sqrt(d*d + 2.0*d*s1 + s0);
    g=fabs(s0 - s1*s1);

    if (sqrt(g)>1.0e-6)
    {
//        molifier=1.0-exp(-0.1*g/pow(RF,2.0));
        molifier=1.0-exp(-500000.0*g/pow(RF,2.0));
        s2=(maga1*(d+s1)-s1*f)/(f*g*maga1);
        dummyVec= SVP(molifier*s2/(-4.0*pi),Cross(a1,l));
    }

    return dummyVec;
}

struct MyVector calcSourceUnitVelocity_New(struct MyVector FPoint, struct MyVector *SPanel, int m)
{
    int i,ip1;
    double Z, P1, P2;
    struct MyVector Area, Center, n_nrm, l_nrm, m_nrm, C;

    struct MyVector s[m], r[m];
    double d[m], R[m], Nl[m], Dl[m], Nt[m], Dt[m];
    double SL[m], SM[m], RL[m], RM[m], A[m];

    struct MyVector J[m];
    struct MyVector res;

    res.X=res.Y=res.Z=0.0;

    if (m==4)
    {
        Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));
        Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
        n_nrm=Normalized(Area);
        l_nrm=Normalized(Minus(SVP(0.5,Sum(SPanel[1],SPanel[2])),Center));
        m_nrm=Cross(n_nrm, l_nrm);
    }
    else
    {
        Area=SVP(0.5,Cross(Minus(SPanel[1],SPanel[0]),Minus(SPanel[2],SPanel[0])));
        Center=SVP(1.0/3.0,Sum(SPanel[0],Sum(SPanel[1],SPanel[2])));
        n_nrm=Normalized(Area);
        l_nrm=Normalized(Minus(SVP(0.5,Sum(SPanel[1],SPanel[2])),Center));
        m_nrm=Cross(n_nrm, l_nrm);
    }

    C=Minus(FPoint,Center);
    Z=Dot(n_nrm,C);

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;

        r[i]=Minus(FPoint,SPanel[i]);
        R[i]=Mag(r[i]);
        s[i]=Minus(SPanel[ip1],SPanel[i]);
        d[i]=Mag(s[i]);

        SL[i]=Dot(s[i],l_nrm);
        SM[i]=Dot(s[i],m_nrm);
        RL[i]=Dot(r[i],l_nrm);
        RM[i]=Dot(r[i],m_nrm);

        A[i]=RM[i]*SL[i] - RL[i]*SM[i];
    }

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;

        P1=Z*Z*SL[i] + A[i]*RM[i];
        P2=Z*Z*SL[i] + A[i]*RM[ip1];

        Nl[i]=R[ip1]+R[i]+d[i];
        Dl[i]=R[ip1]+R[i]-d[i]+1.e-6;

        Nt[i]=Z*SM[i]*(R[ip1]*P1 - R[i]*P2);
        Dt[i]=P1*P2 + Z*Z*R[i]*R[ip1]*SM[i]*SM[i]+1.e-6;
    }

    for (i=0;i<m;i++)
    {
        ip1=i+1;
        if(ip1==m)
            ip1=0;

        J[i]=Sum(SVP(log(fabs(Nl[i]/Dl[i]))/d[i] , Minus(SVP(SM[i],l_nrm),SVP(SL[i],m_nrm))),SVP(atan(Nt[i]/Dt[i]),n_nrm));
    }

    for (i=0;i<m;i++)
        res=Sum( res, J[i]);

    return res;
}

struct MyVector calcSourceUnitVelocity(struct MyVector FPoint, struct MyVector *SPanel)
{
    int i;
    int prev_idx;
    double d,z,m,e1,e2,r1,r2,h1,h2,u,v,delta_theta;
    struct MyVector node_a;
    struct MyVector node_b;
    struct MyVector Dummy,res;

    res.X=0.0;res.Y=0.0;res.Z=0.0;

    double MMM[3][3];
    struct MyVector F;
    struct MyVector N[4];

    struct MyVector Center=SVP(0.25,Sum(Sum(SPanel[0],SPanel[1]),Sum(SPanel[2],SPanel[3])));
    struct MyVector Area=SVP(0.5,Cross(Minus(SPanel[2],SPanel[0]),Minus(SPanel[3],SPanel[1])));

    struct MyVector u_Vec=Normalized(Minus(SPanel[1],SPanel[0]));
    struct MyVector n_Vec=Normalized(Area);
    struct MyVector o_Vec=Cross(n_Vec,u_Vec);

    MMM[0][0]=u_Vec.X;
    MMM[0][1]=u_Vec.Y;
    MMM[0][2]=u_Vec.Z;
    MMM[1][0]=o_Vec.X;
    MMM[1][1]=o_Vec.Y;
    MMM[1][2]=o_Vec.Z;
    MMM[2][0]=n_Vec.X;
    MMM[2][1]=n_Vec.Y;
    MMM[2][2]=n_Vec.Z;

    N[0]=MatrixDotVector(MMM,Minus(SPanel[0],Center));
    N[1]=MatrixDotVector(MMM,Minus(SPanel[1],Center));
    N[2]=MatrixDotVector(MMM,Minus(SPanel[2],Center));
    N[3]=MatrixDotVector(MMM,Minus(SPanel[3],Center));

    F=MatrixDotVector(MMM,Minus(FPoint,Center));
    z = F.Z;

    if (fabs(z)<1.e-7)
    {
        for ( i = 0; i < 4; i++)
        {
            if (i == 0)
                prev_idx = 3;
            else
                prev_idx = i - 1;

            node_a = N[prev_idx];
            node_b = N[i];

            d = sqrt(pow(node_b.X - node_a.X, 2.0) + pow(node_b.Y - node_a.Y, 2.0));

            if (d < 1.0e-7)
            {
                Dummy.X=0.0;Dummy.Y=0.0;Dummy.Z=0.0;
                res=Sum(res,Dummy);
            }
            else
            {
                e1 = pow(F.X - node_a.X, 2.0) + pow(z, 2.0);
                e2 = pow(F.X - node_b.X, 2.0) + pow(z, 2.0);
                r1 = sqrt(e1 + pow(F.Y - node_a.Y, 2.0));
                r2 = sqrt(e2 + pow(F.Y - node_b.Y, 2.0));
                Dummy.X=(node_b.Y-node_a.Y)/d * log((r1 + r2 - d) / (r1 + r2 + d));
                Dummy.Y=(node_a.X-node_b.X)/d * log((r1 + r2 - d) / (r1 + r2 + d));
                Dummy.Z=0.0;
                res=Sum(res,Dummy);
            }
        }
    }
    else
    {
        for ( i = 0; i < 4; i++)
        {
            if (i == 0)
                prev_idx = 3;
            else
                prev_idx = i - 1;

            node_a = N[prev_idx];
            node_b = N[i];
            d = sqrt(pow(node_b.X - node_a.X, 2.0) + pow(node_b.Y - node_a.Y, 2.0));

            if (d < 1.0e-7)
            {
                Dummy.X=0.0;Dummy.Y=0.0;Dummy.Z=0.0;
                res=Sum(res,Dummy);
            }
            else
            {
                m = (node_b.Y - node_a.Y) / (node_b.X - node_a.X);
                e1 = pow(F.X - node_a.X, 2.0) + pow(z, 2.0);
                e2 = pow(F.X - node_b.X, 2.0) + pow(z, 2.0);
                r1 = sqrt(e1 + pow(F.Y - node_a.Y, 2.0));
                r2 = sqrt(e2 + pow(F.Y - node_b.Y, 2.0));
                h1 = (F.X - node_a.X) * (F.Y - node_a.Y);
                h2 = (F.X - node_b.X) * (F.Y - node_b.Y);
                u = (m * e1 - h1) / (z * r1);
                v = (m * e2 - h2) / (z * r2);

                if (u==v)
                    delta_theta = 0.0;
                else
                    delta_theta = atan2(u - v, 1.0 + u * v);

                Dummy.X=(node_b.Y-node_a.Y)/d * log((r1 + r2 - d) / (r1 + r2 + d));
                Dummy.Y=(node_a.X-node_b.X)/d * log((r1 + r2 - d) / (r1 + r2 + d));
                Dummy.Z=delta_theta;
                res=Sum(res,Dummy);
            }
        }
    }
    MMM[0][0]=u_Vec.X;
    MMM[0][1]=o_Vec.X;
    MMM[0][2]=n_Vec.X;
    MMM[1][0]=u_Vec.Y;
    MMM[1][1]=o_Vec.Y;
    MMM[1][2]=n_Vec.Y;
    MMM[2][0]=u_Vec.Z;
    MMM[2][1]=o_Vec.Z;
    MMM[2][2]=n_Vec.Z;
    Dummy=MatrixDotVector(MMM,res);
    res=SVP(-1.0/(4.0*pi),Dummy);
    return res;
}

void MakePhiOldAndPhiBarOld(struct Grid *Grids,double **fi,double **fi_old,double **fibar,double **fibar_old,int n)
{
    int i;
    struct MyVector rVec,Vrb,V_ref;
    double VdotR;

    for(i=0;i<n;i++)
    {
        (*fi_old)[i]=(*fi)[i];
        (*fibar_old)[i]=(*fibar)[i];
    }
}

double FindGCoeff(struct Grid *Grids,int PCell,int NCell)
{
    double R[(Grids->Cells[NCell].NodePerCell)+1],NCellAMag;
    int i;
    struct MyVector PCellCen,NCellCen;
    PCellCen=Grids->Cells[PCell].CellCenter;
    NCellCen=Grids->Cells[NCell].CellCenter;
    NCellAMag=Mag(Grids->Cells[NCell].Area);
    R[0]=1.0/Mag(Minus(PCellCen,NCellCen));

//    for (i=0;i<(Grids->Cells[NCell].NodePerCell);i++)
//        R[i+1]=1.0/Mag(Minus(PCellCen,Grids->Cells[NCell].EdgeCoincide[i]));

//    if ((Grids->Cells[NCell].NodePerCell)==3)
//        return (NCellAMag*VDotV3(R,WTRI));
//    else
//        return (NCellAMag*VDotV4(R,WQUA));

    return (NCellAMag*R[0]);
}

double FindHCoeff(struct Grid *Grids,int PCell,int NCell)
{
    double R1[(Grids->Cells[NCell].NodePerCell)+1],
          R2[(Grids->Cells[NCell].NodePerCell)+1],
          R3[(Grids->Cells[NCell].NodePerCell)+1],
          R4[(Grids->Cells[NCell].NodePerCell)+1],
          NCellAMag;
    int i;
    struct MyVector PCellCen,NCellNor,NCellCen;

    PCellCen=Grids->Cells[PCell].CellCenter;
    NCellNor=Normalized(Grids->Cells[NCell].Area);
    NCellCen=Grids->Cells[NCell].CellCenter;
    NCellAMag=Mag(Grids->Cells[NCell].Area);

    R1[0]=Mag(Minus(PCellCen,NCellCen));
    R2[0]=(-1.0*Dot(Minus(PCellCen,NCellCen),NCellNor))/R1[0];
    R3[0]=1.0/pow(R1[0],2);
    R4[0]=R3[0]*R2[0];

//    for (i=0;i<(Grids->Cells[NCell].NodePerCell);i++)
//    {
//        R1[i+1]=1.0/pow(Mag(Minus(PCellCen,Grids->Cells[NCell].EdgeCoincide[i])),3.0);
//        R2[i+1]=Dot(Minus(Grids->Cells[NCell].EdgeCoincide[i],PCellCen),NCellNor);
//        R3[i+1]=R1[i+1]*R2[i+1];
//    }

//    if ((Grids->Cells[NCell].NodePerCell)==3)
//        return (NCellAMag*VDotV3(R3,WTRI));
//    else
//        return (NCellAMag*VDotV4(R3,WQUA));

    return (NCellAMag*R4[0]);
}

void Discretize(struct Grid *Grids,double **phi,double ***miu,struct MatrixCoefficient *a,
                struct MatrixCoefficient *b,double **src,double **Rhs,double dt, int SStep)
{
    struct MyVector PCellNor;
    double **Un;
    int pcell,i,j,k,s,npc,NC,NTE,cnt,cnt1;
    double J,K,sum;
    double A_Source,B_Dipole;
    struct MyVector FieldPoint, SourcePoint;
    struct MyVector *SourcePanel;
    SourcePanel= (struct MyVector*) calloc(4,sizeof(struct MyVector));
    struct MyVector CG;
    struct MyVector rVec,vBodyTotal,vRES;

    struct MyVector Node0,Node1,Node2,Node3;
    struct MyVector *SP;
    struct MyVector DummyVec;
    SP=(struct MyVector*) calloc(4,sizeof(struct MyVector));

//    Un=AllocMatrix(1,Grids->NC,1,1);
    CG=Grids->DynaVariables.MassCenter;
//    GetMassCenter(&CG,Grids);
    NC=Grids->NC;
//    NST=Grids->NST;
    NTE=Grids->NTE;

    for (i=0;i<NC;i++)
        (*Rhs)[i] = 0.0;


    /// Computing free wake influence (RHS)
    for (pcell=0;pcell<NC;pcell++)
    {
        FieldPoint=Grids->Cells[pcell].CellCenter;
        sum=0.0;
        for (i =0;i<NTE;i++)
        {
            for (j=Grids->NST[i];j<SStep*Grids->NST[i] ;j++)
            {
                if (SimilarityOfWakeAndBodyNormals==1)
                {
                    Node0=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[0]].Pos;
                    Node1=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[1]].Pos;
                    Node2=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[2]].Pos;
                    Node3=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[3]].Pos;
                }
                else
                {
                    Node0=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[0]].Pos;
                    Node1=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[3]].Pos;
                    Node2=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[2]].Pos;
                    Node3=Grids->Wakes[i].WNodes[Grids->Wakes[i].WCells[j].NodeList[1]].Pos;
                }

                SP[0]=Node0; SP[1]=Node1; SP[2]=Node2; SP[3]=Node3;
//            calcSourceDipole_Quad_APAME(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//            calcSourceDipole_Quad_Vortexje_Original(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//            calcSourceDipole_Quad_Std(&A_Source,&B_Dipole, FieldPoint, SourcePanel);

                calcDipole_Delhommau(&B_Dipole, FieldPoint, SP, 4);
//                calcSourceDipole_Quad_Vortexje(&A_Source,&B_Dipole, FieldPoint, SP);
                sum+= (*miu)[j][i] * B_Dipole;
            }
        }
        (*Rhs)[pcell]= sum;
    }



    /// Computing matrices of influence coefficients.
    // Influence coefficients between all non-wake surfaces:
    for (i = 0; i < NC; i++)
    {
        FieldPoint=Grids->Cells[i].CellCenter;
        for (j = 0; j < NC; j++)
        {
            PCellNor=Normalized(Grids->Cells[j].Area);
            SourcePoint=Grids->Cells[j].CellCenter;
//            rVec=Minus(SourcePoint,Grids->DynaVariables.MassCenter);
//            vRES=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
            vRES=getPointKinematicVelocity(Grids, FieldPoint);
            vBodyTotal=Minus(vRES,inf_Vel);
            npc=Grids->Cells[j].NodePerCell;
            if (npc==3)
            {
                SourcePanel[0]=Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos;
                SourcePanel[1]=Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos;
                SourcePanel[2]=Grids->Nodes[Grids->Cells[j].NodeList[2]].Pos;
            }
            else
            {
                SourcePanel[0]=Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos;
                SourcePanel[1]=Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos;
                SourcePanel[2]=Grids->Nodes[Grids->Cells[j].NodeList[2]].Pos;
                SourcePanel[3]=Grids->Nodes[Grids->Cells[j].NodeList[3]].Pos;
            }
//            calcSourceDipole_Quad_APAME(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//            calcSourceDipole_Quad_Vortexje_Original(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//            calcSourceDipole_Quad_Std(&A_Source,&B_Dipole, FieldPoint, SourcePanel);

              calcSourceDipole_Delhommau(&A_Source,&B_Dipole, FieldPoint, SourcePanel,npc);
//              calcSourceDipole_Quad_Vortexje(&A_Source,&B_Dipole, FieldPoint, SourcePanel);

            (*Rhs)[i]-=A_Source*Dot(PCellNor, vBodyTotal);

            if (i==j)
                b->Elem[i][j]=0.5;

            else
                b->Elem[i][j]=-B_Dipole;
        }
    }

    /// The influence of the new wake panels:
    int pa,pb;
    for (i=0;i<NC;i++)
    {
        FieldPoint=Grids->Cells[i].CellCenter;
        for (k=0;k<NTE;k++)
            for (j=0;j<Grids->NST[k];j++)
            {
                if (SimilarityOfWakeAndBodyNormals==1)
                {
                    SourcePanel[0]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[0]].Pos;
                    SourcePanel[1]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[1]].Pos;
                    SourcePanel[2]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[2]].Pos;
                    SourcePanel[3]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[3]].Pos;
                }
                else
                {
                    SourcePanel[0]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[0]].Pos;
                    SourcePanel[1]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[3]].Pos;
                    SourcePanel[2]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[2]].Pos;
                    SourcePanel[3]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[1]].Pos;
                }


//            calcSourceDipole_Quad_APAME(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//            calcSourceDipole_Quad_Vortexje_Original(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//            calcSourceDipole_Quad_Std(&A_Source,&B_Dipole, FieldPoint, SourcePanel);

            calcSourceDipole_Delhommau(&A_Source,&B_Dipole, FieldPoint, SourcePanel,npc);
//            calcSourceDipole_Quad_Vortexje(&A_Source,&B_Dipole, FieldPoint, SourcePanel);

                pa = Grids->Wakes[k].UpperCells[j];
                pb = Grids->Wakes[k].LowerCells[j];

                b->Elem[i][pa]-=B_Dipole;
                b->Elem[i][pb]+=B_Dipole;
            }
    }





//    if (SStep==4)
//    {
//        FILE *fid;
//
//        fid=fopen("./A.txt","w");
//        for(i=0;i<NC;i++)
//            for(j=0;j<NC;j++)
//                fprintf(fid,"%f\n",b->Elem[i][j]);
//        fclose(fid);
//
//        fid=fopen("./B.txt","w");
//        for(i=0;i<NC;i++)
//            fprintf(fid,"%f\n",(*Rhs)[i]);
//        fclose(fid);
//    }








//    ///Computing doublet distribution.
//    for (i=0;i<NC;i++)
//        for (j=0;j<NC;j++)
//            (*Rhs)[i]+=a->Elem[i][j]* (*src)[j];

}

//void Discretize1(struct Grid *Grids,double **phi,double ***miu,struct MatrixCoefficient *g,double **Rhs,double dt, int SStep)
//{
//    struct MyVector PCellNor;
//    double **Un;
//    int i,j,k,s,npc,NC,NST,NTE,cnt,cnt1;
//    double J,K,sum;
//    double A_Source,B_Dipole;
//    struct MyVector FieldPoint;
//    struct MyVector *SourcePanel;
//    SourcePanel= (struct MyVector*) malloc(4*sizeof(struct MyVector));
//    struct MyVector CG;
//    struct MyVector rVec,vBodyTotal,vRES;
//
////    Un=AllocMatrix(1,Grids->NC,1,1);
//    CG=Grids->DynaVariables.MassCenter;
////    GetMassCenter(&CG,Grids);
//    NC=Grids->NC;
//    NST=Grids->NST;
//    NTE=Grids->NTE;
//
//    for (i=0;i<(NC+NTE*NST);i++)
//        (*Rhs)[i] = 0.0;
//
//    ///Inflence of Body Elements (Including Source & Dipole) Over Body Elements
//    for (i=0;i<NC;i++)
//        for (j=0;j<NC;j++)
//        {
//            FieldPoint=Grids->Cells[i].CellCenter;
//            npc=Grids->Cells[j].NodePerCell;
//            if (npc==3)
//            {
//                SourcePanel[0]=Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos;
//                SourcePanel[1]=Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos;
//                SourcePanel[2]=Grids->Nodes[Grids->Cells[j].NodeList[2]].Pos;
//            }
//            else
//            {
//                SourcePanel[0]=Grids->Nodes[Grids->Cells[j].NodeList[0]].Pos;
//                SourcePanel[1]=Grids->Nodes[Grids->Cells[j].NodeList[1]].Pos;
//                SourcePanel[2]=Grids->Nodes[Grids->Cells[j].NodeList[2]].Pos;
//                SourcePanel[3]=Grids->Nodes[Grids->Cells[j].NodeList[3]].Pos;
//            }
//
////            calcSourceDipole(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//            calcSourceDipole_Delhommau(&A_Source,&B_Dipole, FieldPoint, SourcePanel,npc);
////            if (i==j)
////                g->Elem[i][j]=1.0-B_Dipole;
////            else
//                g->Elem[i][j]=-B_Dipole;
//
//
//
//            /// //////////////////////////////////////////////////
//            /// //////////////////////////////////////////////////
//            /// //////////////////////////////////////////////////
////            printf("(%d,%d)>>A=%f , B=%f\n",i,j,A_Source,B_Dipole);
////            getchar();
//            /// //////////////////////////////////////////////////
//            /// //////////////////////////////////////////////////
//            /// //////////////////////////////////////////////////
//
//            PCellNor=Normalized(Grids->Cells[j].Area);
//            rVec=Minus(Grids->Cells[j].CellCenter,CG);
//            vBodyTotal=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
//            vRES=Minus(vBodyTotal,inf_Vel);
//            (*Rhs)[i]+=A_Source*Dot(MinusVec(vRES),PCellNor);
//            ///ehteal darad VbodyTotal ra manfi dar nazar begirim!!!!!
////            free(SourcePanel);
//        }
//
//    ///Inflence of Free Vortex Elements Over Body Elements
//    for (i=0;i<NC;i++)
//    {
//        FieldPoint=Grids->Cells[i].CellCenter;
//        for (j=0;j<NTE;j++)
//        {
//            for (k=NST;k<SStep*NST;k++)
//            {
//                SourcePanel[0]=Grids->Wakes[j].WNodes[Grids->Wakes[j].WCells[k].NodeList[0]].Pos;
//                SourcePanel[1]=Grids->Wakes[j].WNodes[Grids->Wakes[j].WCells[k].NodeList[1]].Pos;
//                SourcePanel[2]=Grids->Wakes[j].WNodes[Grids->Wakes[j].WCells[k].NodeList[2]].Pos;
//                SourcePanel[3]=Grids->Wakes[j].WNodes[Grids->Wakes[j].WCells[k].NodeList[3]].Pos;
////                    calcSourceDipole(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//                calcSourceDipole_Delhommau(&A_Source,&B_Dipole, FieldPoint, SourcePanel,4);
//                (*Rhs)[i]+=B_Dipole*(*miu)[k][j];
//            }
//        }
//    }
//
//
//
//    ///Inflence of Kutta Strip Elements Over Body Elements
//    for (i=0;i<NC;i++)
//    {
//        cnt=0;
//        for (k=0;k<NTE;k++)
//            for (j=0;j<NST;j++)
//            {
//                FieldPoint=Grids->Cells[i].CellCenter;
//                SourcePanel[0]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[0]].Pos;
//                SourcePanel[1]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[1]].Pos;
//                SourcePanel[2]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[2]].Pos;
//                SourcePanel[3]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[j].NodeList[3]].Pos;
////                calcSourceDipole(&A_Source,&B_Dipole, FieldPoint, SourcePanel);
//                calcSourceDipole_Delhommau(&A_Source,&B_Dipole, FieldPoint, SourcePanel,4);
//                g->Elem[i][NC+cnt]=-B_Dipole;
//                cnt++;
//            }
//    }
//
//    for (i=0;i<NC;i++)
//        g->Elem[i][i]=1.0;
///// tafavote khotote bala yeki az jalebtarin nokat dar BEM bood. eshare be maghaleye triantafyllou ke
///// dar anja terme dipole ra i=1..Nb , j!=i gosaste karde bood.
//
//
//
//    cnt=0;
//    for (i=0;i<NTE;i++)
//        for (j=0;j<NST;j++)
//        {
//            g->Elem[cnt+NC][(Grids->Wakes[i].UpperCells[j])]=-1.0;
//            g->Elem[cnt+NC][(Grids->Wakes[i].LowerCells[j])]=1.0;
//            g->Elem[cnt+NC][cnt+NC]=1.0;
//            cnt++;
//        }
//
//
//
//
//
//
//
//
////            if (i==j)
////            {
////                G[i][j]= 1.0;
//////                H[i][j]= 3.0*pi;//+ 1.0/dt;
//////                printf("\r%f percent completed...",(double)i *100./(double)NC);
////            }
////            else
////            {
////                FindMatrixCoeff(&J,&K,Grids,i-1,j-1);
////                G[i][j]=J;
////                H[i][j]=K;
////                G[i][j]=FindGCoeff(Grids,i-1,j-1);
////                H[i][j]=FindHCoeff(Grids,i-1,j-1);
////                printf("\r%f percent completed...",(double)i *100./(double)NC);
////            }
//
////    for (i=1;i<=NC;i++)
////    {
//////        H[i][i]+=2.0*pi;
////        PCellNor=Normalized(Grids->Cells[i-1].Area);
////        OR=Minus(Grids->Cells[i-1].CellCenter,CG);
////        Un[i][1]=-Dot(Grids->DynaVariables.LinVel,PCellNor)-Dot(Grids->DynaVariables.AngVel,Cross(OR,PCellNor));;
////        if (fabs(Mag(LinVel))>geps)
////            Un[i][1]+=Dot(LinVel,PCellNor);
////
////        if (fabs(Mag(AngVel))>geps)
////            Un[i][1]+=Dot(AngVel,Cross(OR,PCellNor));
////    }
//
////    MatDotVec(G,Un,RHS,NC);
////    for (i=1;i<=NC;i++)
////        RHS[i][1] += Phi[i][1] / dt;
//
////    FreeMatrix(Un,1,NC,1,1);
//}

void MakeMiuOnKuttaStrip(struct Grid *Grids,double ***miu,double **phi)
{
    int i,j,pa,pb;
    int NTE;

    NTE=Grids->NTE;

//
//    for(i=0;i<NTE;i++)
//        for(j=(SStep-1)*Grids->NST[i]-1;j>=0;j--)
//            (*miu)[j+Grids->NST[i]][i]=(*miu)[j][i];

    for(i=0;i<NTE;i++)
        for(j=0;j<Grids->NST[i];j++)
        {
            pa = Grids->Wakes[i].UpperCells[j];
            pb = Grids->Wakes[i].LowerCells[j];
            (*miu)[j][i]=(*phi)[pa]-(*phi)[pb];
        }

}

void ShiftMiuOnFreeWake(struct Grid *Grids,double ***miu, int SStep)
{
    int i,j,pa,pb;
    int NST,NTE,NC;

    NTE=Grids->NTE;
    for(i=0;i<NTE;i++)
        for(j = SStep * Grids->NST[i] - 1; j >= 0; j--)
            (*miu)[j+Grids->NST[i]][i]=(*miu)[j][i];


}


void CalcPhiBar(struct Grid *Grids,double **fi,double **fibar)
{
    int i;
    struct MyVector rVec,Vrb,V_ref;
    double VdotR;

    for(i=0;i<Grids->NC;i++)
    {
//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        V_ref=Minus(Vrb,inf_Vel);
        VdotR=Dot(V_ref,rVec);
        (*fibar)[i]=-(*fi)[i]+VdotR;
    }
}

void calcMeanVelocityVortexNodes(struct Grid *Grids,double **phi,double ***miu,
                                 double ***iu,double ***iv,double ***iw, int SStep)
{
    int i,j,k,l,m,cnt, strip_number;
    int NTE=Grids->NTE; /// No of Trailing Edges
    //int NST=Grids->NST; /// No of Segment per Trailing Edges
    int NC=Grids->NC; /// No of Elements over Body

    double max_allowable_velocity;
    double mag_induced_velocity;
    struct MyVector normalized_induced_velocity;

    double AMag, sigma, MagQP3;
    struct MyVector QP, QP_SrcVelocity;

    struct MyVector INDUCED_VELOCITY,Bij_Body,Dij_Body,Dij_Kuta,Dij_Vrtx;
    struct MyVector PCellNor, rVec,vBodyTotal,vRES;
    struct MyVector ControlPoint;
    struct MyVector Node0,Node1,Node2,Node3;
    struct MyVector *SP,*SPanel;
    struct MyVector DummyVec;
    double rFilter=0.0;

    SP      =(struct MyVector*) malloc(2*sizeof(struct MyVector));
    SPanel  =(struct MyVector*) malloc(4*sizeof(struct MyVector));


    for (i=0;i<NTE;i++)
    {
        cnt=0;
        for (j=(Grids->NST[i]+1);j<(Grids->NST[i]+1)*(SStep+1);j++)
        {
            ControlPoint=Grids->Wakes[i].WNodes[j].Pos;
            Bij_Body.X=     Bij_Body.Y=     Bij_Body.Z=0.0;
            Dij_Body.X=     Dij_Body.Y=     Dij_Body.Z=0.0;
            Dij_Kuta.X=     Dij_Kuta.Y=     Dij_Kuta.Z=0.0;
            Dij_Vrtx.X=     Dij_Vrtx.Y=     Dij_Vrtx.Z=0.0;

            for (k=0;k<NC;k++)
            {
                rFilter=0.05*D_CH;


//                ///1-calculation of Bij_Body
//                Bij_Body=Sum(Bij_Body,SVP(src[k],calcSourceUnitVelocity(ControlPoint,SPanel)));

                ///1-calculation of Bij_Body
                PCellNor=Normalized(Grids->Cells[k].Area);
                AMag=Mag(Grids->Cells[k].Area);
//                rVec=Minus(Grids->Cells[k].CellCenter,Grids->DynaVariables.MassCenter);
//                vRES=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
                vRES=getPointKinematicVelocity(Grids, Grids->Cells[k].CellCenter);
                vBodyTotal=Minus(vRES,inf_Vel);
                sigma=Dot(vBodyTotal,PCellNor);
                QP=Minus(Grids->Cells[k].CellCenter,ControlPoint);
                MagQP3=pow(Mag(QP),3.0);
                Bij_Body=Sum(Bij_Body , SVP(sigma*AMag/(4.0*pi*MagQP3),QP));

                ///2-calculation of Dij_Body
                if(Grids->Cells[k].NodePerCell==4)
                {
                    Node0=Grids->Nodes[Grids->Cells[k].NodeList[0]].Pos;
                    Node1=Grids->Nodes[Grids->Cells[k].NodeList[1]].Pos;
                    Node2=Grids->Nodes[Grids->Cells[k].NodeList[2]].Pos;
                    Node3=Grids->Nodes[Grids->Cells[k].NodeList[3]].Pos;

//                    SPanel[0]=Node0;    SPanel[1]=Node1;    SPanel[2]=Node2;    SPanel[3]=Node3;
//                    QP_SrcVelocity=calcSourceUnitVelocity_New(ControlPoint,SPanel,4);

                    DummyVec.X=DummyVec.Y=DummyVec.Z=0.0;
                    SP[0]=Node0; SP[1]=Node1;
                    DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP, rFilter));
                    SP[0]=Node1; SP[1]=Node2;
                    DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP, rFilter));
                    SP[0]=Node2; SP[1]=Node3;
                    DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP, rFilter));
                    SP[0]=Node3; SP[1]=Node0;
                    DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP, rFilter));
                }
                else
                {
                    Node0=Grids->Nodes[Grids->Cells[k].NodeList[0]].Pos;
                    Node1=Grids->Nodes[Grids->Cells[k].NodeList[1]].Pos;
                    Node2=Grids->Nodes[Grids->Cells[k].NodeList[2]].Pos;

//                    SPanel[0]=Node0;    SPanel[1]=Node1;    SPanel[2]=Node2;
//                    QP_SrcVelocity=calcSourceUnitVelocity_New(ControlPoint,SPanel,3);

                    DummyVec.X=DummyVec.Y=DummyVec.Z=0.0;
                    SP[0]=Node0; SP[1]=Node1;
                    DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP, rFilter));
                    SP[0]=Node1; SP[1]=Node2;
                    DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP, rFilter));
                    SP[0]=Node2; SP[1]=Node0;
                    DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP, rFilter));
                }
                Dij_Body=Sum(Dij_Body,SVP((*phi)[k],DummyVec));
//                Bij_Body=Sum(Bij_Body , SVP(sigma/(4.0*pi), QP_SrcVelocity));///ok
            }


            ///3-calculation of Dij_Wake
//            #pragma omp parallel for private(l,Node0,Node1,Node2,Node3,DummyVec,*SP) reduction(+:Dij_Kuta) schedule(static, 1)
//            {

                for (k=0;k<NTE;k++)
                {
                    for (l=0;l<SStep*Grids->NST[k];l++)
                    {
                        strip_number=1+ceil(l/Grids->NST[k]);
                        rFilter=2.0*sqrt(Artificial_Visc*strip_number);

                        if (SimilarityOfWakeAndBodyNormals==1)
                        {
                            Node0=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[0]].Pos;
                            Node1=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[1]].Pos;
                            Node2=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[2]].Pos;
                            Node3=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[3]].Pos;
                        }
                        else
                        {
                            Node0=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[0]].Pos;
                            Node1=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[3]].Pos;
                            Node2=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[2]].Pos;
                            Node3=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[1]].Pos;
                        }

                        DummyVec.X=DummyVec.Y=DummyVec.Z=0.0;
                        SP[0]=Node0; SP[1]=Node1;
                        DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP,rFilter));
                        SP[0]=Node1; SP[1]=Node2;
                        DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP,rFilter));
                        SP[0]=Node2; SP[1]=Node3;
                        DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP,rFilter));
                        SP[0]=Node3; SP[1]=Node0;
                        DummyVec=Sum(DummyVec,calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP,rFilter));
                        Dij_Kuta=Sum(Dij_Kuta,SVP((*miu)[l][k],DummyVec));
                    }
                }
//            }

            ///4-ADDING THE EFFECT OF STARTING VORTEX............
            for (k=0;k<NTE;k++)
            {
                    for (l=0;l<=(SStep-1)*Grids->NST[k];l=l+Grids->NST[k])
                    {
                        strip_number=1+ceil(l/Grids->NST[k]);
                        rFilter=2.*sqrt(Artificial_Visc*strip_number);

                        SP[0]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[0]].Pos;
                        SP[1]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[3]].Pos;
                        Dij_Vrtx=Sum(Dij_Vrtx,SVP((*miu)[l][k],calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP,rFilter)));
                    }
                    for (l=(Grids->NST[k]-1);l<=(SStep*Grids->NST[k]-1);l=l+Grids->NST[k])
                    {
                        strip_number=1+ceil(l/Grids->NST[k]);
                        rFilter=2.*sqrt(Artificial_Visc*strip_number);

                        SP[0]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[2]].Pos;
                        SP[1]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[1]].Pos;
                        Dij_Vrtx=Sum(Dij_Vrtx,SVP((*miu)[l][k],calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP,rFilter)));
                    }
                    for (l=(SStep-1)*Grids->NST[k];l<=(SStep*Grids->NST[k]-1);l++)
                    {
                        strip_number=1+ceil(l/Grids->NST[k]);
                        rFilter=2.*sqrt(Artificial_Visc*strip_number);

                        SP[0]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[3]].Pos;
                        SP[1]=Grids->Wakes[k].WNodes[Grids->Wakes[k].WCells[l].NodeList[2]].Pos;
                        Dij_Vrtx=Sum(Dij_Vrtx,SVP((*miu)[l][k],calcInducedVelocity_LineIntegral_Intelligent_Filter(ControlPoint,SP,rFilter)));
                    }
            }

            if (SimilarityOfWakeAndBodyNormals==1)
                INDUCED_VELOCITY=Sum(inf_Vel,Sum(Sum(Bij_Body,Dij_Body),Sum(Dij_Kuta,MinusVec(Dij_Vrtx))));
            else
                INDUCED_VELOCITY=Sum(inf_Vel,Sum(Sum(Bij_Body,Dij_Body),Sum(Dij_Kuta,Dij_Vrtx)));

            if(FilterWakeVelocity==1)
            {
                max_allowable_velocity=1.5*Mag(inf_Vel);
                mag_induced_velocity=Mag(INDUCED_VELOCITY);
                if(mag_induced_velocity>max_allowable_velocity)
                {
                    normalized_induced_velocity=Normalized(INDUCED_VELOCITY);
                    INDUCED_VELOCITY=SVP(max_allowable_velocity,normalized_induced_velocity);
                }
            }

            (*iu)[cnt][i]=INDUCED_VELOCITY.X;
            (*iv)[cnt][i]=INDUCED_VELOCITY.Y;
            (*iw)[cnt][i]=INDUCED_VELOCITY.Z;
            cnt++;
        }
    }
}

void CalcPressure(struct Grid *Grids,double **fibar,double **fibar_old,double **U,double **V,double **W,double **P, int ss)
{
    int i;
    double Diff_fi,t1, t2;
    double V_ref_dot_phi;
    struct MyVector CG,rVec,Vrb,Vref;
    struct MyVector Purt_Vel;

    for (i=0;i<Grids->NC;i++)
    {
        Purt_Vel.X=(*U)[i];Purt_Vel.Y=(*V)[i];Purt_Vel.Z=(*W)[i];

//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        Vref=Minus(inf_Vel,Vrb);

        t1=pow(Mag(Vref),2.0);
        t2=pow(Mag(Sum(Vref,Purt_Vel)),2.0);

        if (ss==1)
            Diff_fi=0.0;
        else
            Diff_fi=((*fibar)[i] - (*fibar_old)[i])/dt;

        (*P)[i]=Pinf+ 0.5*Rho*(t1-t2)- Rho*Diff_fi;
    }
}

void CalcPressureCoefficient_Cp(struct Grid *Grids,double **fibar,double **fibar_old,double **U,double **V,double **W,double **P, int ss)
{
    int i;
    double Diff_fi, t1, t2;
    struct MyVector rVec,Vrb,Vref;
    struct MyVector Purt_Vel;
    struct MyVector n_Vec, Area, VdotN;

    for (i=0;i<Grids->NC;i++)
    {
        Area=Grids->Cells[i].Area;
        n_Vec=Normalized(Area);
        Purt_Vel.X=(*U)[i];Purt_Vel.Y=(*V)[i];Purt_Vel.Z=(*W)[i];

//        rVec=Minus(Grids->Cells[i].CellCenter,Grids->DynaVariables.MassCenter);
//        Vrb=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,rVec));
        Vrb=getPointKinematicVelocity(Grids, Grids->Cells[i].CellCenter);
        Vref=Minus(Vrb, inf_Vel);

        VdotN=MinusVec(SVP(Dot(n_Vec, Vref),n_Vec));

        t1=pow(Mag(Vref),2.0);
        t2=pow(Mag(Sum3(Vref, Purt_Vel, VdotN)),2.0);

        if (ss==1)
            Diff_fi=0.0;
        else
            Diff_fi=((*fibar)[i] - (*fibar_old)[i])/dt;

        (*P)[i]=1.0 - (t2 + 2.0*Diff_fi)/t1;
    }
}

void CalcAddedMassCoeff(struct Grid *Grids,double **Phi,double miu[6])
{
    int i;
    struct MyVector CG,Area,OR,OCA;
    miu[0]=miu[1]=miu[2]=miu[3]=miu[4]=miu[5]=0.0;
//    GetMassCenter(&CG,Grids);
    CG=Grids->DynaVariables.MassCenter;
    for(i=0;i<Grids->NC;i++)
    {
        Area=Grids->Cells[i].Area;
        OR=Minus(Grids->Cells[i].CellCenter,CG);
        OCA=Cross(OR,Area);
        miu[0] += Phi [i+1][1]*(Area.X);
        miu[1] += Phi [i+1][1]*(Area.Y);
        miu[2] += Phi [i+1][1]*(Area.Z);
        miu[3] += Phi [i+1][1]*(OCA.X);
        miu[4] += Phi [i+1][1]*(OCA.Y);
        miu[5] += Phi [i+1][1]*(OCA.Z);
    }
}
