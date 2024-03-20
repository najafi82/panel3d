///put coordinates of an arbitrary NACA profile here
///attention: centre {(xte+xle)/2} should be set on origin
double NACA6412[36][2]=
{
    {0.0,-0.00126},
    {0.007596123,0.016007034},
    {0.03015369,0.03591518},
    {0.066987298,0.057298469},
    {0.116977778,0.078159041},
    {0.178606195,0.096325434},
    {0.25,0.109711967},
    {0.328989928,0.116680439},
    {0.413175911,0.116214119},
    {0.5,0.109995113},
    {0.586824089,0.09958919},
    {0.671010072,0.08581977},
    {0.75,0.069773176},
    {0.821393805,0.052695668},
    {0.883022222,0.035953685},
    {0.933012702,0.020957022},
    {0.96984631,0.00904999},
    {0.992403877,0.001389879},
    {1.0,-0.00126},
    {0.992403877,-0.000887097},
    {0.96984631,0.00019003},
    {0.933012702,0.001820803},
    {0.883022222,0.003760818},
    {0.821393805,0.00559671},
    {0.75,0.006874796},
    {0.671010072,0.007174463},
    {0.586824089,0.006259193},
    {0.5,0.004145584},
    {0.413175911,0.001209804},
    {0.328989928,-0.002977889},
    {0.25,-0.009113016},
    {0.178606195,-0.015600682},
    {0.116977778,-0.020750762},
    {0.066987298,-0.022991739},
    {0.03015369,-0.02101839},
    {0.007596123,-0.014012397}
};


double NACA16[30][2]=
{
    {0.125,0.0},
    {0.1146725,0.002305},
    {0.097005,0.00532},
    {0.0758075,0.008625},
    {0.052705,0.01186},
    {0.028775,0.0147925},
    {0.0048325,0.0172275},
    {-0.0184575,0.018985},
    {-0.040535,0.0199025},
    {-0.0609025,0.019855},
    {-0.079125,0.0187725},
    {-0.094795,0.016665},
    {-0.10755,0.013605},
    {-0.1170475,0.00972},
    {-0.1229625,0.00515},
    {-0.125,0.0},
    {-0.1229625,-0.00515},
    {-0.1170475,-0.00972},
    {-0.10755,-0.013605},
    {-0.094795,-0.016665},
    {-0.079125,-0.0187725},
    {-0.0609025,-0.019855},
    {-0.040535,-0.0199025},
    {-0.0184575,-0.018985},
    {0.0048325,-0.0172275},
    {0.028775,-0.0147925},
    {0.052705,-0.01186},
    {0.0758075,-0.008625},
    {0.097005,-0.00532},
    {0.1146725,-0.002305}
};

double NACA10[30][2]=
{
    {0.125,0.0},
    {0.1146725,0.00144},
    {0.097005,0.003325},
    {0.0758075,0.00539},
    {0.052705,0.0074125},
    {0.028775,0.009245},
    {0.0048325,0.0107675},
    {-0.0184575,0.011865},
    {-0.040535,0.01244},
    {-0.0609025,0.01241},
    {-0.079125,0.0117325},
    {-0.094795,0.010415},
    {-0.10755,0.0085025},
    {-0.1170475,0.006075},
    {-0.1229625,0.0032175},
    {-0.125,0.0},
    {-0.1229625,-0.0032175},
    {-0.1170475,-0.006075},
    {-0.10755,-0.0085025},
    {-0.094795,-0.010415},
    {-0.079125,-0.0117325},
    {-0.0609025,-0.01241},
    {-0.040535,-0.01244},
    {-0.0184575,-0.011865},
    {0.0048325,-0.0107675},
    {0.028775,-0.009245},
    {0.052705,-0.0074125},
    {0.0758075,-0.00539},
    {0.097005,-0.003325},
    {0.1146725,-0.00144}
};

const int TRI_NODES[3/*No of Edges*/][3/*Max nodes per Edge + 1*/]=
{
    {2,0,1},//Edge 0
    {2,1,2},//Edge 1
    {2,2,0},//Edge 2
};

const int QUA_NODES[4/*No of Edges*/][3/*Max nodes per Edge + 1*/]=
{
    {2,0,1},//Edge 0
    {2,1,2},//Edge 1
    {2,2,3},//Edge 2
    {2,3,0},//Edge 3
};

const double pi=3.14159265358979323846;

const double eps =1.e-10;

//const int NR_END=1;

//const double WTRI[4]={3./6. , 1./6. , 1./6. , 1./6.};
//const double WQUA[5]={3./7. , 1./7. , 1./7. , 1./7. , 1./7.};
const double WTRI[4]={-27./48. , 25./48 , 25./48. , 25./48.};
const double WQUA[5]={-27./73. , 25./73 , 25./73. , 25./73. , 25./73.};

struct MyVector
{
    double X,Y,Z;
};

struct MyVector LinVel;
struct MyVector AngVel;
double Rho=1.0;
double Visc=0.00001002;
///dynamic viscosity of air @ 20`C =0.0000186, density=1.225
///dynamic viscosity of fresh water @ 20`C =0.0010020, density=1000
double Pinf=0.0; //101325.0;

struct Cell
{
    int Type;  /// Cell Type==> 3:Tri   2:Quad
    int BType; /// Bdry Type==> 1:Body  2:Bottom  3:Farfield  4:Free Surface
    int NodePerCell;
    int *NodeList;   ///Local Numbering
    int *Ngb;  ///Global Neighbours Number
    struct MyVector *EdgeNormal;
    struct MyVector *EdgeCenter;
    struct MyVector CellCenter;
    struct MyVector Area;
    char Side;  /// 'b': Back Side     'f': Face Side
};

struct Node
{
    struct MyVector Pos;
    int NumOfSharedCells;
    int *SharedCells;
};

struct WCell
{
    int *NodeList;
    struct MyVector CellCenter;
    struct MyVector Area;
};

struct WNode
{
    struct MyVector Pos;
    int NumOfSharedCells;
};

//struct Row
//{
//  struct WNode *WNodes;
//  struct WCell *WCells;
//};

struct Wake
{
//    struct Row *Rows;
    struct WNode *WNodes;
    struct WCell *WCells;
    int *UpperCells;
    int *LowerCells;
};

struct DynamicsVariables
{
    double Mass;
    double Inertia[3][3];
    struct MyVector MassCenter,NewMassCenter;
    struct MyVector BodyOr,    NewBodyOr;

    struct MyVector LinAccel,  NewLinAccel;
    struct MyVector AngAccel,  NewAngAccel;
    struct MyVector LinVel,    NewLinVel;
    struct MyVector AngVel,    NewAngVel;
};

struct Zone
{
    int NC;
    int *Cells;
};
struct TECells
{
    int NC;
    int *CellNo;
    int *EdgeNo;
};
struct Grid
{
    int NC;
    int NN;
    int NTE;  /// Number of Trailing Edges
    int *NST;    /// Number of Segments in each Trailing Edge
    struct TECells *TE;
    struct Node *Nodes;   ///Global Numbering
    struct Cell *Cells;
    struct Wake *Wakes;
    struct DynamicsVariables DynaVariables;
    int NZ;
    struct Zone *Zones;
};

struct MatrixCoefficient
{
    double **Elem;
};

int SimStep;

int MAX_SIM_STEP;
int NoIter, timeStepNo=0 ,saveCounter=1;

struct MatrixCoefficient G;
struct MatrixCoefficient B;
struct MatrixCoefficient A;
double *Phi, *Phi_Old, *PhiBar, *PhiBar_Old,*Phi_Nodal,*RHS,*SRC_COEF;
double *BdyPurt_U,*BdyPurt_V,*BdyPurt_W,*Press;
double **VrtxInd_U,**VrtxInd_V,**VrtxInd_W,**Miu;


double error,
        Time,
        StartTime=0.0,
        FinalTime=5.24, //1.97, //4.0,
        dt=0.01818;//0.0109525;//0.001080247;

char InPutFile[80];
struct MyVector Grad;
struct MyVector Gravity={0.0,0.0,0.0};
//struct MyVector inf_Vel={9.0,0.0,0.0};
struct MyVector inf_Vel={-1.0,0.0,0.0};
struct MyVector omega={0.0, 0.0, 0.0};
struct MyVector ref_p={0.0, 0.0, 0.0};
struct Grid Grids;
double Far_Field_Diameter=5.0;
int SimilarityOfWakeAndBodyNormals=0;
int Lifting_Problem=1;
int FilterWakeVelocity=1;
double translateMeshinYdir=0.;
double translateMeshinXdir=0.;

double D_CH;
double Artificial_Visc;


