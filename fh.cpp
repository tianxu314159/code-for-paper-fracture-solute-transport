#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <iostream>



#include <vector>
#include <cmath>




using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::D3Q19Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor
#define ADYNAMICS AdvectionDiffusionBGKdynamics
#define NSDYNAMICS BGKdynamics
//全局定义nx ny nz 
   const plint  nx = 200;
   const plint  ny = 200;
   const plint  nz = 60;
   const T hotTemperature = 0; //in 初始浓度
   const T coldTemperature = 1;//in K,入口浓度
   
   const T dx= 1*1e-5 ;//in ph
   const T dt=1e-6;//in ph
   const T deltap=4.62115309*8*1.00180324584252E-06;//压力梯度800pa/m

   T TAUF= 0.53021;
   T TAUTF=0.56;
   plint temperatureTimeFactor= 1000;
  

  T nsOmega = 1/TAUF;//流动omega
  T adOmegaf = 1/TAUTF;//液体omega
  plint data[nx][ny][nz];



class PressureGradient {
public:
    PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
    { }
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
    {
        velocity.resetToZero();
        density = (T)1 - deltaP*NSDESCRIPTOR<T>::invCs2 / (T)(nx-1) * (T)iX;

    }
private:
    T deltaP;
    plint nx;
};



//几何读取
void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField3D<int>& geometry)
{
    
    Box3D sliceBox(0,0,0,ny-1,0,nz-1);
    std::auto_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
    plb_ifstream geometryFile(fNameIn.c_str());
   for (plint iX=0; iX<=nx-1; iX++) {    
        geometryFile >> *slice;
        copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
    }
    //pcout<<geometry<<endl;
    {
        VtkImageOutput3D<T> vtkOut("geo", 1.0);
        vtkOut.writeData<float>(*copyConvert<int,T>(geometry, geometry.getBoundingBox()), "tag", 1.0);
    }   
};



template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY , plint iZ) const {
        return density;
    }
private:
    T density;
};
//boundary setup

void boundarySetup (
        MultiBlockLattice3D<T, NSDESCRIPTOR>& nsLattice,
        MultiBlockLattice3D<T, ADESCRIPTOR>& adLattice,
        OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>& nsBoundaryCondition,
        OnLatticeAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>& adBoundaryCondition,
        MultiScalarField3D<int>& geometry,
        T deltaP
         )
{
    Array<T,3> zero((T) 0, (T) 0, (T) 0);
    
    Box3D bottom(0,nx-1,0,ny-1,0,0);
    Box3D top(0,nx-1,0,ny-1,nz-1,nz-1);

    Box3D lateral1(0,nx-1,0,0,0,nz-1);
    Box3D lateral2(0,nx-1,ny-1,ny-1,0,nz-1);

    Box3D inlet(0,0,1,ny-2,1,nz-2) ;
    Box3D outlet(nx-1,nx-1,1,ny-2,1,nz-2) ;


    nsBoundaryCondition.addPressureBoundary0N(inlet,nsLattice);
    setBoundaryDensity(nsLattice, inlet, (T) 1.);

    //setAD_NLDboundaryDynamics(adLattice, inlet, boundary::dirichlet);

    pcout<<"inflie the geodata to geofile"<<endl;
    ifstream infile;//定义读取文件流，相对于程序来说是in
    infile.open("geoSMOOTH.dat");//打开文件
    for (plint i=0; i<nx; i++)//x
    {
       // pcout<<"the i is "<<i<<endl;
        for (plint j=0; j<ny; j++)//y
        {
               for (plint k=0;k<nz;k++){//z
               infile >> data[i][j][k];
               }
        }
    }
    infile.close();//读取完成之后关闭文件
   pcout<<"define the inlet and outlet condition!"<<endl;
//--------------------defien the inlet tempreture boundary---------------------------
     for (plint j=0; j<ny; j++) {//y      
               for (plint k=0;k<nz;k++){//z
                 Box3D intp(0,0,j,j,k,k);  
                 Box3D outp(nx-1,nx-1,j,j,k,k);
                if (data[0][j][k]==0)   { 
                setAD_NLDboundaryDynamics(adLattice, intp, boundary::dirichlet);
                setBoundaryDensity(adLattice,intp,coldTemperature);}
               /* else {  
                setAD_NLDboundaryDynamics(adLattice,intp,boundary::neumann);                  
                setBoundaryDensity(adLattice, intp,coldTemperature); }*/
                if (data[nx-1][j][k]==0)
                {
                	nsBoundaryCondition.addPressureBoundary0P(outp,nsLattice);
                    setBoundaryDensity(nsLattice, outp, (T) 1.- deltaP*NSDESCRIPTOR<T>::invCs2);
                    setAD_NLDboundaryDynamics(adLattice, outp, boundary::neumann);
                }
                            }}
 //----------------------------------------------------------------------------------------
    instantiateOuterNLDboundary(adLattice, adLattice.getBoundingBox());   
   // setAD_NLDboundaryDynamics(adLattice, bottom,   boundary::dirichlet);
   // setAD_NLDboundaryDynamics(adLattice, top,      boundary::dirichlet);
    //setAD_NLDboundaryDynamics(adLattice, lateral1, boundary::dirichlet);
   // setAD_NLDboundaryDynamics(adLattice, lateral2, boundary::dirichlet);
   // setBoundaryDensity(adLattice, bottom,  hotTemperature );
    //setBoundaryDensity(adLattice, top,     hotTemperature );
    //setBoundaryDensity(adLattice, lateral1,  hotTemperature );
    //setBoundaryDensity(adLattice, lateral2,  hotTemperature);

   //setAD_NLDboundaryDynamics(adLattice, outlet, boundary::neumann);
   
   
    //integrateProcessingFunctional(new FluidPressureOutlet3D<T,NSDESCRIPTOR,0,+1>(), outlet, nsLattice, 0);
    
    initializeAtEquilibrium( nsLattice, nsLattice.getBoundingBox(), PressureGradient(deltaP, nx) );
    
    initializeAtEquilibrium(adLattice, adLattice.getBoundingBox(), hotTemperature, Array<T,3>((T) 0., (T) 0, (T) 0));
    
    nsLattice.initialize();

    adLattice.initialize();

    pcout<<"define the boundary condition done!"<<endl;
}

void writeVTK(MultiBlockLattice3D<T,NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice3D<T,ADESCRIPTOR>& adLattice,
              plint iter)
{
    VtkImageOutput3D<T> vtkOut(createFileName("TV", iter, 8), dx);
    //vtkOut.writeData<float>(*computeVelocityNorm(nsLattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(nsLattice), "velocity", dx/dt);
    vtkOut.writeData<float>(*computeDensity(adLattice), "temperature", T (0.01));
}

void writeGif(MultiBlockLattice3D<T,NSDESCRIPTOR> nsLattice,
              MultiBlockLattice3D<T,ADESCRIPTOR> adLattice,int iT)
{
    const plint imSize = 600;
   
    Box3D slice(0, nx-1, (ny-1)/2, (ny-1)/2, 0, nz-1);
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("u", iT, 6),
                               *computeVelocityNorm(nsLattice, slice),
                               imSize, imSize);
    // Temperature is the order-0 moment of the advection-diffusion model. It can 
    //    therefore be computed with the function "computeDensity".
    imageWriter.writeScaledGif(createFileName("Concentrition", iT, 6),
                               *computeDensity(adLattice, slice),
                               imSize, imSize);
}

// main loop
int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::IOpolicy().activateParallelIO(true);

    global::directories().setOutputDir("./tmp/");
    pcout<<"reading the geo data"<<endl;
    MultiScalarField3D<int> geometry(nx,ny,nz);
    readGeometry("geoSMOOTH.dat", "tmp/", geometry);

    MultiBlockLattice3D<T, NSDESCRIPTOR> nsLattice (
            nx,ny,nz,new NSDYNAMICS<T, NSDESCRIPTOR>(nsOmega) );
    nsLattice.periodicity().toggleAll(false);  
   // nsLattice.periodicity().toggle(0, true); 

    MultiBlockLattice3D<T, ADESCRIPTOR> adLattice (
            nx,ny,nz,new ADYNAMICS<T, ADESCRIPTOR>(adOmegaf) );

    defineDynamics(nsLattice, geometry, new BounceBack<T,NSDESCRIPTOR>(), 1);//流动矩阵中为1的地方反弹
    defineDynamics(nsLattice, geometry, new NoDynamics<T,NSDESCRIPTOR>(), 2);

   defineDynamics(adLattice, geometry, new BounceBack<T,ADESCRIPTOR>(), 1);
   defineDynamics(adLattice, geometry, new NoDynamics<T,ADESCRIPTOR>(), 2);//1的地方固体改变omega


    adLattice.periodicity().toggleAll(false);  

   // adLattice.periodicity().toggle(1,true);
     //外边界
    OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>*
        nsBoundaryCondition = createLocalBoundaryCondition3D<T,NSDESCRIPTOR>();      
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>*
        adBoundaryCondition = createLocalRegularizedAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>();   

    nsLattice.toggleInternalStatistics(false);
    adLattice.toggleInternalStatistics(false);
    pcout<<"define the BoundaryCondition"<<endl;
   boundarySetup(nsLattice, adLattice,*nsBoundaryCondition, *adBoundaryCondition, geometry, deltap);
//耦合
    integrateProcessingFunctional(
            new LatticeToPassiveAdvDiff3D<T,NSDESCRIPTOR,ADESCRIPTOR>((T) TimeFactor),
            nsLattice.getBoundingBox(), nsLattice, adLattice, 1);

    plint iT = 0;
    plint maxT = 10000000000000;
    plint startIter = 1000;
    plint saveIter = 100;
     // Main loop over time iterations.
    for (iT = 0; iT <= maxT; ++iT) 
    {   
        if (iT % (saveIter*TimeFactor) == 0)
        {
            pcout <<"  Writing pic out at time :" << iT*dt << "s"<< endl;
            writeGif(nsLattice,adLattice,iT);
            writeVTK(nsLattice,adLattice,iT);
        }
         nsLattice.collideAndStream();
        // Lattice Boltzmann iteration step.
         if (iT >=startIter && iT % TimeFactor == 0) {
            adLattice.collideAndStream();
        }
       // adLattice.collideAndStream();
       
    }   
    
}

