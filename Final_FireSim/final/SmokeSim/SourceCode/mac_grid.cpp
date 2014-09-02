// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>
using namespace std;


// Globals:
MACGrid target;

std::vector<vec3> zeroCell;
std::vector<vec3> accBand;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 


MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
   mPhi = orig.mPhi;
   mN_x = orig.mN_x;
   mN_y = orig.mN_y;
   mN_z = orig.mN_z;
   mY = orig.mY;    
   mStatus = orig.mStatus;
   for(int i=0;i<frameNo;i++) phiFrame[i] = orig.phiFrame[i];
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   
   mPhi = orig.mPhi;
   mN_x = orig.mN_x; 
   mN_y = orig.mN_y;
   mN_z = orig.mN_z; 
   mY = orig.mY;  
   mStatus = orig.mStatus;
   for(int i=0;i<frameNo;i++) phiFrame[i] = orig.phiFrame[i];
   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);
   mPhi.initialize(0.0);

   mN_x.initialize();
   mN_y.initialize();
   mN_z.initialize();
   mY.initialize(0.0);     
   mStatus.initialize();
   initalizePhi();
   for(int i=0;i<frameNo;i++) phiFrame[i].initialize();

   setUpAMatrix();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{   
    // TODO: Set initial values for density, temperature, and velocity.
	
	mT(center[0],center[1],center[2]) = Tmax;
	mD(center[0],center[1],center[2]) = 1.0;
	
	FOR_EACH_CELL{
		double Phix, Phiy, Phiz;
		if(i==0) Phix = (mPhi(i+1,j,k)-mPhi(i,j,k))/theCellSize;
		else if(i==theDim[0]-1) Phix = (mPhi(i,j,k)-mPhi(i-1,j,k))/theCellSize;
		else Phix = (mPhi(i+1,j,k)-mPhi(i-1,j,k))/2/theCellSize;

		if(j==0) Phiy = (mPhi(i,j+1,k)-mPhi(i,j,k))/theCellSize;
		else if(j==theDim[1]-1) Phiy = (mPhi(i,j,k)-mPhi(i,j-1,k))/theCellSize;
		else Phiy = (mPhi(i,j+1,k)-mPhi(i,j-1,k))/2/theCellSize;

		if(k==0) Phiz = (mPhi(i,j,k+1)-mPhi(i,j,k))/theCellSize;
		else if(k==theDim[2]-1) Phiz = (mPhi(i,j,k)-mPhi(i,j,k-1))/theCellSize;
		else Phiz = (mPhi(i,j,k+1)-mPhi(i,j,k-1))/2/theCellSize;

		double Phi_mag = vec3(Phix, Phiy, Phiz).Length();
		//if(Phi_mag == 0.0) Phi_mag = 1e-10;
		mN_x(i,j,k) = Phix/Phi_mag;
		mN_y(i,j,k) = Phiy/Phi_mag;
		mN_z(i,j,k) = Phiz/Phi_mag;
		//if(i==10&&k==5) cout<<j<<": "<<mN_y(i,j,k)<<endl;
		//printf("%lf\n", mPhi(i,j,k)/theCellSize);
		cout<<"Gradient magnitude of Phi: "<<Phi_mag<<endl;
		//cout<<"BEFORE AD: " << mPhi(i,j,k)/theCellSize << endl;
	}
}

void MACGrid::initalizePhi(){
	FOR_EACH_CELL{
		double dist = theCellSize * std::sqrt(((double)i-center[0])*((double)i-center[0])+((double)j-center[1])*((double)j-center[1])+((double)k-center[2])*((double)k-center[2]));   //distance to center
		//printf("%lf\n",dist);
		if(std::abs(dist-Radius)<0.000001) mPhi(i,j,k)=0; 
		else mPhi(i,j,k)=(Radius-dist);
		
		if(dist<Radius) mY(i,j,k)=1.0;
		/*if(j>=20) mPhi(i,j,k) = 1;
		else mPhi(i,j,k) = -1;*/
	}
}
void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.
	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) {

			vec3 posG = vec3(i*theCellSize, (j+0.5)*theCellSize, (k+0.5)*theCellSize);
			double newPhi = mPhi.interpolate(posG);
			vec3 oldPos = posG - getVelocity(posG)*dt;
			double oldPhi = mPhi.interpolate(oldPos);

			if(newPhi>0 && oldPhi<=0 ){
				double uf = mU.interpolate(posG - getVelocity(posG)*dt);
				double nx = mN_x.interpolate(posG - getVelocity(posG)*dt);
				//double ny = mN_y.interpolate(posG - getVelocity(posG)*dt);
				//double nz = mN_z.interpolate(posG - getVelocity(posG)*dt);
				//vec3 n = vec3(nx, ny, nz);
				//double Vf = uf*nx;
				//double Vh = Vf + (rho_f / rho_h -1)*S;
				double uh=(rho_f / rho_h -1)*S*nx+uf;
				//vec3 uh = Vh*n + vec3(uf,0,0) - Vf*n;
				target.mU(i,j,k) = uh;
			}else if(newPhi<=0 && oldPhi>0){
				double uf = mU.interpolate(posG - getVelocity(posG)*dt);
				//double nx = mN_x.interpolate(posG - getVelocity(posG)*dt);
				double nx = mN_x.interpolate(posG - getVelocity(posG)*dt);
				//double nz = mN_z.interpolate(posG - getVelocity(posG)*dt);
				//vec3 n = vec3(nx, ny, nz);
				//double Vf = uf*ny;
				//double Vh = Vf + (rho_f / rho_h -1)*S;
				double uh=-(rho_f / rho_h -1)*S*nx+uf;
				//vec3 uh = Vh*n + vec3(0,uf,0) - Vf*n;
				target.mU(i,j,k) = uh;
			 
			 
			}else{
				target.mU(i,j,k) = mU.interpolate(oldPos);		
				
			}
			 
		 }
	  }
	}

	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 vec3 posG = vec3((i+0.5)*theCellSize, j*theCellSize, (k+0.5)*theCellSize);
			 double newPhi = mPhi.interpolate(posG);
			 vec3 oldPos = posG - getVelocity(posG)*dt;
			 double oldPhi = mPhi.interpolate(oldPos);
			 if(newPhi>0 && oldPhi<=0 ){
				double uf = mV.interpolate(posG - getVelocity(posG)*dt);
				//double nx = mN_x.interpolate(posG - getVelocity(posG)*dt);
				double ny = mN_y.interpolate(posG - getVelocity(posG)*dt);
				//double nz = mN_z.interpolate(posG - getVelocity(posG)*dt);
				//vec3 n = vec3(nx, ny, nz);
				//double Vf = uf*ny;
				//double Vh = Vf + (rho_f / rho_h -1)*S;
				double uh=(rho_f / rho_h -1)*S*ny+uf;
				//vec3 uh = Vh*n + vec3(0,uf,0) - Vf*n;
				target.mV(i,j,k) = uh;
				//if(target.mV(i,j,k)!=0.0) cout<<target.mV(i,j,k)<<" "<<i<<" "<<j<<" "<<k<<"asdffasds"<<endl;
			 }else if(newPhi<=0 && oldPhi>0){
				double uf = mV.interpolate(posG - getVelocity(posG)*dt);
				//double nx = mN_x.interpolate(posG - getVelocity(posG)*dt);
				double ny = mN_y.interpolate(posG - getVelocity(posG)*dt);
				//double nz = mN_z.interpolate(posG - getVelocity(posG)*dt);
				//vec3 n = vec3(nx, ny, nz);
				//double Vf = uf*ny;
				//double Vh = Vf + (rho_f / rho_h -1)*S;
				double uh=-(rho_f / rho_h -1)*S*ny+uf;
				//vec3 uh = Vh*n + vec3(0,uf,0) - Vf*n;
				target.mV(i,j,k) = uh;
			 
			 
			}else{
				target.mV(i,j,k) = mV.interpolate(oldPos);	
				
			}
		 }
	  }
	}

	for(int k = 0; k < theDim[MACGrid::Z]+1; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 vec3 posG = vec3((i+0.5)*theCellSize, (j+0.5)*theCellSize, k*theCellSize);
			 double newPhi = mPhi.interpolate(posG);
			 vec3 oldPos = posG - getVelocity(posG)*dt;
			 double oldPhi = mPhi.interpolate(oldPos);
			 if(newPhi>0 && oldPhi<=0 ){
				double uf = mW.interpolate(posG - getVelocity(posG)*dt);
				//double nx = mN_x.interpolate(posG - getVelocity(posG)*dt);
				//double ny = mN_y.interpolate(posG - getVelocity(posG)*dt);
				double nz = mN_z.interpolate(posG - getVelocity(posG)*dt);
				//vec3 n = vec3(nx, ny, nz);
				//double Vf = uf*nz;
				//double Vh = Vf + (rho_f / rho_h -1)*S;
				double uh=(rho_f / rho_h -1)*S*nz+uf;
				//vec3 uh = Vh*n + vec3(0,0,uf) - Vf*n;
				target.mW(i,j,k) = uh;
			}else if(newPhi<=0 && oldPhi>0){
				double uf = mW.interpolate(posG - getVelocity(posG)*dt);
				//double nx = mN_x.interpolate(posG - getVelocity(posG)*dt);
				double nz = mN_z.interpolate(posG - getVelocity(posG)*dt);
				//double nz = mN_z.interpolate(posG - getVelocity(posG)*dt);
				//vec3 n = vec3(nx, ny, nz);
				//double Vf = uf*ny;
				//double Vh = Vf + (rho_f / rho_h -1)*S;
				double uh=-(rho_f / rho_h -1)*S*nz+uf;
				//vec3 uh = Vh*n + vec3(0,uf,0) - Vf*n;
				target.mW(i,j,k) = uh;
			 
			 
			}else{
				target.mW(i,j,k) = mW.interpolate(oldPos);			
			}

		 }
	  }
	}
    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectY(double dt) {       
	double k=1.0;                                         // the constant in equation 16
    // TODO: Calculate new temp and store in target.
	FOR_EACH_CELL{
		vec3 pos = getCenter(i,j,k);
		vec3 vel = getVelocity(pos);
		vec3 oldPos = pos-vel*dt;
			
		double Ymid = mY.interpolate(oldPos)-k;
		double Ynew = -k*dt+Ymid;
				
		if(mPhi(i,j,k)<=0)
		target.mY(i,j,k) = Ynew;
		else target.mY(i,j,k) = 1.0;
	}
		
    // Then save the result to our object.
    mY = target.mY;
}

void MACGrid::advectTemperature(double dt)         
{
    // TODO: Calculate new temp and store in target.
          FOR_EACH_CELL{
				vec3 pos = getCenter(i,j,k);
				vec3 vel = getVelocity(pos);
				vec3 oldPos = pos-vel*dt;
				//target.mT(i,j,k) = mT.interpolate(oldPos);
				//cout<<mY(9,6,2)<<endl;
				double T=Tignition;
				double elapse = 1-mY(i,j,k);
				//cout<<elapse<<endl;
				if(mPhi(i,j,k)>0) T = Tmax;
			    else  
				{
					/*
					if(elapse<0) T=Tignition;
					if(elapse>=0&&elapse<=1) T = Tignition+(Tmax-Tignition)*elapse*5;
					if(elapse>1&&elapse<2) T=Tmax;
					
					if(elapse>=2) T=mT.interpolate(oldPos)-CoolT*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax);
					*/

					T=mT.interpolate(oldPos)-CoolT*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax);

					/*
					if(mY(i,j,k)>=0.8&&mY(i,j,k)<1) T = -(Tmax-Tignition)*(Y-0.8)*5;
					if(mY(i,j,k)<0.8) T=mT.interpolate(oldPos)-CoolT*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax)*(mT(i,j,k)/Tmax);
					if(mY(i,j,k)>=1)  T=Tignition;
					*/
						
				} 
				target.mT(i,j,k)=T;
				
			    }
	
    // Then save the result to our object.
		 mT = target.mT;
	/*FOR_EACH_CELL {
		cout<<mT(i,j,k)<<endl;
	}*/
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target.
	FOR_EACH_CELL{
		vec3 pos = getCenter(i,j,k);
		vec3 vel = getVelocity(pos);
		vec3 oldPos = pos-vel*dt;
		target.mD(i,j,k) = mD.interpolate(oldPos);
	}
			
    // Then save the result to our object.
    mD = target.mD;
}

void MACGrid::advectPhi(double dt){
	
	FOR_EACH_CELL{	
		vec3 pos = getCenter(i,j,k); 
		vec3 uf = getVelocity(pos);
		vec3 w = uf + S*(vec3(mN_x(i,j,k), mN_y(i,j,k), mN_z(i,j,k)));///////////////////////////////////////////notice
		double Phi_x, Phi_y, Phi_z;
		if(w[0] > 0){
			if(i==0) Phi_x = 0;// mPhi(i,j,k)/theCellSize;//(mPhi(i+1,j,k) - mPhi(i,j,k))/theCellSize;
			else Phi_x = (mPhi(i,j,k) - mPhi(i-1,j,k))/theCellSize;
		}else{
			if(i == theDim[0] - 1) Phi_x = 0;//-mPhi(i,j,k)/theCellSize;//(mPhi(i,j,k) - mPhi(i-1,j,k))/theCellSize;
			Phi_x = (mPhi(i+1,j,k) - mPhi(i,j,k))/theCellSize;
		}
		if(w[1] > 0){
			if(j==0) Phi_y = 0;//mPhi(i,j,k)/theCellSize;//(mPhi(i,j+1,k) - mPhi(i,j,k))/theCellSize;
			else Phi_y = (mPhi(i,j,k) - mPhi(i,j-1,k))/theCellSize;
		}else{
			if(j == theDim[1] - 1) Phi_y = 0;//-mPhi(i,j,k)/theCellSize;//(mPhi(i,j,k) - mPhi(i,j-1,k))/theCellSize;
			else Phi_y = (mPhi(i,j+1,k) - mPhi(i,j,k))/theCellSize;
		}
		if(w[2] > 0){
			if(k==0) Phi_z = 0;//mPhi(i,j,k)/theCellSize;//(mPhi(i,j,k+1) - mPhi(i,j,k))/theCellSize;
			else Phi_z = (mPhi(i,j,k) - mPhi(i,j,k-1))/theCellSize;
		}else{
			if(k == theDim[2] - 1) Phi_z = 0;//-mPhi(i,j,k)/theCellSize;//(mPhi(i,j,k) - mPhi(i,j,k-1))/theCellSize;
			else Phi_z = (mPhi(i,j,k+1) - mPhi(i,j,k))/theCellSize;
		}
		
		target.mPhi(i,j,k) = mPhi(i,j,k) - dt*(w[0]*Phi_x + w[1]*Phi_y + w[2]*Phi_z);

		//cout<<"AFTER AD: " << target.mPhi(i,j,k)/theCellSize << endl;

	}

	mPhi = target.mPhi;
}

void MACGrid::computeSurface(){ 
	// Interpolate to compute the Phi=0 isocontour.
	
	// vector<vec3> surface;
	FOR_EACH_CELL{
		if(mPhi(i,j,k)==0) {    // equal to 0
			if(mStatus(i,j,k)!=1.0){
				mStatus(i,j,k) = 1.0;
				accBand.push_back(vec3(i,j,k));
			}
		}

		if(i!=theDim[0]-1){   // for direction i
			if(mPhi(i,j,k)*mPhi(i+1,j,k)<0){  // negative sign
				if(mStatus(i,j,k)!=1.0){
					mStatus(i,j,k) = 1.0;
					accBand.push_back(vec3(i,j,k));
				}
				if(mStatus(i+1,j,k)!=1.0){
					mStatus(i+1,j,k) = 1.0;
					accBand.push_back(vec3(i+1,j,k));
				}
			}		
		}

		if(j!=theDim[1]-1){   // for direction j
			if(mPhi(i,j,k)*mPhi(i,j+1,k)<0){  // negative sign
				if(mStatus(i,j,k)!=1.0){
					mStatus(i,j,k) = 1.0;
					accBand.push_back(vec3(i,j,k));
				}
				if(mStatus(i,j+1,k)!=1.0){
					mStatus(i,j+1,k) = 1.0;
					accBand.push_back(vec3(i,j+1,k));
				}
			}		
		}

		if(k!=theDim[2]-1){   // for direction k
			if(mPhi(i,j,k)*mPhi(i,j,k+1)<0){  // negative sign
				if(mStatus(i,j,k)!=1.0){
					mStatus(i,j,k) = 1.0;
					accBand.push_back(vec3(i,j,k));
				}
				if(mStatus(i,j,k+1)!=1.0){
					mStatus(i,j,k+1) = 1.0;
					accBand.push_back(vec3(i,j,k+1));
				}
			}		
		}
	}
	//FOR_EACH_CELL{
	//	//cout<<mPhi(i,j,k)/theCellSize<<endl;
	//}
	/*FOR_EACH_CELL{
		cout<<"AFTER SURF: " << mPhi(i,j,k)/theCellSize << endl;
	}*/
				
}

struct tentGrid{  // struct in the heap
	vec3 coord;
	double tentPhi;
	tentGrid(vec3 v, double p) : coord(v), tentPhi(p) {}
};

struct tentGrid_greater_than {  // for sorting the heap
    bool operator()(tentGrid const& a, tentGrid const& b) const {
        return abs(a.tentPhi) > abs(b.tentPhi);
    }
};

int MACGrid::hasAccNeighbor(int i, int j, int k){  // return 0 if no acc neighbor, return 1 if has positive acc nbor, return -1 if has negative acc nbor
	int posCount = 0;
	int negCount = 0;
	if(i!=0) {
		if(mStatus(i-1,j,k)==1.0) {
			if(mPhi(i-1,j,k)>=0) posCount++;
			if(mPhi(i-1,j,k)< 0) negCount++;
		}
	}
	if(i!=theDim[0]-1) {
		if(mStatus(i+1,j,k)==1.0) {
			if(mPhi(i+1,j,k)>=0) posCount++;
			if(mPhi(i+1,j,k)< 0) negCount++;
		}
	}
	if(j!=0) {
		if(mStatus(i,j-1,k)==1.0) {
			if(mPhi(i,j-1,k)>=0) posCount++;
			if(mPhi(i,j-1,k)< 0) negCount++;
		}
	}
	if(j!=theDim[1]-1) {
		if(mStatus(i,j+1,k)==1.0) {
			if(mPhi(i,j+1,k)>=0) posCount++;
			if(mPhi(i,j+1,k)< 0) negCount++;
		}
	}
	if(k!=0) {
		if(mStatus(i,j,k-1)==1.0) {
			if(mPhi(i,j,k-1)>=0) posCount++;
			if(mPhi(i,j,k-1)< 0) negCount++;
		}
	}
	if(k!=theDim[2]-1) {
		if(mStatus(i,j,k+1)==1.0) {
			if(mPhi(i,j,k+1)>=0) posCount++;
			if(mPhi(i,j,k+1)< 0) negCount++;
		}
	}

	if(posCount+negCount==0) return 0;
	else if(posCount>negCount) return 1;
	else return -1;
}

void MACGrid::reinitialization(){
	vector<tentGrid> heap;
	make_heap(heap.begin(), heap.end(), tentGrid_greater_than());
	double tempPhi = 0.0;
	
	// !!!!!INSIDE
	// push tentative band into heap
	FOR_EACH_CELL{
		if(mStatus(i,j,k)==0) { // not accepted or tentative
			if(hasAccNeighbor(i,j,k)==1){  // inside
				tempPhi = computeSinglePhi(i,j,k,true);
				mStatus(i,j,k) = -1.0; // change status to tentative
				mPhi(i,j,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j,k),tempPhi)); 
				push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
	}

	while(heap.size()!=0){ // until all inside cells are marched through
		//cout <<heap.size()<<endl;
		// take min cell from heap, push it to accBand
		tentGrid g = heap.front();
		vec3 index_g = g.coord;
		int i = index_g[0], j = index_g[1], k = index_g[2];
		double phi_g = g.tentPhi;
		pop_heap(heap.begin(),heap.end(), tentGrid_greater_than()); heap.pop_back();
		accBand.push_back(index_g);
		mStatus(i,j,k) = 1.0;

		if(heap.size()==0) break;

		if(i!=0){
			tempPhi = computeSinglePhi(i-1,j,k,true);
			if(mStatus(i-1,j,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i-1,j,k) && tempPhi<it->tentPhi) {
						it->tentPhi = tempPhi;
						mPhi(i-1,j,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
				
			} else if(mStatus(i-1,j,k)==0){ // no value neightbor: compute
				mStatus(i-1,j,k) = -1.0; // change status to tentative
				mPhi(i-1,j,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i-1,j,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(i!=theDim[0]-1){
			tempPhi = computeSinglePhi(i+1,j,k,true);
			if(mStatus(i+1,j,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i+1,j,k) && tempPhi<it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i+1,j,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i+1,j,k)==0){ // no value neightbor: compute
				mStatus(i+1,j,k) = -1.0; // change status to tentative
				mPhi(i+1,j,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i+1,j,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(j!=0){
			tempPhi = computeSinglePhi(i,j-1,k,true);
			if(mStatus(i,j-1,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j-1,k) && tempPhi<it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i,j-1,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j-1,k)==0){ // no value neightbor: compute
				mStatus(i,j-1,k) = -1.0; // change status to tentative
				mPhi(i,j-1,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j-1,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(j!=theDim[1]-1){
			tempPhi = computeSinglePhi(i,j+1,k,true);
			if(mStatus(i,j+1,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j+1,k) && tempPhi<it->tentPhi){
						it->tentPhi = tempPhi;
						mPhi(i,j+1,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j+1,k)==0){ // no value neightbor: compute
				mStatus(i,j+1,k) = -1.0; // change status to tentative
				mPhi(i,j+1,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j+1,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(k!=0){
			tempPhi = computeSinglePhi(i,j,k-1,true);
			if(mStatus(i,j,k-1)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j,k-1) && tempPhi<it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i,j,k-1) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j,k-1)==0){ // no value neightbor: compute
				mStatus(i,j,k-1) = -1.0; // change status to tentative
				mPhi(i,j,k-1) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j,k-1),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(k!=theDim[2]-1){
			tempPhi = computeSinglePhi(i,j,k+1,true);
			if(mStatus(i,j,k+1)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j,k+1) && tempPhi<it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i,j,k+1) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j,k+1)==0){ // no value neightbor: compute
				mStatus(i,j,k+1) = -1.0; // change status to tentative
				mPhi(i,j,k+1) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j,k+1),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
	}

	// !!!!!OUTSIDE
	// push tentative band into heap
	FOR_EACH_CELL{
		if(mStatus(i,j,k)==0) { // not accepted or tentative
			if(hasAccNeighbor(i,j,k)==-1){  // outside
				tempPhi = -computeSinglePhi(i,j,k,false);
				mStatus(i,j,k) = -1.0; // change status to tentative
				mPhi(i,j,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
	}


	while(heap.size()!=0){ // until all inside cells are marched through
		//cout <<heap.size()<<endl;
		// take min cell from heap, push it to accBand
		tentGrid g = heap.front();
		vec3 index_g = g.coord;
		int i = index_g[0], j = index_g[1], k = index_g[2];
		double phi_g = g.tentPhi;
		pop_heap(heap.begin(),heap.end(), tentGrid_greater_than()); heap.pop_back();
		accBand.push_back(index_g);
		mStatus(i,j,k) = 1.0;

		if(heap.size()==0) break;

		if(i!=0){
			tempPhi = -computeSinglePhi(i-1,j,k,false);
			if(mStatus(i-1,j,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i-1,j,k) && tempPhi>it->tentPhi){ 
						it->tentPhi = tempPhi;   // NOTICE tempPhi>it->tentPhi
						mPhi(i-1,j,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i-1,j,k)==0){ // no value neightbor: compute
				mStatus(i-1,j,k) = -1.0; // change status to tentative
				mPhi(i-1,j,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i-1,j,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(i!=theDim[0]-1){
			tempPhi = -computeSinglePhi(i+1,j,k,false);
			if(mStatus(i+1,j,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i+1,j,k) && tempPhi>it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i+1,j,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i+1,j,k)==0){ // no value neightbor: compute
				mStatus(i+1,j,k) = -1.0; // change status to tentative
				mPhi(i+1,j,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i+1,j,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(j!=0){
			tempPhi = -computeSinglePhi(i,j-1,k,false);
			if(mStatus(i,j-1,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j-1,k) && tempPhi>it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i,j-1,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j-1,k)==0){ // no value neightbor: compute
				mStatus(i,j-1,k) = -1.0; // change status to tentative
				mPhi(i,j-1,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j-1,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(j!=theDim[1]-1){
			tempPhi = -computeSinglePhi(i,j+1,k,false);
			if(mStatus(i,j+1,k)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j+1,k) && tempPhi>it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i,j+1,k) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j+1,k)==0){ // no value neightbor: compute
				mStatus(i,j+1,k) = -1.0; // change status to tentative
				mPhi(i,j+1,k) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j+1,k),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(k!=0){
			tempPhi = -computeSinglePhi(i,j,k-1,false);
			if(mStatus(i,j,k-1)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j,k-1) && tempPhi>it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i,j,k-1) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j,k-1)==0){ // no value neightbor: compute
				mStatus(i,j,k-1) = -1.0; // change status to tentative
				mPhi(i,j,k-1) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j,k-1),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
		if(k!=theDim[2]-1){
			tempPhi = -computeSinglePhi(i,j,k+1,false);
			if(mStatus(i,j,k+1)==-1.0){  // tentative neighbor: update			
				for(vector<tentGrid>::iterator it=heap.begin();it!=heap.end();it++){
					if(it->coord==vec3(i,j,k+1) && tempPhi>it->tentPhi){ 
						it->tentPhi = tempPhi;
						mPhi(i,j,k+1) = tempPhi;
					}
				}
				make_heap(heap.begin(), heap.end(), tentGrid_greater_than()); // reconstruct heap
			} else if(mStatus(i,j,k+1)==0){ // no value neightbor: compute
				mStatus(i,j,k+1) = -1.0; // change status to tentative
				mPhi(i,j,k+1) = tempPhi;
				heap.push_back(tentGrid(vec3(i,j,k+1),tempPhi)); push_heap(heap.begin(),heap.end(), tentGrid_greater_than());
			}
		}
	}


	//FOR_EACH_CELL{
	//	cout<<"AFTER REIN: " << i<<","<<j<<","<<k<<","<<mPhi(i,j,k)/theCellSize << endl;
	//}
	
	// remember to empty out accBand
	accBand.clear();
	mStatus.initialize(0.0);

	//if(flag==0) flag = 1;
	//else{
	//	FOR_EACH_CELL{
	//		mPhi(i,j,k) += 0.01;
	//	}
	//	flag = 0;
	//}

	//mPhi.initialize();
	////Px(1,0,1) = 13;
	//mPhi(2,1,0) = 8;
	//mPhi(1,2,0) = 8.2;
	////mStatus(1,0,1) = 1;
	//mStatus(1,2,0) = 1;
	//mStatus(2,1,0) = 1;

	//cout<<computeSinglePhi(1,1,0, true);

}

double MACGrid::computeSinglePhi(int i, int j, int k, bool sign){
	
	GridData status;
	status.initialize();
	status = mStatus;
	GridData Phi;
	Phi.initialize();
	Phi = mPhi;

	double unknownPhi;
	int flag;
	double p1, p2, p3;
	if(sign == true){
		if(i!=0) {
			if(status(i-1,j,k)==1 && Phi(i-1,j,k)<0) status(i-1,j,k) = 0;
		}
		if(i!=theDim[0]-1) {
			if(status(i+1,j,k)==1 && Phi(i+1,j,k)<0) status(i+1,j,k) = 0;
		}
		if(j!=0) {
			if(status(i,j-1,k)==1 && Phi(i,j-1,k)<0) status(i,j-1,k) = 0;
		}
		if(j!=theDim[1]-1) {
			if(status(i,j+1,k)==1 && Phi(i,j+1,k)<0) status(i,j+1,k) = 0;
		}
		if(k!=0) {
			if(status(i,j,k-1)==1 && Phi(i,j,k-1)<0) status(i,j,k-1) = 0;
		}
		if(k!=theDim[2]-1) {
			if(status(i,j,k+1)==1 && Phi(i,j,k+1)<0) status(i,j,k+1) = 0;
		}	
	}else{
		if(i!=0) {
			if(status(i-1,j,k)==1) {
				if(Phi(i-1,j,k)>0) status(i-1,j,k) = 0;
				else Phi(i-1,j,k) = -Phi(i-1,j,k);
			}
		}
		if(i!=theDim[0]-1) {
			if(status(i+1,j,k)==1) {
				if(Phi(i+1,j,k)>0) status(i+1,j,k) = 0;
				else Phi(i+1,j,k) = -Phi(i+1,j,k);
			}
		}
		if(j!=0) {
			if(status(i,j-1,k)==1) {
				if(Phi(i,j-1,k)>0) status(i,j-1,k) = 0;
				else Phi(i,j-1,k) = -Phi(i,j-1,k);
			}
		}
		if(j!=theDim[1]-1) {
			if(status(i,j+1,k)==1) {
				if(Phi(i,j+1,k)>0) status(i,j+1,k) = 0;
				else Phi(i,j+1,k) = -Phi(i,j+1,k);
			}
		}
		if(k!=0) {
			if(status(i,j,k-1)==1) {
				if(Phi(i,j,k-1)>0) status(i,j,k-1) = 0;
				else Phi(i,j,k-1) = -Phi(i,j,k-1);
			}
		}
		if(k!=theDim[2]-1) {
			if(status(i,j,k+1)==1) {
				if(Phi(i,j,k+1)>0) status(i,j,k+1) = 0;
				else Phi(i,j,k+1) = -Phi(i,j,k+1);
			}
		}
	}
	// if sign is +, compute only with the possitive value in neighbors, if sign is - , compute with negative values in neighbors, and convert them into possitive;

		bool xDir, yDir, zDir;
		double Phix, Phiy, Phiz;
		if(i==0){
			xDir = (status(i+1,j,k)==1);
			if(xDir) Phix = Phi(i+1,j,k);
		}else if(i==theDim[0]-1){
			xDir = (status(i-1,j,k)==1);
			if(xDir) Phix = Phi(i-1,j,k);
		}else{
			xDir = (status(i-1,j,k)==1 || status(i+1,j,k)==1);
			if((status(i-1,j,k)==1 && status(i+1,j,k)==1)) Phix = min( Phi(i-1,j,k), Phi(i+1,j,k));
			else Phix = (status(i-1,j,k)==1)?Phi(i-1,j,k):Phi(i+1,j,k);
		}

		if(j==0){
			yDir = (status(i,j+1,k)==1);
			if(yDir) Phiy = Phi(i,j+1,k);
		}else if(j==theDim[1]-1){
			yDir = (status(i,j-1,k)==1);
			if(yDir) Phiy = Phi(i,j-1,k);
		}else{
			yDir = (status(i,j-1,k)==1 || status(i,j+1,k)==1);
			if((status(i,j-1,k)==1 && status(i,j+1,k)==1)) Phiy = min( Phi(i,j-1,k), Phi(i,j+1,k));
			else Phiy = (status(i,j-1,k)==1)? Phi(i,j-1,k): Phi(i,j+1,k);
		}

		if(k==0){
			zDir = (status(i,j,k+1)==1);
			if(zDir) Phiz = Phi(i,j,k+1);
		}else if(k==theDim[2]-1){
			zDir = (status(i,j,k-1)==1);
			if(zDir) Phiz = Phi(i,j,k-1);
		}else{
			zDir = (status(i,j,k-1)==1 || status(i,j,k+1)==1);
			if((status(i,j,k-1)==1 && status(i,j,k+1)==1)) Phiz = min( Phi(i,j,k-1), Phi(i,j,k+1));
			else Phiz = (status(i,j,k-1)==1)? Phi(i,j,k-1): Phi(i,j,k+1);
		}

		if(xDir) assert(Phix >= 0);
		if(yDir) assert(Phiy >= 0);
		if(zDir) assert(Phiz >= 0);

		//cout<<xDir<<" "<<yDir<<" "<<zDir<<endl;
		if(  xDir && yDir && zDir ){  // 3 components in the equation
			//cout<<"3"<<endl;
			/*if(status(i-1,j,k)==1 && status(i+1,j,k)==1){
				Phix = min( Phi(i-1,j,k), Phi(i+1,j,k));
			}else{
				Phix = status(i-1,j,k)==1?  Phi(i-1,j,k):Phi(i+1,j,k);
			}
			if(status(i,j-1,k)==1 && status(i,j+1,k)==1){
				Phiy = min( Phi(i,j-1,k), Phi(i,j+1,k));
			}else{
				Phiy = status(i,j-1,k)==1?  Phi(i,j-1,k):Phi(i,j+1,k);
			}
			if(status(i,j,k-1)==1 && status(i,j,k+1)==1){
				Phiz = min( Phi(i,j,k-1), Phi(i,j,k+1));
			}else{
				Phiz = status(i,j,k-1)==1?  Phi(i,j,k-1):Phi(i,j,k+1);
			}*/
			double Phis, Phim, Phil;
			vector<double> threePhi;
			threePhi.push_back(Phix);
			threePhi.push_back(Phiy);
			threePhi.push_back(Phiz);
			sort(threePhi.begin(), threePhi.end());
			Phil = threePhi[2];
			Phim = threePhi[1];
			Phis = threePhi[0];

			double P = ((Phil-Phix)*(Phil-Phix) + (Phil-Phiy)*(Phil-Phiy) + (Phil-Phiz)*(Phil-Phiz))/theCellSize/theCellSize;
			if(P <= 1) {
				flag = 3;	
				p1 = Phil; p2 = Phim; p3 = Phis;
			}else{
				P = ((Phim-Phim)*(Phim-Phim) + (Phim-Phis)*(Phim-Phis))/theCellSize/theCellSize;
				if(P <= 1){
					flag = 2;
					p1 = Phis; p2 = Phim;
				}else{
					flag = 1;
					p1 = Phis;
				}

			}
				
		}else if( (xDir && yDir && !zDir) || (xDir && !yDir && zDir) || (!xDir && yDir && zDir)){ // 2 components in the equation
			double Phim, Phis;
			double P;
			if(xDir && yDir && !zDir){
				if(Phix > Phiy) {
					Phim = Phix;
					Phis = Phiy;
				}else{
					Phim = Phiy;
					Phis = Phix;
				}
				P = ((Phim-Phix)*(Phim-Phix) + (Phim-Phiy)*(Phim-Phiy))/theCellSize/theCellSize;
			}else if(xDir && !yDir && zDir){
				if(Phix > Phiz) {
					Phim = Phix;
					Phis = Phiz;
				}else{
					Phim = Phiz;
					Phis = Phix;
				}
				P = ((Phim-Phix)*(Phim-Phix) + (Phim-Phiz)*(Phim-Phiz))/theCellSize/theCellSize;
			}else{
				if(Phiz > Phiy) {
					Phim = Phiz;
					Phis = Phiy;
				}else{
					Phim = Phiy;
					Phis = Phiz;
				}
				P = ((Phim-Phiy)*(Phim-Phiy) + (Phim-Phiz)*(Phim-Phiz))/theCellSize/theCellSize;
			}
			//cout<<"2"<<endl;
			if(P <= 1){
				flag = 2;
				p1 = Phis; p2 = Phim;
			}else{
				flag = 1;
				p1 = Phis;
			}
			
		}else {  // only 1 component in the equation
			//cout<<"1"<<endl;
			
			if(xDir && !yDir && !zDir){
				p1 = Phix;
			}else if(!xDir && yDir && !zDir){
				p1 = Phiy;
			}else{
				p1 = Phiz;
			}
			flag = 1;		
		}
		if(flag == 3){
			//cout<<"FLAG 33333"<<endl;
			
			unknownPhi = (2*(p1+p2+p3) + sqrt(4*(p1+p2+p3)*(p1+p2+p3) - 12*(p1*p1 + p2*p2 + p3*p3 - theCellSize*theCellSize)))/6; // larger solution + sign
			//if(i==5&&j==16&&k==4) cout<<i<<" "<<j<<" "<<k<<"FIND IT "<<p1<<","<<p2<<","<<p3<<endl;
			//if(i==4&&j==16&&k==4) cout<<i<<" "<<j<<" "<<k<<"FIND IT "<<p1<<","<<p2<<","<<p3<<endl;
		}else if(flag == 2){
			//cout<<"FLAG 222"<<endl;
			unknownPhi = (2*(p1+p2) + sqrt(4*(p1+p2)*(p1+p2) - 8*(p1*p1 + p2*p2 - theCellSize*theCellSize)))/4; // larger solution + sign
			//if(i==5&&j==16&&k==4) cout<<i<<" "<<j<<" "<<k<<"FIND IT "<<p1<<","<<p2<<endl;
			//if(i==4&&j==16&&k==4) cout<<i<<" "<<j<<" "<<k<<"FIND IT "<<p1<<","<<p2<<endl;
		}else{
			//cout<<"FLAG 1"<<endl;
			unknownPhi = p1 + theCellSize;
			//if(i==5&&j==16&&k==4) cout<<i<<" "<<j<<" "<<k<<"FIND IT "<<p1<<endl;
			//if(i==4&&j==16&&k==4) cout<<i<<" "<<j<<" "<<k<<"FIND IT "<<p1<<endl;
		}
	
	
	//assert(unknownPhi >= 0);
	return unknownPhi;
}


void MACGrid::computeBouyancy(double dt)                  
{
	// TODO: Calculate bouyancy and store in target.
	vec3 Fbuoy;
	double alpha=0,beta=0.1;
	double Tair = 0;
	FOR_EACH_CELL{
		//if(mPhi(i,j,k)>=0) beta=1;       //outer flame beta (bigger)
		//else beta=1;                     //inner gas beta

		if(j==0) Fbuoy = vec3(0, -alpha*(mD(i,j,k))/2.0+beta*(mT(i,j,k)-Tair)/2.0, 0);
		else Fbuoy = vec3(0, -alpha*(mD(i,j,k)+mD(i,j-1,k))/2.0+beta*(mT(i,j,k)+mT(i,j-1,k)-2*Tair)/2.0, 0);
		target.mV(i,j,k) = target.mV(i,j,k) + dt*Fbuoy[1];
	}
		
	// Then save the result to our object.
    mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	// Apply the forces to the current velocity and store the result in target.
	double epsilo=2.0;
	GridData omegaX, omegaY, omegaZ;
	GridData Uc, Vc, Wc;
	omegaX.initialize();
	omegaY.initialize();
	omegaZ.initialize();
	Uc.initialize();
	Vc.initialize();
	Wc.initialize();
	GridData omegaLen;
	omegaLen.initialize();
	vec3 omegaGrad;
	vec3 n;
	GridData FvortX, FvortY, FvortZ;
	FvortX.initialize();
	FvortY.initialize();
	FvortZ.initialize();

	// Calculate the center velocities
	FOR_EACH_CELL{
		vec3 pos = getCenter(i,j,k);
		vec3 vel = getVelocity(pos);
		Uc(i,j,k) = vel[0];
		Vc(i,j,k) = vel[1];
		Wc(i,j,k) = vel[2];
	}
	
	// Calculate omega
	FOR_EACH_CELL{
		double x1,x2;
		double y1,y2;
		double z1,z2;
				
		if(j==0) x1 = Wc(i,j+1,k) - 0;
		else if(j==theDim[1]-1) x1 = 0 - Wc(i,j-1,k);
		else x1 = Wc(i,j+1,k) - Wc(i,j-1,k);

		if(k==0) x2 = Vc(i,j,k+1) - 0;
		else if(k==theDim[2]-1) x2 = 0 - Vc(i,j,k-1);
		else x2 = Vc(i,j,k+1) - Vc(i,j,k-1);

		if(k==0) y1 = Uc(i,j,k+1) - 0;
		else if(k==theDim[2]-1) y1 = 0 - Uc(i,j,k-1);
		else y1 = Uc(i,j,k+1) - Uc(i,j,k-1);

		if(i==0) y2 = Wc(i+1,j,k) - 0;
		else if(i==theDim[0]-1) y2 = 0 - Wc(i-1,j,k);
		else y2 = Wc(i+1,j,k) - Wc(i-1,j,k);

		if(i==0) z1 = Vc(i+1,j,k) - 0;
		else if(i==theDim[0]-1) z1 = 0 - Vc(i-1,j,k);
		else z1 = Vc(i+1,j,k) - Vc(i-1,j,k);

		if(j==0) z2 = Uc(i,j+1,k) - 0;
		else if(j==theDim[1]-1) z2 = 0 - Uc(i,j-1,k);
		else z2 = Uc(i,j+1,k) - Uc(i,j-1,k);

		omegaX(i,j,k) = (x1-x2)/2.0/theCellSize;
		omegaY(i,j,k) = (y1-y2)/2.0/theCellSize;
		omegaZ(i,j,k) = (z1-z2)/2.0/theCellSize;

		omegaLen(i,j,k) = vec3(omegaX(i,j,k), omegaY(i,j,k), omegaZ(i,j,k)).Length();
	}
		
	// Calculate the gradient of |omega|, N, F, and the velocities
	FOR_EACH_CELL{
		double x,y,z;
		if(i==0) x = omegaLen(i+1,j,k)- 0;
		else if(i==theDim[0]-1) x = 0 -  omegaLen(i-1,j,k);
		else x = omegaLen(i+1,j,k) - omegaLen(i-1,j,k);

		if(j==0) y = omegaLen(i,j+1,k) - 0;
		else if(j==theDim[1]-1) y = 0 - omegaLen(i,j-1,k);
		else y = omegaLen(i,j+1,k) - omegaLen(i,j-1,k);

		if(k==0) z = omegaLen(i,j,k+1) - 0;
		else if(k==theDim[2]-1) z = 0 - omegaLen(i,j,k-1);
		else z = omegaLen(i,j,k+1) - omegaLen(i,j,k-1);

		omegaGrad = vec3(x/2.0/theCellSize, y/2.0/theCellSize, z/2.0/theCellSize);

		n = omegaGrad/(omegaGrad.Length() + 1e-20);

		vec3 f = epsilo*theCellSize* n.Cross(vec3(omegaX(i,j,k),omegaY(i,j,k),omegaZ(i,j,k) ));
		FvortX(i,j,k) = f[0];
		FvortY(i,j,k) = f[1];
		FvortZ(i,j,k) = f[2];

		if(i!=0) target.mU(i,j,k) = target.mU(i,j,k) + dt * FvortX.interpolate(vec3( i*theCellSize, (j+0.5)*theCellSize, (k+0.5)*theCellSize));
		if(j!=0) target.mV(i,j,k) = target.mV(i,j,k) + dt * FvortY.interpolate(vec3( (i+0.5)*theCellSize, (j)*theCellSize, (k+0.5)*theCellSize));
		if(k!=0) target.mW(i,j,k) = target.mW(i,j,k) + dt * FvortZ.interpolate(vec3( (i+0.5)*theCellSize, (j+0.5)*theCellSize, (k)*theCellSize));			
	}
		
	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   //computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
	// TODO: Solve Ap = d for pressure.
	// 1. Construct d
	// 2. Construct A
	// 3. Solve for p
	// Subtract pressure from our velocity and save in target.


	// FIRST compute pressure inside
	setUpAMatrix();
	GridData d;
	d.initialize();
	
	FOR_EACH_CELL {
		double U1,U2,V1,V2,W1,W2;
				
		if(i==0) U1=0;
		else U1 = mU(i,j,k);
		if(i==(theDim[0]-1)) U2=0;
		else U2 = mU(i+1,j,k);

		if(j==0) V1=0;
		else V1 = mV(i,j,k);
		if(j==(theDim[1]-1)) V2=0;
		else V2 = mV(i,j+1,k);

		if(k==0) W1=0;
		else W1 = mW(i,j,k);
		if(k==(theDim[2]-1)) W2=0;
		else W2 = mW(i,j,k+1);
				
		d(i,j,k) = (U2-U1+V2-V1+W2-W1)*(-theCellSize)*rho_f /dt;
	}
	
	//cout<<"before first"<<endl;
	conjugateGradient(AMatrix, target.mP, d, 1000, 1e-5);

	//SECOND compute pressure outside using inside surface and equation 3
	GridData mFlag; // used to identify inside cells, surface and outside cells
	mFlag.initialize(-1.0);
	FOR_EACH_CELL{
		if (mPhi(i,j,k)>0) mFlag(i,j,k) = 1.0; 
	}

	FOR_EACH_CELL{
		double curp = 0.0;
		int  count = 0;
		bool surf = false;

		/*double Nx = std::abs(mN_x(i,j,k));
		double Ny = std::abs(mN_y(i,j,k));
		double Nz = std::abs(mN_z(i,j,k));*/

		double Nx = -mN_x(i,j,k);
		double Ny = -mN_y(i,j,k);
		double Nz = -mN_z(i,j,k);


		if(mPhi(i,j,k)<0){
			if(i!=0 && mPhi(i-1,j,k)>0) {
				count++;
				surf = true;

				vec3 pos1 = getCenter(i,j,k);
				double vu = mU.interpolate(pos1);
				double vv = mV.interpolate(pos1);
				double vw = mW.interpolate(pos1);
				double Vh = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));   //dotProduct

				vec3 pos2 = getCenter(i-1,j,k);
				vu = mU.interpolate(pos2);
				vv = mV.interpolate(pos2);
				vw = mW.interpolate(pos2);
				double Vf = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));

				//cout<<curp<<"+";
				double p = rho_f*S*S + target.mP(i-1,j,k) - rho_h*(Vh-Vf+S)*(Vh-Vf+S);	
				curp += p;
				//cout<<p<<"="<<curp<<" ";
			}
			if(i!=theDim[0]-1 && mPhi(i+1,j,k)>0) {
				count++;
				surf = true;

				vec3 pos1 = getCenter(i,j,k);
				double vu = mU.interpolate(pos1);
				double vv = mV.interpolate(pos1);
				double vw = mW.interpolate(pos1);
				double Vh = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));   //dotProduct

				vec3 pos2 = getCenter(i+1,j,k);
				vu = mU.interpolate(pos2);
				vv = mV.interpolate(pos2);
				vw = mW.interpolate(pos2);
				double Vf = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));

				//cout<<curp<<"+";
				double p = rho_f*S*S + target.mP(i+1,j,k) - rho_h*(Vh-Vf+S)*(Vh-Vf+S);
				curp += p;
				//cout<<p<<"="<<curp<<" ";
			}
			if(j!=0 && mPhi(i,j-1,k)>0) {
				count++;
				surf = true;
							
				vec3 pos1 = getCenter(i,j,k);
				double vu = mU.interpolate(pos1);
				double vv = mV.interpolate(pos1);
				double vw = mW.interpolate(pos1);
				double Vh = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));   //dotProduct

				vec3 pos2 = getCenter(i,j-1,k);
				vu = mU.interpolate(pos2);
				vv = mV.interpolate(pos2);
				vw = mW.interpolate(pos2);
				double Vf = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));

				//cout<<curp<<"+";
				double p = rho_f*S*S + target.mP(i,j-1,k) - rho_h*(Vh-Vf+S)*(Vh-Vf+S);
				curp += p;
				//cout<<p<<"="<<curp<<" ";
			}
			if(j!=theDim[1]-1 && mPhi(i,j+1,k)>0) {
				count++;
				surf = true;

				vec3 pos1 = getCenter(i,j,k);
				double vu = mU.interpolate(pos1);
				double vv = mV.interpolate(pos1);
				double vw = mW.interpolate(pos1);
				double Vh = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));   //dotProduct

				vec3 pos2 = getCenter(i,j+1,k);
				vu = mU.interpolate(pos2);
				vv = mV.interpolate(pos2);
				vw = mW.interpolate(pos2);
				double Vf = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));

				//cout<<curp<<"+";
				double p = rho_f*S*S + target.mP(i,j+1,k) - rho_h*(Vh-Vf+S)*(Vh-Vf+S);
				curp += p;
				//cout<<p<<"="<<curp<<" ";
			}
			if(k!=0 && mPhi(i,j,k-1)>0) {
				count++;
				surf = true;

				vec3 pos1 = getCenter(i,j,k);
				double vu = mU.interpolate(pos1);
				double vv = mV.interpolate(pos1);
				double vw = mW.interpolate(pos1);
				double Vh = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));   //dotProduct

				vec3 pos2 = getCenter(i,j,k-1);
				vu = mU.interpolate(pos2);
				vv = mV.interpolate(pos2);
				vw = mW.interpolate(pos2);
				double Vf = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));

				//cout<<curp<<"+";
				double p = rho_f*S*S + target.mP(i,j,k-1) - rho_h*(Vh-Vf+S)*(Vh-Vf+S);
				curp += p;
				//cout<<p<<"="<<curp<<" ";
			}
			if(k!=theDim[2]-1 && mPhi(i,j,k+1)>0) {
				count++;
				surf = true;
							
				vec3 pos1 = getCenter(i,j,k);
				double vu = mU.interpolate(pos1);
				double vv = mV.interpolate(pos1);
				double vw = mW.interpolate(pos1);
				double Vh = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));   //dotProduct

				vec3 pos2 = getCenter(i,j,k+1);
				vu = mU.interpolate(pos2);
				vv = mV.interpolate(pos2);
				vw = mW.interpolate(pos2);
				double Vf = Dot(vec3(vu,vv,vw),vec3(Nx,Ny,Nz));

				//cout<<curp<<"+";
				double p = rho_f*S*S + target.mP(i,j,k+1) - rho_h*(Vh-Vf+S)*(Vh-Vf+S);
				curp += p;
				//cout<<p<<"="<<curp<<" ";
			}
		}

		if(surf==true){
			mFlag(i,j,k) = 0;
			target.mP(i,j,k) = curp/count;
			//cout<<": "<<curp<<"/"<<count<<" = "<<target.mP(i,j,k)<<endl;
		}
	}

	// calculate for outside pressure
	AMatrix.diag.initialize();
	AMatrix.plusI.initialize();
	AMatrix.plusJ.initialize();
	AMatrix.plusK.initialize();
	setUpAMatrixOutside(mFlag);  // modified A matrix
	d.initialize();
	GridData mPout;
	mPout.initialize();

	FOR_EACH_CELL {
		double U1,U2,V1,V2,W1,W2;

		if(mFlag(i,j,k)==1.0) d(i,j,k) = target.mP(i,j,k);
		else if(mFlag(i,j,k)==0.0) d(i,j,k) = target.mP(i,j,k);
		else{
			if(i==0) U1=0;
			else U1 = mU(i,j,k);
			if(i==(theDim[0]-1)) U2=0;	
			else U2 = mU(i+1,j,k);

			if(j==0) V1=0;
			else V1 = mV(i,j,k);
			if(j==(theDim[1]-1)) V2=0;
			else V2 = mV(i,j+1,k);

			if(k==0) W1=0;
			else W1 = mW(i,j,k);
			if(k==(theDim[2]-1)) W2=0;
			else W2 = mW(i,j,k+1);
			
			/*if(i!=0 && mFlag(i-1,j,k)>0) d(i,j,k) += target.mP(i-1,j,k);
			if(i!=theDim[0]-1 && mFlag(i+1,j,k)>0) d(i,j,k) += target.mP(i+1,j,k);
			if(j!=0 && mFlag(i,j-1,k)>0) d(i,j,k) += target.mP(i,j-1,k);
			if(j!=theDim[1]-1 && mFlag(i,j+1,k)>0) d(i,j,k) += target.mP(i,j+1,k);
			if(k!=0 && mFlag(i,j,k-1)>0) d(i,j,k) += target.mP(i,j,k-1);
			if(k!=theDim[2]-1 && mFlag(i,j,k+1)>0) d(i,j,k) += target.mP(i,j,k+1);*/

			d(i,j,k) = (U2-U1+V2-V1+W2-W1)*(-theCellSize)*rho_h /dt;
		}
	}
	//cout<<"before second"<<endl;
	conjugateGradient(AMatrix, mPout, d, 1000, 1e-5);

	FOR_EACH_CELL{
		if(mFlag(i,j,k)==-1.0)
			target.mP(i,j,k) = mPout(i,j,k);
	}

	for(int k = 0; k < theDim[MACGrid::Z]+1; k++) {
		for(int j = 0; j < theDim[MACGrid::Y]+1; j++) {
			for(int i = 0; i < theDim[MACGrid::X]+1; i++) {
			if(i==0||i==theDim[MACGrid::X]) target.mU(i,j,k) = 0;
			else target.mU(i,j,k) -= dt*(target.mP(i,j,k)-target.mP(i-1,j,k))/theCellSize;
			if(j==0||j==theDim[MACGrid::Y]) target.mV(i,j,k) = 0;
			else target.mV(i,j,k) -= dt*(target.mP(i,j,k)-target.mP(i,j-1,k))/theCellSize;
			if(k==0||k==theDim[MACGrid::Z]) target.mW(i,j,k) = 0;
			else target.mW(i,j,k) -= dt*(target.mP(i,j,k)-target.mP(i,j,k-1))/theCellSize;
			}
		}
	}


	// Then save the result to our object
	
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;

	/*FOR_EACH_CELL{
		bool surf = false;
		if(i!=0 && mPhi(i,j,k)*mPhi(i-1,j,k)<0) surf = true;
		if(i!=theDim[0]-1 && mPhi(i,j,k)*mPhi(i+1,j,k)<0) surf = true;
		if(j!=0 && mPhi(i,j,k)*mPhi(i,j-1,k)<0) surf = true;
		if(j!=theDim[1]-1 && mPhi(i,j,k)*mPhi(i,j+1,k)<0) surf = true;
		if(k!=0 && mPhi(i,j,k)*mPhi(i,j,k-1)<0) surf = true;
		if(k!=theDim[2]-1 && mPhi(i,j,k)*mPhi(i,j,k+1)<0) surf = true;
		if(surf){
			vec3 n = -vec3(mN_x(i,j,k),mN_y(i,j,k),mN_z(i,j,k));
			vec3 pos = getCenter(i,j,k);
			double vu = mU.interpolate(pos);
			double vv = mV.interpolate(pos);
			double vw = mW.interpolate(pos);
			double Vh = Dot(vec3(vu,vv,vw),n);
			mU(i,j,k) = mU(i,j,k) + Vh*n[0];
			mU(i+1,j,k) = mU(i+1,j,k) + Vh*n[0];
			mV(i,j,k) = mV(i,j,k) + Vh*n[1];
			mV(i,j+1,k) = mV(i,j+1,k) + Vh*n[1];
			mW(i,j,k) = mW(i,j,k) + Vh*n[2];
			mW(i,j,k+1) = mW(i,j,k+1) + Vh*n[2];
		}
	}*/


	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());

	/*for(int i=0; i< theDim[0]-1; i++){
		for(int j=0; j<theDim[1]-1;j++){
			for(int k=0; k<theDim[2]-1;k++){
				double divergence =  mU(i+1,j,k)-mU(i,j,k)+mV(i,j+1,k)-mV(i,j,k)+mW(i,j,k+1)-mW(i,j,k);
				assert(divergence < 1e-5);
			}
		}
	}*/
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

void MACGrid::setUpAMatrixOutside(GridData mFlag){
	FOR_EACH_CELL {
		if(mFlag(i,j,k)==0.0) AMatrix.diag(i,j,k) = 1;
		else if(mFlag(i,j,k)==1.0) AMatrix.diag(i,j,k) = 1;
		else {
			int numFluidNeighbors = 0;
			if (i-1 >= 0) {
				//if(mFlag(i-1,j,k)!=0.0 && mFlag(i-1,j,k)!=1.0){
					AMatrix.plusI(i-1,j,k) = -1;
					numFluidNeighbors++;
				//}
			}
			if (i+1 < theDim[MACGrid::X]) {
				//if(mFlag(i+1,j,k)!=0.0 && mFlag(i+1,j,k)!=1.0){
					AMatrix.plusI(i,j,k) = -1;
					numFluidNeighbors++;
				//}
			}
			if (j-1 >= 0) {
				//if(mFlag(i,j-1,k)!=0.0 && mFlag(i,j-1,k)!=1.0){
					AMatrix.plusJ(i,j-1,k) = -1;
					numFluidNeighbors++;
				//}
			}
			if (j+1 < theDim[MACGrid::Y]) {
				//if(mFlag(i,j+1,k)!=0.0 && mFlag(i,j+1,k)!=1.0){
					AMatrix.plusJ(i,j,k) = -1;
					numFluidNeighbors++;
				//}
			}
			if (k-1 >= 0) {
				//if(mFlag(i,j,k-1)!=0.0 && mFlag(i,j,k-1)!=1.0){
					AMatrix.plusK(i,j,k-1) = -1;
					numFluidNeighbors++;
				//}
			}
			if (k+1 < theDim[MACGrid::Z]) {
				//if(mFlag(i,j,k+1)!=0.0 && mFlag(i,j,k+1)!=1.0){
					AMatrix.plusK(i,j,k) = -1;
					numFluidNeighbors++;
				//}
			}
			// Set the diagonal:
			AMatrix.diag(i,j,k) = numFluidNeighbors;
		}
	}
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

double MACGrid::getPhi(const vec3& pt)
{
  return mPhi.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}


/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:
	z = r;//(z=Mr)

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		// For now, just bypass the preconditioner:
		z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}





/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
	 

	// accband and outside neighbor
	/*if(mStatus(i,j,k)==1){
		return vec4(1.0,0.0,0.0,1.0);
	}
	else if(mStatus(i,j,k)!=1 && hasAccNeighbor(i,j,k)==-1) return vec4(1.0,1.0,1.0,1.0);
	else return vec4(1.0,0.0,0.0,0.0); */

	// accband and inside neighbor
	/*if(mStatus(i,j,k)==1){
		if(mPhi(i,j,k)>0) return vec4(1.0,0.0,0.0,1.0);
		else return vec4(0.0,0.0,1.0,1.0);
	}
	else if(mStatus(i,j,k)!=1 && hasAccNeighbor(i,j,k)==1) {
		return vec4(1.0,1.0,1.0,1.0);
	}
	else return vec4(1.0,0.0,0.0,0.0);*/

	// accband
	/*if(mStatus(i,j,k)==1){
		return vec4(1.0,0.0,0.0,1.0);
	}
	else return vec4(1.0,0.0,0.0,0.0); */

	// accneighbor
	/*if(mStatus(i,j,k)!=1){
		if(hasAccNeighbor(i,j,k)==1) return vec4(1.0,0.0,0.0,1.0);
		else if(hasAccNeighbor(i,j,k)==-1) return vec4(1.0,1.0,1.0,1.0);
		else return vec4(1.0,0.0,0.0,0.0); 
	}
	else return vec4(1.0,0.0,0.0,0.0); */
	/*if(mStatus(i,j,k)==1) return vec4(1.0,0.0,0.0,1.0);
	else if(mStatus(i,j,k) == -1) return vec4(1.0,1.0,1.0,1.0);
	else return vec4(0.0,0.0,1.0,0.0);

	double t = mT(i,j,k);
	if(t<20) return vec4(0.0,0.0,0.0,value);
	if(t>40) return vec4(1.0,0.0,0.0,value);*/
	
	//if(t<5) return vec4(0.0,0.0,0.0,value);
	/*
	double phi=mPhi(i,j,k);
	if(phi<=0) return vec4(1.0,1.0,1.0,value);
	if(phi>0)  return vec4(1.0,0.0,0.0,value);
    */

	/*if(k==center[2]){
		if(mPhi(i,j,k)>0) return vec4(mPhi(i,j,k)/3.0,0.0,0.2,1.0);
		else return vec4(-mPhi(i,j,k)/20.0,0.2,0.0,1.0);
	}else return vec4(0.0,0.0,1.0,0.0);*/

	if(k==center[2]){
		if(mPhi(i,(j+2.5)/1.5+2,k)>0.4) return vec4(1.0,1.0,1.0,0.5);
		else if(mPhi(i,(j+2.5)/1.5+2,k)>0.2) return vec4(1.0,1.0,0.5,0.5);
		else if(mPhi(i,(j+2.5)/1.6+2,k)>0) return vec4(1.0,0.5,0.25,0.5);
		else if(mPhi(i,(j+2.5)/1.6+2,k)>-0.15) return vec4(1.0,0.1,0.1,0.3);
		else if(mPhi(i,(j+2.5)/1.7+2,k)>-0.3) return vec4(1.0,0.0,0.0,0.1);
		else if(mPhi(i,(j+2.5)/2+2,k)>-0.5) return vec4(1.0,0.0,0.0,0.05);
		else return vec4(0.0,0.0,1.0,0.0);
	}else return vec4(0.0,0.0,1.0,0.0);

	//return vec4(1.0, 0.4, 0.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    /*double value = getDensity(pt); 
	if(pt[2]>theCellSize*center[2]&&pt[2]<theCellSize*(center[2]+1)){
		if(mPhi.interpolate(pt)>0) return vec4(mPhi.interpolate(pt)/3.0,0.0,0.2,1.0);
		else return vec4(-mPhi.interpolate(pt)/20.0,0.2,0.0,1.0);
	}else return vec4(0.0,0.0,1.0,0.0);*/

	if(pt[2]>theCellSize*center[2]&&pt[2]<theCellSize*(center[2]+1)){
		if(mPhi.interpolate(vec3(pt[0],(pt[1]+1.25)/1.5+1,pt[2]))>0.4) return vec4(1.0,1.0,1.0,0.5);
		else if(mPhi.interpolate(vec3(pt[0],(pt[1]+1.25)/1.5+1,pt[2]))>0.2) return vec4(1.0,1.0,0.5,0.5);
		else if(mPhi.interpolate(vec3(pt[0],(pt[1]+1.25)/1.6+1,pt[2]))>0) return vec4(1.0,0.5,0.25,0.5);
		else if(mPhi.interpolate(vec3(pt[0],(pt[1]+1.25)/1.6+1,pt[2]))>-0.15) return vec4(1.0,0.1,0.1,0.3);
		else if(mPhi.interpolate(vec3(pt[0],(pt[1]+1.25)/1.7+1,pt[2]))>-0.3) return vec4(1.0,0.0,0.0,0.1);
		else if(mPhi.interpolate(vec3(pt[0],(pt[1]+1.25)/2+1,pt[2]))>-0.5) return vec4(1.0,0.0,0.0,0.05);
		else return vec4(0.0,0.0,1.0,0.0);
	}else return vec4(0.0,0.0,1.0,0.0);

	/*if(mStatus.interpolate(pt)>0.5) return vec4(1.0,0.0,0.0,1.0);
	else if(mStatus.interpolate(pt)<-0.5) return vec4(1.0,1.0,1.0,1.0);
	else return vec4(0.0,0.0,1.0,0.0);
	
	double t = getTemperature(pt);
	if(t<20) return vec4(0.0,0.0,0.0,value);
	if(t>40) return vec4(1.0,0.0,0.0,value);
    */

	//if(t<5) return vec4(0.0,0.0,0.0,value);
	/*
	double phi=getPhi(pt);
	if(phi<=0) return vec4(1.0,1.0,1.0,value);
	if(phi>0)  return vec4(1.0,0.0,0.0,value);
	*/

	//return vec4(1.0, 0.4, 0.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}


void MACGrid::savePhi(int i){
	phiFrame[i] = mPhi;
}

void MACGrid::setPhi(int i){
	mPhi = phiFrame[i];
}

void MACGrid::phiInterpolate(double r){
	double interval = 1.0/(frameNo-1);
	for(int i=0;i<frameNo;i++){
		if(fabs(r-i*interval)<1.0e-8) {
			mPhi = phiFrame[i];
			//cout<<i<<endl;
			return;
		}
	}
	int index=0;
	while(r>index*interval){
		index++;
	}
	assert(r>(index-1)*interval);
	assert(r<index*interval);
	FOR_EACH_CELL{
		double a = (r-(index-1)*interval)/interval;
		double b = (index*interval-r)/interval;
		mPhi(i,j,k) = b*phiFrame[index-1](i,j,k) + a*phiFrame[index](i,j,k);
	}
	return;
}

void MACGrid::computePhi(double &a, double &l, double &r, int &flag, bool &leftDir){
	if(flag<8) {
		if(flag>=4) savePhi(flag-4);
		if(flag==7) {
			r = 0.95;
			a = 1.0;
			leftDir = false;
		}
	} else{
		if(a>r && leftDir==false){
			double f = (double)rand()/RAND_MAX;
			l = f*r;
			leftDir = true;
		} 
		if(a<l && leftDir==true){
			double f = (double)rand()/RAND_MAX;
			r = l+f*(1.0-l); 
			leftDir = false;
		}

		if(leftDir) a -= 0.05;
		else a += 0.05;
		phiInterpolate(a);
	}
	flag++;
}