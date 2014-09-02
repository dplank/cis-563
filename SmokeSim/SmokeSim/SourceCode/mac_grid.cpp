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


// Globals:
MACGrid target;


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
   mD2 = orig.mD2; // ADDED!
   mT = orig.mT;
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
   mD2 = orig.mD2; // ADDED!
   mT = orig.mT;   

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
   mD2.initialize(); // ADDED!
   mT.initialize(0.0);

   setUpAMatrix();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO: Set initial values for density, temperature, and velocity.

	for(int i=0;i<1;i++){
		for(int j=0;j<theDim[2];j++){
			// First smoke. Red.
			mD(i,0,j) = 1.0;
			mT(i,0,j) = 30.0;
			mV(i,1,j) = 20.0;
			// Second smoke. White.
			mD2(theDim[0]-1-i,0,j) = 1.0;
			mT(theDim[0]-1-i,0,j) = 30.0;
			mV(theDim[0]-1-i,1,j) = 20.0;
		}
	}
}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.
	target.mU = mU;
    target.mV = mV;
    target.mW = mW;

	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) {
			 vec3 posG = vec3(i*theCellSize, (j+0.5)*theCellSize, (k+0.5)*theCellSize);
			 target.mU(i,j,k) = mU.interpolate(posG - getVelocity(posG)*dt);
		 }
	  }
	}

	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 vec3 posG = vec3((i+0.5)*theCellSize, j*theCellSize, (k+0.5)*theCellSize);
			 target.mV(i,j,k) = mV.interpolate(posG - getVelocity(posG)*dt);
		 }
	  }
	}

	for(int k = 0; k < theDim[MACGrid::Z]+1; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 vec3 posG = vec3((i+0.5)*theCellSize, (j+0.5)*theCellSize, k*theCellSize);
			 target.mW(i,j,k) = mW.interpolate(posG - getVelocity(posG)*dt);
		 }
	  }
	}

    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
	target.mT = mT;
	FOR_EACH_CELL{
		vec3 posG = getCenter(i,j,k);
		target.mT(i,j,k) = mT.interpolate(posG - getVelocity(posG)*dt);
	}
    // Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target.
	target.mD = mD;
	FOR_EACH_CELL{
		vec3 posG = getCenter(i,j,k);
		target.mD(i,j,k) = mD.interpolate(posG - getVelocity(posG)*dt);
	}
    // Then save the result to our object.
    mD = target.mD;


	// ADDED! Advect density for second smoke
	target.mD2 = mD2;
	FOR_EACH_CELL{
		vec3 posG = getCenter(i,j,k);
		target.mD2(i,j,k) = mD2.interpolate(posG - getVelocity(posG)*dt);
	}
    // Then save the result to our object.
    mD2 = target.mD2;
}

void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target.
	target.mV = mV;
	double a = 3.0, b = 3.0; // to tune
	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 1; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 double s = (mD(i,j,k) + mD(i,j-1,k) + mD2(i,j,k) + mD2(i,j-1,k))/4; // CHANGED! Original: double s = (mD(i,j,k) + mD(i,j-1,k))/2;
			 double T = (mT(i,j,k)+mT(i,j-1,k))/2;
			 double fbuoy = -a*s+b*T;
			 target.mV(i,j,k) += dt*fbuoy;
		 }
	  }
	}
    // Then save the result to our object.
    mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	// Apply the forces to the current velocity and store the result in target.
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	double eps = 8.0; // to tune

	GridData velu,velv,velw,wu,wv,ww,wlen,fu,fv,fw;
	velu.initialize();
	velv.initialize();
	velw.initialize();
	wu.initialize();
	wv.initialize();
	ww.initialize();
	wlen.initialize();
	fu.initialize();
	fv.initialize();
	fw.initialize();

	FOR_EACH_CELL {
		vec3 pos = getCenter(i,j,k);
		vec3 vel = getVelocity(pos);
		velu(i,j,k) = vel[0];
		velv(i,j,k) = vel[1];
		velw(i,j,k) = vel[2];
    }

	FOR_EACH_CELL {
		double u=0.0, v=0.0, w=0.0;
		if(j!=theDim[1]) { u += velw(i,j+1,k); w -= velu(i,j+1,k);}
		if(j!=0) { u -= velw(i,j-1,k); w += velu(i,j-1,k);}
		if(k!=theDim[2]) { u -= velv(i,j,k+1); v += velu(i,j,k+1);}
		if(k!=0) { u += velu(i,j,k-1); v -= velu(i,j,k-1);}
		if(i!=theDim[0]) { v -= velw(i+1,j,k); w += velv(i+1,j,k);}
		if(i!=0) { v += velw(i-1,j,k); w -= velv(i-1,j,k);}

		wu(i,j,k) = u/(2*theCellSize);
		wv(i,j,k) = v/(2*theCellSize);
		ww(i,j,k) = w/(2*theCellSize);
		vec3 wvec = vec3(wu(i,j,k),wv(i,j,k),ww(i,j,k));
		wlen(i,j,k) = wvec.Length();
	}

	FOR_EACH_CELL {
		double u=0.0, v=0.0, w=0.0;
		if(j!=theDim[1]) v += wlen(i,j+1,k);
		if(j!=0) v -= wlen(i,j-1,k);
		if(k!=theDim[2]) w += wlen(i,j,k+1);
		if(k!=0) w -= wlen(i,j,k-1);
		if(i!=theDim[0]) u += wlen(i+1,j,k);
		if(i!=0) u -= wlen(i-1,j,k);

		double wgradlen=0.0;
		vec3 wgradvec = vec3(u,v,w)/(2*theCellSize);
		wgradlen = wgradvec.Length();

		vec3 N = wgradvec/(wgradlen + 1e-20);
		vec3 wvec = vec3(wu(i,j,k),wv(i,j,k),ww(i,j,k));
		vec3 fvec = eps*theCellSize*(N.Cross(wvec));
		fu(i,j,k) = fvec[0];
		fv(i,j,k) = fvec[1];
		fw(i,j,k) = fvec[2];
	}

	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 1; i < theDim[MACGrid::X]; i++) {
			 target.mU(i,j,k) += dt*(fu(i,j,k)+fu(i-1,j,k))/2;
		 }
	  }
	}

	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 1; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 target.mV(i,j,k) += dt*(fv(i,j,k)+fu(i,j-1,k))/2;
		 }
	  }
	}

	for(int k = 1; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 target.mW(i,j,k) += dt*(fw(i,j,k)+fu(i,j,k-1))/2;
		 }
	  }
	}
	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
	// TODO: Solve Ap = d for pressure.
	// 1. Construct d
	// 2. Construct A
	// 3. Solve for p
	// Subtract pressure from our velocity and save in target.
	target.mP = mP;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	setUpAMatrix(); // not needed
	GridData d;
	d.initialize();
	FOR_EACH_CELL {
		if(i==0) target.mU(i,j,k) = 0;
		if(i==theDim[0]-1) target.mU(i+1,j,k) = 0;
		if(j==0) target.mV(i,j,k) = 0;
		if(j==theDim[1]-1) target.mV(i,j+1,k) = 0;
		if(k==0) target.mW(i,j,k) = 0;
		if(k==theDim[2]-1) target.mW(i,j,k+1) = 0;

		double udiv = target.mU(i+1,j,k) - target.mU(i,j,k) + target.mV(i,j+1,k) - target.mV(i,j,k) + target.mW(i,j,k+1) - target.mW(i,j,k);

		d(i,j,k) = -theCellSize/dt * udiv;		
	}

	conjugateGradient(AMatrix,target.mP,d,1000,0.00001);

	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 1; i < theDim[MACGrid::X]; i++) {
			 target.mU(i,j,k) -= dt*(target.mP(i,j,k)-target.mP(i-1,j,k))/theCellSize;
		 }
	  }
	}

	for(int k = 0; k < theDim[MACGrid::Z]; k++) {
      for(int j = 1; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 target.mV(i,j,k) -= dt*(target.mP(i,j,k)-target.mP(i,j-1,k))/theCellSize;
		 }
	  }
	}

	for(int k = 1; k < theDim[MACGrid::Z]; k++) {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) {
         for(int i = 0; i < theDim[MACGrid::X]; i++) {
			 target.mW(i,j,k) -= dt*(target.mP(i,j,k)-target.mP(i,j,k-1))/theCellSize;
		 }
	  }
	}
	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
	// check divergence
	/*FOR_EACH_CELL {
		double udiv = target.mU(i+1,j,k) - target.mU(i,j,k) + target.mV(i,j+1,k) - target.mV(i,j,k) + target.mW(i,j,k+1) - target.mW(i,j,k);
		printf("%lf\n",udiv/theCellSize);	
	}*/
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

double MACGrid::getDensity2(const vec3& pt) // ADDED!
{
   return mD2.interpolate(pt);
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
	z = r;
	//applyPreconditioner(r,z);

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
		//applyPreconditioner(r,z);

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

void MACGrid::applyPreconditioner(const GridData & r, GridData & z){
	double tau = 0.97;
	GridData precon, q;
	precon.initialize();
	q.initialize();
	FOR_EACH_CELL {
		double Aii = (i==0)?0:AMatrix.plusI(i-1,j,k);
		double Aij = (j==0)?0:AMatrix.plusI(i,j-1,k);
		double Aik = (k==0)?0:AMatrix.plusI(i,j,k-1);
		double Aji = (i==0)?0:AMatrix.plusJ(i-1,j,k);
		double Ajj = (j==0)?0:AMatrix.plusJ(i,j-1,k);
		double Ajk = (k==0)?0:AMatrix.plusJ(i,j,k-1);
		double Aki = (i==0)?0:AMatrix.plusK(i-1,j,k);
		double Akj = (j==0)?0:AMatrix.plusK(i,j-1,k);
		double Akk = (k==0)?0:AMatrix.plusK(i,j,k-1);
		
		double pri = (i==0)?0:precon(i-1,j,k);
		double prj = (j==0)?0:precon(i,j-1,k);
		double prk = (k==0)?0:precon(i,j,k-1);
		double e = AMatrix.diag(i,j,k) - (Aii*pri)*(Aii*pri) - (Ajj*prj)*(Ajj*prj) - (Akk*prk)*(Akk*prk) - tau*(Aii*(Aji+Aki)*pri*pri + Ajj*(Aij+Akj)*prj*prj + Akk*(Aik+Ajk)*prk*prk);
		precon(i,j,k) = 1/sqrt(e + 10e-30);	
	}

	FOR_EACH_CELL {
		double Aii = (i==0)?0:AMatrix.plusI(i-1,j,k);
		double Ajj = (j==0)?0:AMatrix.plusJ(i,j-1,k);
		double Akk = (k==0)?0:AMatrix.plusK(i,j,k-1);
		
		double pri = (i==0)?0:precon(i-1,j,k);
		double prj = (j==0)?0:precon(i,j-1,k);
		double prk = (k==0)?0:precon(i,j,k-1);

		double qi = (i==0)?0:q(i-1,j,k);
		double qj = (j==0)?0:q(i,j-1,k);
		double qk = (k==0)?0:q(i,j,k-1);
		double t = r(i,j,k) - Aii*pri*qi - Ajj*prj*qj - Akk*prk*qk;
		q(i,j,k) = t * precon(i,j,k);	
	}

	FOR_EACH_CELL_REVERSE {
		double zi = (i==theDim[0]-1)?0:z(i+1,j,k);
		double zj = (j==theDim[1]-1)?0:q(i,j+1,k);
		double zk = (k==theDim[2]-1)?0:q(i,j,k+1);
		double t = q(i,j,k) - precon(i,j,k) * (AMatrix.plusI(i,j,k)*zi - AMatrix.plusJ(i,j,k)*zj - AMatrix.plusK(i,j,k)*zk);
		z(i,j,k) = t * precon(i,j,k);
	}


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
    return vec4(1.0, 0.0, 0.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(1.0, 0.0, 0.0, value);

}

vec4 MACGrid::getRenderColor2(int i, int j, int k) // ADDED!
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD2(i, j, k); 
    return vec4(1.0, 1.0, 1.0, value);

}

vec4 MACGrid::getRenderColor2(const vec3& pt) // ADDED!
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity2(pt); 
    return vec4(1.0, 1.0, 1.0, value);

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

			/* ADDED! */
			vec4 color1_ = getRenderColor2(pos1);
            vec4 color2_ = getRenderColor2(pos2);

			double v = 1/(color1[3] + color1_[3]);
			color1 = v*(color1[3]*color1 + color1_[3]*color1_);
			color1[3] = 1/v;

			v = 1/(color2[3] + color2_[3]);
			color2 = v*(color2[3]*color2 + color2_[3]*color2_);
			color2[3] = 1/v;
			/* ADDED! END*/
			
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
			
			/* ADDED! */
			vec4 color1_ = getRenderColor2(pos1);
            vec4 color2_ = getRenderColor2(pos2);

			double v = 1/(color1[3] + color1_[3]);
			color1 = v*(color1[3]*color1 + color1_[3]*color1_);
			color1[3] = 1/v;

			v = 1/(color2[3] + color2_[3]);
			color2 = v*(color2[3]*color2 + color2_[3]*color2_);
			color2[3] = 1/v;
			/* ADDED! END*/

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
			
			/* ADDED! */
			vec4 color1_ = getRenderColor2(pos1);
            vec4 color2_ = getRenderColor2(pos2);

			double v = 1/(color1[3] + color1_[3]);
			color1 = v*(color1[3]*color1 + color1_[3]*color1_);
			color1[3] = 1/v;

			v = 1/(color2[3] + color2_[3]);
			color2 = v*(color2[3]*color2 + color2_[3]*color2_);
			color2[3] = 1/v;
			/* ADDED! END*/

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
			
			/* ADDED! */
			vec4 color1_ = getRenderColor2(pos1);
            vec4 color2_ = getRenderColor2(pos2);

			double v = 1/(color1[3] + color1_[3]);
			color1 = v*(color1[3]*color1 + color1_[3]*color1_);
			color1[3] = 1/v;

			v = 1/(color2[3] + color2_[3]);
			color2 = v*(color2[3]*color2 + color2_[3]*color2_);
			color2[3] = 1/v;
			/* ADDED! END*/

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

	  /* ADDED! */
	  vec4 color = getRenderColor(i,j,k);
	  vec4 color_ = getRenderColor2(i,j,k);
	  double v = 1/(color[3] + color_[3]);
	  color = v*(color[3]*color + color_[3]*color_);
      cube.color = color;
	  /* ADDED! END*/

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
