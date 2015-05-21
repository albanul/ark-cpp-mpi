#include <stdio.h>
#include <cmath>
#include <string.h>
#include "../headers/constants.h"
#include "../headers/bulk.h"
#include "../headers/forces.h"
#include "../headers/stress.h"
#include "../headers/BoundaryConditions.h"
#include "../headers/parallel.h"
#include "../headers/Arr2d.h"

void Input();
void InitializeData();
void TimeStepSize();
void Phase1();
void Phase2();
void StressTensor();
void UseForces();
void WriteDataParaView();
void FreeMemory();
void WriteEnergy();

double QCriterion(int i, int j, int k);
void AverageValues();
void WriteAverageValues(double);

double min3d(double, double, double);
double max3d(double, double, double);

static void swap8(void *);

void wr(int n) {
	if (rank == 0) {
		char filename[100];
		sprintf(filename, "%sout_%d.txt", dirPath, n);

		FILE *fd = fopen(filename, "w");

		for (int i = 0; i < n1 + 1; ++i)
		{
			for (int j = 0; j < n2 + 1; ++j)
			{
				fprintf(fd, "%10.6f", ronCon->elem(i, j, 1));
			}
			fprintf(fd, "\n");
		}
		fprintf(fd, "\n");

		fclose(fd);
	}
}

int main(int argc, char** argv) {
	dims[0] = 0; dims[1] = 0; dims[2] = 1;
	x1Period = false; x2Period = false; x3Period = false;
	int periods[3] = {x1Period ? 1 : 0, x2Period ? 1 : 0, x3Period ? 1 : 0};

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Dims_create(nproc, 3, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &newComm);
	MPI_Comm_rank(newComm, &rank);
	MPI_Cart_coords(newComm, rank, 3, coords);

	MPI_Cart_shift(newComm, 0, 1, &left, &right);
	MPI_Cart_shift(newComm, 1, 1, &down, &up);
	MPI_Cart_shift(newComm, 2, 1, &bottom, &top);

	double startTime = 0.;

	double time = MPI_Wtime();

	Input();
	InitializeData();

	//first output
	WriteDataParaView();
	// WriteEnergy();

	do
	{
		TimeStepSize();
		nStep++;
		TIME += 0.5 * dt;
		Phase1();
		StressTensor();
		UseForces();
		Phase2();
		Phase1();
		UseForces();
		TIME += 0.5 * dt;

		// if (TIME >= 60.) {
		// 	if (startTime < 60.) startTime = TIME;
		// 	AverageValues();
		// }

		if (nStep >= 70000 && nStep % nPrint == 0)
		{
			WriteDataParaView();
			// WriteEnergy();
			if (rank == 0) printf("step: %d dt:%E time:%8.4f\n", nStep, dt, TIME);
		}
	} while (nStep < nStop);

	// if (TIME > 60.0) WriteAverageValues(startTime);

	time = MPI_Wtime() - time;
	double maxTime;
	MPI_Allreduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, newComm);
	if (rank == 0) printf("Program time: %8.3f\n", maxTime);

	FreeMemory();
	MPI_Finalize();
	return 0;
}

void Input() {
	// Little and Big Endian
	needSwap = true;

	// Use Tecplot or ParaView for output
	useTecplot = false;

	// geometry index
	l = 1;

	// Project directory path
	dirPath = "../out/";

	// total number of grid nodes along the x1 axis
	n1_g = 201;
	// total number of grid nodes along the X2 axis
	n2_g = 61;
	// total number of grid nodes along the X3 axis
	n3_g = 2;

	px = dims[0];
	py = dims[1];
	pz = dims[2];

	coordX = coords[0];
	coordY = coords[1];
	coordZ = coords[2];

	// number of grid nodes along the x1 axis to 1 processor
	n1 = (n1_g - 1) / px + 1;
	// number of grid nodes along the X2 axis to 1 processor
	n2 = (n2_g - 1) / py + 1;
	// number of grid nodes along the X3 axis to 1 processor
	n3 = (n3_g - 1) / pz + 1;

	// coordinates of west plane
	x1_w = 0.;
	// coordinates of east plane
	x1_e = 20.;

	// coordinates of south plane
	x2_s = 0.;
	// coordinates of north plane
	x2_n = 6.;

	// coordinates of bottom plane
	x3_b = 0.;
	// coordinates of top plane
	x3_t = 0.1;

	// total number of steps
	nStop = 100000;
	// print interval
	nPrint = 1000;

	// Courant number
	CFL = 0.2;

	// kinematic viscosity
	VIS = 1. / 1000.; // 0.5
	// initial temperature
	t0 = 1.;

	// initial velocity along the x1 axis
	u10 = 1.;
	// initial velocity along the X2 axis
	u20 = 0.;
	// initial velocity along the X3 axis
	u30 = 0.;
	// sound velocity
	sound = 10.;

	// unperturbed density of the liquid
	ro0_g = 1.;

	// #####################################################
	//              block of arrays allocation
	// #####################################################

	// coordinates of grid nodes along the all of the axises
	x1 = new double[n1 + 2];
	x2 = new double[n2 + 2];
	x3 = new double[n3 + 2];

	// variables on the current time step
	roCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	tCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	u1Con = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	u2Con = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	u3Con = new Arr3d(n1 + 1, n2 + 1, n3 + 1);

	// variables on the next time step
	ronCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	tnCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	u1nCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	u2nCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	u3nCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);

	// forces
	f1 = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	f2 = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
	f3 = new Arr3d(n1 + 1, n2 + 1, n3 + 1);

	// solid/liqud condition
	condition = new Arr3d(n1 + 1, n2 + 1, n3 + 1);

	// average values
	averageU1 = new Arr3d(n1 + 1, n2 + 2, n3 + 1);
	averageU2 = new Arr3d(n1 + 1, n2 + 2, n3 + 1);
	averageU3 = new Arr3d(n1 + 1, n2 + 2, n3 + 1);

	// variables perpendicular to the axis x1
	ro1 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	t1 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u11 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u21 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u31 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	p1 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

	// variables perpendicular to the axis X2
	ro2 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	t2 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u12 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u22 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u32 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	p2 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

	// variables perpendicular to the axis X3
	ro3 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	t3 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u13 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u23 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	u33 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	p3 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

	// NMAX = MAX(n1, n2, n3)
	int nmax = (int)max3d(n1, n2, n3);

	// additional buffers for phase 2
	rBuf    = new double[nmax + 2];
	qBuf    = new double[nmax + 1];
	tfBuf   = new double[nmax + 2];
	tbBuf   = new double[nmax + 2];
	u2fBuf  = new double[nmax + 2];
	u2bBuf  = new double[nmax + 2];
	u3fBuf  = new double[nmax + 2];
	u3bBuf  = new double[nmax + 2];

	// friction stress
	sigm11 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	sigm21 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	sigm31 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

	sigm12 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	sigm22 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	sigm32 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

	sigm13 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	sigm23 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
	sigm33 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
}

void InitializeData() {
	nStep = 0;
	TIME = 0.0;

	// grid step along the X1 axis
	dx1 = (x1_e - x1_w) / (n1_g - 1);
	// grid step along the X2 axis
	dx2 = (x2_n - x2_s) / (n2_g - 1);
	// grid step along the X3 axis
	dx3 = (x3_t - x3_b) / (n3_g - 1);

	x1[0] = x1_w + ((n1 - 1) * dx1 * coordX) - dx1;
	x2[0] = x2_s + ((n2 - 1) * dx2 * coordY) - dx2;
	x3[0] = x3_b + ((n3 - 1) * dx3 * coordZ) - dx3;

	// #####################################################
	//          block of arrays initialization
	// #####################################################

	// along the X1 axis
	for (int i = 0; i < n1 + 1 ; ++i) {
		x1[i + 1] = x1[i] + dx1;
	}

	// along the X2 axis
	for (int j = 0; j < n2 + 1; ++j) {
		x2[j + 1] = x2[j] + dx2;
	}

	// along the X3 axis
	for (int k = 0; k < n3 + 1; ++k) {
		x3[k + 1] = x3[k] + dx3;
	}

	for (int i = 0; i <= n1; ++i) {
		for (int j = 0; j <= n2; ++j) {
			for (int k = 0; k <= n3; ++k) {
				double x = 0.5 * (x1[i] + x1[i + 1]);
				double y = 0.5 * (x2[j] + x2[j + 1]);
				//u1nCon->elem(i, j, k) = (y <= 3) && (y >= 2) && (x >= 5) ? 0. : u10;
				u1nCon->elem(i, j, k) = y <= 1 ? 0. : u10;
				u2nCon->elem(i, j, k) = u20;
				u3nCon->elem(i, j, k) = u30;
				ronCon->elem(i, j, k) = ro0_g;
				tnCon->elem(i, j, k) = t0;

				averageU1->elem(i, j, k) = 0.;
				averageU2->elem(i, j, k) = 0.;
				averageU3->elem(i, j, k) = 0.;
			}
		}
	}

	for (int i = 1; i <= n1; ++i) {
		for (int j = 1; j <= n2; ++j) {
			for (int k = 1; k <= n3; ++k) {
				double y = x2[j], y_c = 0.5 * (x2[j] + x2[j + 1]);
				double x = x1[i], x_c = 0.5 * (x1[i] + x1[i + 1]);

				p1->elem(i, j, k) = p2->elem(i, j, k) = p3->elem(i, j, k) = 0.;

				// u11->elem(i, j, k) = (y_c <= 3) && (y_c >= 2) && (x_c >= 5) ? 0. : u10;
				// u12->elem(i, j, k) = (y <= 3) && (y >= 2) && (x >= 5) ? 0. : u10;
				// u13->elem(i, j, k) = (y_c <= 3) && (y_c >= 2) && (x_c >= 5) ? 0. : u10;

				u11->elem(i, j, k) = y_c <= 1 ? 0. : u10;
				u12->elem(i, j, k) = y <= 1 ? 0. : u10;
				u13->elem(i, j, k) = y_c <= 1 ? 0. : u10;

				u21->elem(i, j, k) = u22->elem(i, j, k) = u23->elem(i, j, k) = u20;
				u31->elem(i, j, k) = u32->elem(i, j, k) = u33->elem(i, j, k) = u30;
				ro1->elem(i, j, k) = ro2->elem(i, j, k) = ro3->elem(i, j, k) = ro0_g;
			}
		}
	}

	for (int i = 0; i < n1 + 1; ++i) {
		for (int j = 0; j < n2 + 1; ++j) {
			for (int k = 0; k < n3 + 1; ++k) {
				double x = 0.5 * (x1[i] + x1[i + 1]), y = 0.5 * (x2[j] + x2[j + 1]);
				// condition->elem(i, j, k) = ((x >= 5) && (x <= 6) && (y >= 2) && (y <= 3)) ? 1 : 0;
				condition->elem(i, j, k) = (x <= 5 && y <= 1) ? 1 : 0;
			}
		}
	}

	for (int i = 0; i < n1 + 1; ++i) {
		for (int k = 0; k < n3 + 1; ++k) {
			if (coordY == 0) condition->elem(i, 0, k) = 1;
			if (coordY == py - 1) condition->elem(i, n2, k) = 1;
		}
	}
}

void TimeStepSize() {
	double u1_c, u2_c, u3_c, dtu1, dtu2, dtu3, dtu, dtv1, dtv2, dtv3, dtv;

	dt = pow(10, 8);

	for (int i = 1; i < n1; ++i) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k < n3; ++k) {
				u1_c = u1Con->elem(i, j, k);
				u2_c = u2Con->elem(i, j, k);
				u3_c = u3Con->elem(i, j, k);

				dtu1 = CFL * dx1 / (sound + fabs(u1_c));
				dtu2 = CFL * dx2 / (sound + fabs(u2_c));
				dtu3 = CFL * dx3 / (sound + fabs(u3_c));

				// DTU = MIN(DTU1, DTU2, DTU3)
				dtu = min3d(dtu1, dtu2, dtu3);

				if (VIS > pow(10, -16)) {
					dtv1 = CFL * dx1 * dx1 / (2.*VIS);
					dtv2 = CFL * dx2 * dx2 / (2.*VIS);
					dtv3 = CFL * dx3 * dx3 / (2.*VIS);

					// DTV = MIN (DTV1, DTV2, DTV3)
					dtv = min3d(dtv1, dtv2, dtv3);
				} else {
					dtv = pow(10, 16);
				}

				// DT = MIN(DT, DTU, DTV)
				dt = min3d(dt, dtu, dtv);
			}
		}
	}

	double dttemp;

	MPI_Allreduce(&dt, &dttemp, 1, MPI_DOUBLE, MPI_MIN, newComm);

	if (dt > dttemp) dt = dttemp;
}

void Phase1() {
	// initialization
	for (int i = 1; i < n1; i++)
	{
		for (int j = 1; j < n2; j++)
		{
			for (int k = 1; k < n3; k++)
			{
				u1Con->elem(i, j, k) = u1nCon->elem(i, j, k);
				u2Con->elem(i, j, k) = u2nCon->elem(i, j, k);
				u3Con->elem(i, j, k) = u3nCon->elem(i, j, k);
				roCon->elem(i, j, k) = ronCon->elem(i, j, k);
				tCon->elem(i, j, k) = tnCon->elem(i, j, k);
			}
		}
	}

	// geometric characteristics of the computational cell
	double ds1, ds2, ds3, dvc;

	// velocity, density, temperature and pressure on the eastern plane
	double  u1_e, u2_e, u3_e, ro_e, t_e, p_e,
	        // velocity, density , temperature and pressure on the western plane
	        u1_w, u2_w, u3_w, ro_w, t_w, p_w,
	        // velocity, density , temperature and pressure on the northern plane
	        u1_n, u2_n, u3_n, ro_n, t_n, p_n,
	        // velocity, density , temperature and pressure on the southern plane
	        u1_s, u2_s, u3_s, ro_s, t_s, p_s,
	        // velocity, density , temperature and pressure on the top plane
	        u1_t, u2_t, u3_t, ro_t, t_t, p_t,
	        // velocity, density , temperature and pressure on the bottom plane
	        u1_b, u2_b, u3_b, ro_b, t_b, p_b,
	        // velocity, density and temperature in the cell center
	        u1_c, u2_c, u3_c, ro_c, t_c,
	        // velocity, density and temperature in the cell center on the next time step
	        u1_cn, u2_cn, u3_cn, ro_cn, t_cn;

	// plane squares
	ds1 = dx2 * dx3;
	ds2 = dx1 * dx3;
	ds3 = dx1 * dx2;

	// cell volume
	dvc = dx1 * dx2 * dx3;

	for (int i = 1; i < n1; i++)
	{
		for (int j = 1; j < n2; j++)
		{
			for (int k = 1; k < n3; k++)
			{
				if (condition->elem(i, j, k) < 0.5) {
					// #########################################################
					//  get velocity, density , temperature and pressure values
					// #########################################################

					// east plane
					u1_e = u11->elem(i + 1, j, k);
					u2_e = u21->elem(i + 1, j, k);
					u3_e = u31->elem(i + 1, j, k);
					ro_e = ro1->elem(i + 1, j, k);
					t_e = t1->elem(i + 1, j, k);
					p_e = p1->elem(i + 1, j, k);
					// west plane
					u1_w = u11->elem(i, j, k);
					u2_w = u21->elem(i, j, k);
					u3_w = u31->elem(i, j, k);
					ro_w = ro1->elem(i, j, k);
					t_w = t1->elem(i, j, k);
					p_w = p1->elem(i, j, k);

					// north plane
					u1_n = u12->elem(i, j + 1, k);
					u2_n = u22->elem(i, j + 1, k);
					u3_n = u32->elem(i, j + 1, k);
					ro_n = ro2->elem(i, j + 1, k);
					t_n = t2->elem(i, j + 1, k);
					p_n = p2->elem(i, j + 1, k);
					// south plane
					u1_s = u12->elem(i, j, k);
					u2_s = u22->elem(i, j, k);
					u3_s = u32->elem(i, j, k);
					ro_s = ro2->elem(i, j, k);
					t_s = t2->elem(i, j, k);
					p_s = p2->elem(i, j, k);

					// top plane
					u1_t = u13->elem(i, j, k + 1);
					u2_t = u23->elem(i, j, k + 1);
					u3_t = u33->elem(i, j, k + 1);
					ro_t = ro3->elem(i, j, k + 1);
					t_t = t3->elem(i, j, k + 1);
					p_t = p3->elem(i, j, k + 1);
					// bottom plane
					u1_b = u13->elem(i, j, k);
					u2_b = u23->elem(i, j, k);
					u3_b = u33->elem(i, j, k);
					ro_b = ro3->elem(i, j, k);
					t_b = t3->elem(i, j, k);
					p_b = p3->elem(i, j, k);

					// cell center
					u1_c = u1Con->elem(i, j, k);
					u2_c = u2Con->elem(i, j, k);
					u3_c = u3Con->elem(i, j, k);
					ro_c = roCon->elem(i, j, k);
					t_c = tCon->elem(i, j, k);

					// #####################################################
					//              new values evaluating
					// #####################################################

					// new density
					ro_cn = (ro_c * dvc - 0.5 * dt * (
					             (ro_e * u1_e - ro_w * u1_w) * ds1 +
					             (ro_n * u2_n - ro_s * u2_s) * ds2 +
					             (ro_t * u3_t - ro_b * u3_b) * ds3)) / dvc;

					// new conservative velocity along the X1 axis
					u1_cn = (ro_c * u1_c * dvc - 0.5 * dt * (
					             ((ro_e * u1_e * u1_e + p_e) - (ro_w * u1_w * u1_w + p_w)) * ds1 +
					             (ro_n * u1_n * u2_n - ro_s * u1_s * u2_s) * ds2 +
					             (ro_t * u1_t * u3_t - ro_b * u1_b * u3_b) * ds3)) / (ro_cn * dvc);

					// new conservative velocity along the X2 axis
					u2_cn = (ro_c * u2_c * dvc - 0.5 * dt * (
					             (ro_e * u2_e * u1_e - ro_w * u2_w * u1_w) * ds1 +
					             ((ro_n * u2_n * u2_n + p_n) - (ro_s * u2_s * u2_s + p_s)) * ds2 +
					             (ro_t * u2_t * u3_t - ro_b * u2_b * u3_b) * ds3)) / (ro_cn * dvc);

					// new conservative velocity along the X3 axis
					u3_cn = (ro_c * u3_c * dvc - 0.5 * dt * (
					             (ro_e * u3_e * u1_e - ro_w * u3_w * u1_w) * ds1 +
					             (ro_n * u3_n * u2_n - ro_s * u3_s * u2_s) * ds2 +
					             ((ro_t * u3_t * u3_t + p_t) - (ro_b * u3_b * u3_b + p_b)) * ds3)) / (ro_cn * dvc);

					// new temperature
					t_cn = (ro_c * t_c * dvc - 0.5 * dt * (
					            (ro_e * t_e * u1_e - ro_w * t_w * u1_w) * ds1 +
					            (ro_n * t_n * u2_n - ro_s * t_s * u2_s) * ds2 +
					            (ro_t * t_t * u3_t - ro_b * t_b * u3_b) * ds3)) / (dvc * ro_cn);
				} else {
					u1_cn = 0.;
					u2_cn = 0.;
					u3_cn = 0.;
					ro_cn = ro0_g;
					t_cn = t0;
				}

				// finally
				u1nCon->elem(i, j, k) = u1_cn;
				u2nCon->elem(i, j, k) = u2_cn;
				u3nCon->elem(i, j, k) = u3_cn;
				ronCon->elem(i, j, k) = ro_cn;
				tnCon->elem(i, j, k) = t_cn;
			}
		}
	}
}

void StressTensor() {
	// initialization of friction stress arrays
	for (int i = 0; i <= n1 + 1; ++i) {
		for (int j = 0; j <= n2 + 1; ++j) {
			for (int k = 0; k <= n3 + 1; ++k) {
				sigm11->elem(i, j, k) = 0.0;
				sigm21->elem(i, j, k) = 0.0;
				sigm31->elem(i, j, k) = 0.0;

				sigm12->elem(i, j, k) = 0.0;
				sigm22->elem(i, j, k) = 0.0;
				sigm32->elem(i, j, k) = 0.0;

				sigm13->elem(i, j, k) = 0.0;
				sigm23->elem(i, j, k) = 0.0;
				sigm33->elem(i, j, k) = 0.0;
			}
		}
	}


	// #####################################################
	//                  boudary conditions
	// #####################################################

	// along the X1 axis
	Arr2d* sendBuffLeftU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftU3 = new Arr2d(n2 - 1, n3 - 1);

	Arr2d* sendBuffRightU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightU3 = new Arr2d(n2 - 1, n3 - 1);

	Arr2d* recvBuffLeftU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftU3 = new Arr2d(n2 - 1, n3 - 1);

	Arr2d* recvBuffRightU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightU3 = new Arr2d(n2 - 1, n3 - 1);

	for (int j = 0; j < n2 - 1; ++j) {
		for (int k = 0; k < n3 - 1; ++k) {
			sendBuffLeftU1->elem(j, k) = u1Con->elem(1, j + 1, k + 1);
			sendBuffLeftU2->elem(j, k) = u2Con->elem(1, j + 1, k + 1);
			sendBuffLeftU3->elem(j, k) = u3Con->elem(1, j + 1, k + 1);

			sendBuffRightU1->elem(j, k) = u1Con->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightU2->elem(j, k) = u2Con->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightU3->elem(j, k) = u3Con->elem(n1 - 1, j + 1, k + 1);
		}
	}


	int count = (n2 - 1) * (n3 - 1);
	MPI_Sendrecv(sendBuffRightU1->ref(), count, MPI_DOUBLE, right, 0, recvBuffLeftU1->ref(), count, MPI_DOUBLE, left, 0, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightU2->ref(), count, MPI_DOUBLE, right, 1, recvBuffLeftU2->ref(), count, MPI_DOUBLE, left, 1, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightU3->ref(), count, MPI_DOUBLE, right, 2, recvBuffLeftU3->ref(), count, MPI_DOUBLE, left, 2, newComm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(sendBuffLeftU1->ref(), count, MPI_DOUBLE, left, 3, recvBuffRightU1->ref(), count, MPI_DOUBLE, right, 3, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftU2->ref(), count, MPI_DOUBLE, left, 4, recvBuffRightU2->ref(), count, MPI_DOUBLE, right, 4, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftU3->ref(), count, MPI_DOUBLE, left, 5, recvBuffRightU3->ref(), count, MPI_DOUBLE, right, 5, newComm, MPI_STATUS_IGNORE);

	for (int j = 0; j < n2 - 1; ++j)
	{
		for (int k = 0; k < n3 - 1; ++k)
		{
			u1Con->elem(0, j + 1, k + 1) = recvBuffLeftU1->elem(j, k);
			u2Con->elem(0, j + 1, k + 1) = recvBuffLeftU2->elem(j, k);
			u3Con->elem(0, j + 1, k + 1) = recvBuffLeftU3->elem(j, k);

			u1Con->elem(n1, j + 1, k + 1) = recvBuffRightU1->elem(j, k);
			u2Con->elem(n1, j + 1, k + 1) = recvBuffRightU2->elem(j, k);
			u3Con->elem(n1, j + 1, k + 1) = recvBuffRightU3->elem(j, k);
		}
	}

	delete sendBuffLeftU1;
	delete sendBuffLeftU2;
	delete sendBuffLeftU3;

	delete sendBuffRightU1;
	delete sendBuffRightU2;
	delete sendBuffRightU3;

	delete recvBuffLeftU1;
	delete recvBuffLeftU2;
	delete recvBuffLeftU3;

	delete recvBuffRightU1;
	delete recvBuffRightU2;
	delete recvBuffRightU3;

	// along the X2 axis
	Arr2d* sendBuffDownU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownU3 = new Arr2d(n1 - 1, n3 - 1);

	Arr2d* sendBuffUpU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpU3 = new Arr2d(n1 - 1, n3 - 1);

	Arr2d* recvBuffDownU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownU3 = new Arr2d(n1 - 1, n3 - 1);

	Arr2d* recvBuffUpU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpU3 = new Arr2d(n1 - 1, n3 - 1);

	for (int i = 0; i < n1 - 1; ++i) {
		for (int k = 0; k < n3 - 1; ++k) {
			sendBuffDownU1->elem(i, k) = u1Con->elem(i + 1, 1, k + 1);
			sendBuffDownU2->elem(i, k) = u2Con->elem(i + 1, 1, k + 1);
			sendBuffDownU3->elem(i, k) = u3Con->elem(i + 1, 1, k + 1);

			sendBuffUpU1->elem(i, k) = u1Con->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpU2->elem(i, k) = u2Con->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpU3->elem(i, k) = u3Con->elem(i + 1, n2 - 1, k + 1);
		}
	}

	count = (n1 - 1) * (n3 - 1);
	MPI_Sendrecv(sendBuffDownU1->ref(), count, MPI_DOUBLE, down, 6, recvBuffUpU1->ref(), count, MPI_DOUBLE, up, 6, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownU2->ref(), count, MPI_DOUBLE, down, 7, recvBuffUpU2->ref(), count, MPI_DOUBLE, up, 7, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownU3->ref(), count, MPI_DOUBLE, down, 8, recvBuffUpU3->ref(), count, MPI_DOUBLE, up, 8, newComm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(sendBuffUpU1->ref(), count, MPI_DOUBLE, up, 9, recvBuffDownU1->ref(), count, MPI_DOUBLE, down, 9, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpU2->ref(), count, MPI_DOUBLE, up, 10, recvBuffDownU2->ref(), count, MPI_DOUBLE, down, 10, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpU3->ref(), count, MPI_DOUBLE, up, 11, recvBuffDownU3->ref(), count, MPI_DOUBLE, down, 11, newComm, MPI_STATUS_IGNORE);

	for (int i = 0; i < n1 - 1; ++i) {
		for (int k = 0; k < n3 - 1; ++k) {
			u1Con->elem(i + 1, 0, k + 1) = recvBuffDownU1->elem(i, k);
			u2Con->elem(i + 1, 0, k + 1) = recvBuffDownU2->elem(i, k);
			u3Con->elem(i + 1, 0, k + 1) = recvBuffDownU3->elem(i, k);

			u1Con->elem(i + 1, n2, k + 1) = recvBuffUpU1->elem(i, k);
			u2Con->elem(i + 1, n2, k + 1) = recvBuffUpU2->elem(i, k);
			u3Con->elem(i + 1, n2, k + 1) = recvBuffUpU3->elem(i, k);
		}
	}

	delete sendBuffDownU1;
	delete sendBuffDownU2;
	delete sendBuffDownU3;

	delete sendBuffUpU1;
	delete sendBuffUpU2;
	delete sendBuffUpU3;

	delete recvBuffDownU1;
	delete recvBuffDownU2;
	delete recvBuffDownU3;

	delete recvBuffUpU1;
	delete recvBuffUpU2;
	delete recvBuffUpU3;

	// along the X3 axis
	Arr2d* sendBuffBottomU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomU3 = new Arr2d(n1 - 1, n2 - 1);

	Arr2d* sendBuffTopU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopU3 = new Arr2d(n1 - 1, n2 - 1);

	Arr2d* recvBuffBottomU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomU3 = new Arr2d(n1 - 1, n2 - 1);

	Arr2d* recvBuffTopU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopU3 = new Arr2d(n1 - 1, n2 - 1);

	for (int i = 0; i < n1 - 1; ++i) {
		for (int j = 0; j < n2 - 1; ++j) {
			sendBuffBottomU1->elem(i, j) = u1Con->elem(i + 1, j + 1, 1);
			sendBuffBottomU2->elem(i, j) = u2Con->elem(i + 1, j + 1, 1);
			sendBuffBottomU3->elem(i, j) = u3Con->elem(i + 1, j + 1, 1);

			sendBuffTopU1->elem(i, j) = u1Con->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopU2->elem(i, j) = u2Con->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopU3->elem(i, j) = u3Con->elem(i + 1, j + 1, n3 - 1);
		}
	}

	count = (n1 - 1) * (n2 - 1);
	MPI_Sendrecv(sendBuffBottomU1->ref(), count, MPI_DOUBLE, bottom, 12, recvBuffTopU1->ref(), count, MPI_DOUBLE, top, 12, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomU2->ref(), count, MPI_DOUBLE, bottom, 13, recvBuffTopU2->ref(), count, MPI_DOUBLE, top, 13, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomU3->ref(), count, MPI_DOUBLE, bottom, 14, recvBuffTopU3->ref(), count, MPI_DOUBLE, top, 14, newComm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(sendBuffTopU1->ref(), count, MPI_DOUBLE, top, 15, recvBuffBottomU1->ref(), count, MPI_DOUBLE, bottom, 15, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopU2->ref(), count, MPI_DOUBLE, top, 16, recvBuffBottomU2->ref(), count, MPI_DOUBLE, bottom, 16, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopU3->ref(), count, MPI_DOUBLE, top, 17, recvBuffBottomU3->ref(), count, MPI_DOUBLE, bottom, 17, newComm, MPI_STATUS_IGNORE);

	for (int i = 0; i < n1 - 1; ++i) {
		for (int j = 0; j < n2 - 1; ++j) {
			u1Con->elem(i + 1, j + 1, 0) = recvBuffBottomU1->elem(i, j);
			u2Con->elem(i + 1, j + 1, 0) = recvBuffBottomU2->elem(i, j);
			u3Con->elem(i + 1, j + 1, 0) = recvBuffBottomU3->elem(i, j);

			u1Con->elem(i + 1, j + 1, n3) = recvBuffTopU1->elem(i, j);
			u2Con->elem(i + 1, j + 1, n3) = recvBuffTopU2->elem(i, j);
			u3Con->elem(i + 1, j + 1, n3) = recvBuffTopU3->elem(i, j);
		}
	}

	delete sendBuffBottomU1;
	delete sendBuffBottomU2;
	delete sendBuffBottomU3;

	delete sendBuffTopU1;
	delete sendBuffTopU2;
	delete sendBuffTopU3;

	delete recvBuffBottomU1;
	delete recvBuffBottomU2;
	delete recvBuffBottomU3;

	delete recvBuffTopU1;
	delete recvBuffTopU2;
	delete recvBuffTopU3;

	if (coordX == px - 1) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k < n3; ++k) {
				// on the west plane
				u1Con->elem(n1, j, k) = u1Con->elem(n1 - 1, j, k);
				u2Con->elem(n1, j, k) = u2Con->elem(n1 - 1, j, k);
				u3Con->elem(n1, j, k) = u3Con->elem(n1 - 1, j, k);
			}
		}
	}

	if (coordX == 0) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k < n3; ++k) {
				// on the east plane
				u1Con->elem(0, j, k) = u1Con->elem(1, j, k);
				u2Con->elem(0, j, k) = u2Con->elem(1, j, k);
				u3Con->elem(0, j, k) = u3Con->elem(1, j, k);
			}
		}
	}

	if (coordZ == 0) {
		for (int i = 1; i < n1; ++i) {
			for (int j = 1; j < n2; ++j) {
				// on the top plane
				u1Con->elem(i, j, 0) = u1Con->elem(i, j, 1);
				u2Con->elem(i, j, 0) = u2Con->elem(i, j, 1);
				u3Con->elem(i, j, 0) = u3Con->elem(i, j, 1);
			}
		}
	}

	if (coordZ == pz - 1) {
		for (int i = 1; i < n1; ++i) {
			for (int j = 1; j < n2; ++j) {
				// on the bottom plane
				u1Con->elem(i, j, n3) = u1Con->elem(i, j, n3 - 1);
				u2Con->elem(i, j, n3) = u2Con->elem(i, j, n3 - 1);
				u3Con->elem(i, j, n3) = u3Con->elem(i, j, n3 - 1);
			}
		}
	}

	// #####################################################
	//              bypassing along the faces
	// #####################################################

	double v1, v2, v3;
	double cond_c, cond_cw;
	double u1_c, u1_cw, u2_c, u2_cw, u3_c, u3_cw;

	// bypassing along the face perpendicular to X1
	for (int k = 1; k < n3; ++k) {
		for (int j = 1; j < n2; ++j) {
			for (int i = 1; i <= n1; ++i) {
				// velocity components in cell centers
				u1_c = u1Con->elem(i, j, k);
				u1_cw = u1Con->elem(i - 1, j, k);

				u2_c = u2Con->elem(i, j, k);
				u2_cw = u2Con->elem(i - 1, j, k);

				u3_c = u3Con->elem(i, j, k);
				u3_cw = u3Con->elem(i - 1, j, k);

				cond_c = condition->elem(i, j, k);
				cond_cw = condition->elem(i - 1, j, k);

				if (cond_cw < 0.5 && cond_c < 0.5) {
					v1 = VIS * (u1_c - u1_cw) / dx1;
					v2 = VIS * (u2_c - u2_cw) / dx1;
					v3 = VIS * (u3_c - u3_cw) / dx1;
				} else if (cond_cw < 0.5 && cond_c > 0.5) {
					v1 = 2 * VIS * (-u1_cw) / dx1;
					v2 = 2 * VIS * (-u2_cw) / dx1;
					v3 = 2 * VIS * (-u3_cw) / dx1;
				} else if (cond_cw > 0.5 && cond_c < 0.5) {
					v1 = 2 * VIS * u1_c / dx1;
					v2 = 2 * VIS * u2_c / dx1;
					v3 = 2 * VIS * u3_c / dx1;
				} else {
					v1 = v2 = v3 = 0.;
				}

				// friction stress
				sigm11->elem(i, j, k) = v1;
				sigm21->elem(i, j, k) = v2;
				sigm31->elem(i, j, k) = v3;
			}
		}
	}

	double cond_cs;
	double u1_cs, u2_cs, u3_cs;

	// bypassing along the face perpenditcular to X2
	for (int k = 1; k < n3; ++k) {
		for (int i = 1; i < n1; ++i) {
			for (int j = 1; j <= n2; ++j) {
				// velocity components in cell centers
				u1_c = u1Con->elem(i, j, k);
				u1_cs = u1Con->elem(i, j - 1, k);

				u2_c = u2Con->elem(i, j, k);
				u2_cs = u2Con->elem(i, j - 1, k);

				u3_c = u3Con->elem(i, j, k);
				u3_cs = u3Con->elem(i, j - 1, k);

				cond_c = condition->elem(i, j, k);
				cond_cs = condition->elem(i, j - 1, k);

				if (cond_cs < 0.5 && cond_c < 0.5) {
					v1 = VIS * (u1_c - u1_cs) / dx2;
					v2 = VIS * (u2_c - u2_cs) / dx2;
					v3 = VIS * (u3_c - u3_cs) / dx2;
				} else if (cond_cs < 0.5 && cond_c > 0.5) {
					v1 = 2 * VIS * (-u1_cs) / dx2;
					v2 = 2 * VIS * (-u2_cs) / dx2;
					v3 = 2 * VIS * (-u3_cs) / dx2;
				} else if (cond_cs > 0.5 && cond_c < 0.5) {
					v1 = 2 * VIS * u1_c / dx2;
					v2 = 2 * VIS * u2_c / dx2;
					v3 = 2 * VIS * u3_c / dx2;
				} else {
					v1 = v2 = v3 = 0.;
				}

				// friction stress
				sigm12->elem(i, j, k) = v1;
				sigm22->elem(i, j, k) = v2;
				sigm32->elem(i, j, k) = v3;
			}
		}
	}

	double cond_cb;
	double u1_cb, u2_cb, u3_cb;

	// bypassing along the face perpenditcular to X3
	for (int i = 1; i < n1; ++i) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k <= n3; ++k) {
				// velocity components in the cell centers
				u1_c = u1Con->elem(i, j, k);
				u1_cb = u1Con->elem(i, j, k - 1);

				u2_c = u2Con->elem(i, j, k);
				u2_cb = u2Con->elem(i, j, k - 1);

				u3_c = u3Con->elem(i, j, k);
				u3_cb = u3Con->elem(i, j, k - 1);

				cond_c = condition->elem(i, j, k);
				cond_cb = condition->elem(i, j, k - 1);

				if (cond_cb < 0.5 && cond_c < 0.5) {
					v1 = VIS * (u1_c - u1_cb) / dx3;
					v2 = VIS * (u2_c - u2_cb) / dx3;
					v3 = VIS * (u3_c - u3_cb) / dx3;
				} else if (cond_cb < 0.5 && cond_c > 0.5) {
					v1 = 2 * VIS * (-u1_cb) / dx3;
					v2 = 2 * VIS * (-u2_cb) / dx3;
					v3 = 2 * VIS * (-u3_cb) / dx3;
				} else if (cond_cb > 0.5 && cond_c < 0.5) {
					v1 = 2 * VIS * u1_c / dx3;
					v2 = 2 * VIS * u2_c / dx3;
					v3 = 2 * VIS * u3_c / dx3;
				} else {
					v1 = v2 = v3 = 0.;
				}

				// friction stress
				sigm13->elem(i, j, k) = v1;
				sigm23->elem(i, j, k) = v2;
				sigm33->elem(i, j, k) = v3;
			}
		}
	}

	// #####################################################
	//              friction forces computation
	// #####################################################

	double ds1, ds2, ds3;

	// area of the face perpendicuar to x1
	ds1 = dx2 * dx3;
	ds2 = dx1 * dx3;
	ds3 = dx1 * dx2;

	for (int i = 1; i < n1; ++i) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k < n3; ++k) {
				// friction forces
				if (condition->elem(i, j, k) == 0) {
					double ro_c = roCon->elem(i, j, k);

					f1->elem(i, j, k) =
					    ((sigm11->elem(i + 1, j, k) - sigm11->elem(i, j, k)) * ds1 +
					     (sigm12->elem(i, j + 1, k) - sigm12->elem(i, j, k)) * ds2 +
					     (sigm13->elem(i, j, k + 1) - sigm13->elem(i, j, k)) * ds3) * ro_c;

					f2->elem(i, j, k) =
					    ((sigm21->elem(i + 1, j, k) - sigm21->elem(i, j, k)) * ds1 +
					     (sigm22->elem(i, j + 1, k) - sigm22->elem(i, j, k)) * ds2 +
					     (sigm23->elem(i, j, k + 1) - sigm23->elem(i, j, k)) * ds3) * ro_c;

					f3->elem(i, j, k) =
					    ((sigm31->elem(i + 1, j, k) - sigm31->elem(i, j, k)) * ds1 +
					     (sigm32->elem(i, j + 1, k) - sigm32->elem(i, j, k)) * ds2 +
					     (sigm33->elem(i, j, k + 1) - sigm33->elem(i, j, k)) * ds3) * ro_c;
				} else {
					f1->elem(i, j, k) = 0.;
					f2->elem(i, j, k) = 0.;
					f3->elem(i, j, k) = 0.;
				}
			}
		}
	}
}

void UseForces() {
	double dvc, ro_cn;

	// cell volume
	dvc = dx1 * dx2 * dx3;

	for (int i = 1; i < n1; ++i) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k < n3; ++k) {
				if (condition->elem(i, j, k) == 0) {
					ro_cn = ronCon->elem(i, j, k);

					u1nCon->elem(i, j, k) = (ro_cn * dvc * u1nCon->elem(i, j, k) + 0.5 * dt * f1->elem(i, j, k)) / (dvc * ro_cn);
					u2nCon->elem(i, j, k) = (ro_cn * dvc * u2nCon->elem(i, j, k) + 0.5 * dt * f2->elem(i, j, k)) / (dvc * ro_cn);
					u3nCon->elem(i, j, k) = (ro_cn * dvc * u3nCon->elem(i, j, k) + 0.5 * dt * f3->elem(i, j, k)) / (dvc * ro_cn);
				}
			}
		}
	}
}

void Phase2() {
	// along the X1 axis
	Arr2d* sendBuffLeftU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftU3 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftRo = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftRon = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftT = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftTn = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftU11 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftU21 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftU31 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftRo1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffLeftP1 = new Arr2d(n2 - 1, n3 - 1);

	Arr2d* sendBuffRightU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightU3 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightRo = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightRon = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightT = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightTn = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightU11 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightU21 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightU31 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightRo1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* sendBuffRightP1 = new Arr2d(n2 - 1, n3 - 1);

	Arr2d* recvBuffLeftU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftU3 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftRo = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftRon = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftT = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftTn = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftU11 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftU21 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftU31 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftRo1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffLeftP1 = new Arr2d(n2 - 1, n3 - 1);

	Arr2d* recvBuffRightU1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightU2 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightU3 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightRo = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightRon = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightT = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightTn = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightU11 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightU21 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightU31 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightRo1 = new Arr2d(n2 - 1, n3 - 1);
	Arr2d* recvBuffRightP1 = new Arr2d(n2 - 1, n3 - 1);

	for (int j = 0; j < n2 - 1; ++j) {
		for (int k = 0; k < n3 - 1; ++k) {
			sendBuffLeftU1->elem(j, k) = u1nCon->elem(1, j + 1, k + 1);
			sendBuffLeftU2->elem(j, k) = u2nCon->elem(1, j + 1, k + 1);
			sendBuffLeftU3->elem(j, k) = u3nCon->elem(1, j + 1, k + 1);
			sendBuffLeftRo->elem(j, k) = roCon->elem(1, j + 1, k + 1);
			sendBuffLeftRon->elem(j, k) = ronCon->elem(1, j + 1, k + 1);
			sendBuffLeftT->elem(j, k) = tCon->elem(1, j + 1, k + 1);
			sendBuffLeftTn->elem(j, k) = tnCon->elem(1, j + 1, k + 1);
			sendBuffLeftU11->elem(j, k) = u11->elem(2, j + 1, k + 1);
			sendBuffLeftU21->elem(j, k) = u21->elem(2, j + 1, k + 1);
			sendBuffLeftU31->elem(j, k) = u31->elem(2, j + 1, k + 1);
			sendBuffLeftRo1->elem(j, k) = ro1->elem(2, j + 1, k + 1);
			sendBuffLeftP1->elem(j, k) = p1->elem(2, j + 1, k + 1);

			sendBuffRightU1->elem(j, k) = u1nCon->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightU2->elem(j, k) = u2nCon->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightU3->elem(j, k) = u3nCon->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightRo->elem(j, k) = roCon->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightRon->elem(j, k) = ronCon->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightT->elem(j, k) = tCon->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightTn->elem(j, k) = tnCon->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightU11->elem(j, k) = u11->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightU21->elem(j, k) = u21->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightU31->elem(j, k) = u31->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightRo1->elem(j, k) = ro1->elem(n1 - 1, j + 1, k + 1);
			sendBuffRightP1->elem(j, k) = p1->elem(n1 - 1, j + 1, k + 1);
		}
	}


	int count = (n2 - 1) * (n3 - 1);
	MPI_Sendrecv(sendBuffRightU1->ref(), count, MPI_DOUBLE, right, 18, recvBuffLeftU1->ref(), count, MPI_DOUBLE, left, 18, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightU2->ref(), count, MPI_DOUBLE, right, 19, recvBuffLeftU2->ref(), count, MPI_DOUBLE, left, 19, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightU3->ref(), count, MPI_DOUBLE, right, 20, recvBuffLeftU3->ref(), count, MPI_DOUBLE, left, 20, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightRo->ref(), count, MPI_DOUBLE, right, 21, recvBuffLeftRo->ref(), count, MPI_DOUBLE, left, 21, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightRon->ref(), count, MPI_DOUBLE, right, 100, recvBuffLeftRon->ref(), count, MPI_DOUBLE, left, 100, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightT->ref(), count, MPI_DOUBLE, right, 22, recvBuffLeftT->ref(), count, MPI_DOUBLE, left, 22, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightTn->ref(), count, MPI_DOUBLE, right, 101, recvBuffLeftTn->ref(), count, MPI_DOUBLE, left, 101, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightU11->ref(), count, MPI_DOUBLE, right, 23, recvBuffLeftU11->ref(), count, MPI_DOUBLE, left, 23, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightU21->ref(), count, MPI_DOUBLE, right, 24, recvBuffLeftU21->ref(), count, MPI_DOUBLE, left, 24, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightU31->ref(), count, MPI_DOUBLE, right, 25, recvBuffLeftU31->ref(), count, MPI_DOUBLE, left, 25, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightRo1->ref(), count, MPI_DOUBLE, right, 26, recvBuffLeftRo1->ref(), count, MPI_DOUBLE, left, 26, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffRightP1->ref(), count, MPI_DOUBLE, right, 27, recvBuffLeftP1->ref(), count, MPI_DOUBLE, left, 27, newComm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(sendBuffLeftU1->ref(), count, MPI_DOUBLE, left, 28, recvBuffRightU1->ref(), count, MPI_DOUBLE, right, 28, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftU2->ref(), count, MPI_DOUBLE, left, 29, recvBuffRightU2->ref(), count, MPI_DOUBLE, right, 29, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftU3->ref(), count, MPI_DOUBLE, left, 30, recvBuffRightU3->ref(), count, MPI_DOUBLE, right, 30, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftRo->ref(), count, MPI_DOUBLE, left, 31, recvBuffRightRo->ref(), count, MPI_DOUBLE, right, 31, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftRon->ref(), count, MPI_DOUBLE, left, 102, recvBuffRightRon->ref(), count, MPI_DOUBLE, right, 102, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftT->ref(), count, MPI_DOUBLE, left, 32, recvBuffRightT->ref(), count, MPI_DOUBLE, right, 32, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftTn->ref(), count, MPI_DOUBLE, left, 103, recvBuffRightTn->ref(), count, MPI_DOUBLE, right, 103, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftU11->ref(), count, MPI_DOUBLE, left, 33, recvBuffRightU11->ref(), count, MPI_DOUBLE, right, 33, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftU21->ref(), count, MPI_DOUBLE, left, 34, recvBuffRightU21->ref(), count, MPI_DOUBLE, right, 34, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftU31->ref(), count, MPI_DOUBLE, left, 35, recvBuffRightU31->ref(), count, MPI_DOUBLE, right, 35, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftRo1->ref(), count, MPI_DOUBLE, left, 36, recvBuffRightRo1->ref(), count, MPI_DOUBLE, right, 36, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffLeftP1->ref(), count, MPI_DOUBLE, left, 37, recvBuffRightP1->ref(), count, MPI_DOUBLE, right, 37, newComm, MPI_STATUS_IGNORE);

	for (int j = 0; j < n2 - 1; ++j)
	{
		for (int k = 0; k < n3 - 1; ++k)
		{
			u1nCon->elem(0, j + 1, k + 1) = recvBuffLeftU1->elem(j, k);
			u2nCon->elem(0, j + 1, k + 1) = recvBuffLeftU2->elem(j, k);
			u3nCon->elem(0, j + 1, k + 1) = recvBuffLeftU3->elem(j, k);
			roCon->elem(0, j + 1, k + 1) = recvBuffLeftRo->elem(j, k);
			ronCon->elem(0, j + 1, k + 1) = recvBuffLeftRon->elem(j, k);
			tCon->elem(0, j + 1, k + 1) = recvBuffLeftT->elem(j, k);
			tnCon->elem(0, j + 1, k + 1) = recvBuffLeftTn->elem(j, k);
			u11->elem(0, j + 1, k + 1) = recvBuffLeftU11->elem(j, k);
			u21->elem(0, j + 1, k + 1) = recvBuffLeftU21->elem(j, k);
			u31->elem(0, j + 1, k + 1) = recvBuffLeftU31->elem(j, k);
			ro1->elem(0, j + 1, k + 1) = recvBuffLeftRo1->elem(j, k);
			p1->elem(0, j + 1, k + 1) = recvBuffLeftP1->elem(j, k);

			u1nCon->elem(n1, j + 1, k + 1) = recvBuffRightU1->elem(j, k);
			u2nCon->elem(n1, j + 1, k + 1) = recvBuffRightU2->elem(j, k);
			u3nCon->elem(n1, j + 1, k + 1) = recvBuffRightU3->elem(j, k);
			roCon->elem(n1, j + 1, k + 1) = recvBuffRightRo->elem(j, k);
			ronCon->elem(n1, j + 1, k + 1) = recvBuffRightRon->elem(j, k);
			tCon->elem(n1, j + 1, k + 1) = recvBuffRightT->elem(j, k);
			tnCon->elem(n1, j + 1, k + 1) = recvBuffRightTn->elem(j, k);
			u11->elem(n1 + 1, j + 1, k + 1) = recvBuffRightU11->elem(j, k);
			u21->elem(n1 + 1, j + 1, k + 1) = recvBuffRightU21->elem(j, k);
			u31->elem(n1 + 1, j + 1, k + 1) = recvBuffRightU31->elem(j, k);
			ro1->elem(n1 + 1, j + 1, k + 1) = recvBuffRightRo1->elem(j, k);
			p1->elem(n1 + 1, j + 1, k + 1) = recvBuffRightP1->elem(j, k);
		}
	}

	delete sendBuffLeftU1;
	delete sendBuffLeftU2;
	delete sendBuffLeftU3;
	delete sendBuffLeftRo;
	delete sendBuffLeftRon;
	delete sendBuffLeftT;
	delete sendBuffLeftTn;
	delete sendBuffLeftU11;
	delete sendBuffLeftU21;
	delete sendBuffLeftU31;
	delete sendBuffLeftRo1;
	delete sendBuffLeftP1;

	delete sendBuffRightU1;
	delete sendBuffRightU2;
	delete sendBuffRightU3;
	delete sendBuffRightRo;
	delete sendBuffRightRon;
	delete sendBuffRightT;
	delete sendBuffRightTn;
	delete sendBuffRightU11;
	delete sendBuffRightU21;
	delete sendBuffRightU31;
	delete sendBuffRightRo1;
	delete sendBuffRightP1;

	delete recvBuffLeftU1;
	delete recvBuffLeftU2;
	delete recvBuffLeftU3;
	delete recvBuffLeftRo;
	delete recvBuffLeftRon;
	delete recvBuffLeftT;
	delete recvBuffLeftTn;
	delete recvBuffLeftU11;
	delete recvBuffLeftU21;
	delete recvBuffLeftU31;
	delete recvBuffLeftRo1;
	delete recvBuffLeftP1;

	delete recvBuffRightU1;
	delete recvBuffRightU2;
	delete recvBuffRightU3;
	delete recvBuffRightRo;
	delete recvBuffRightRon;
	delete recvBuffRightT;
	delete recvBuffRightTn;
	delete recvBuffRightU11;
	delete recvBuffRightU21;
	delete recvBuffRightU31;
	delete recvBuffRightRo1;
	delete recvBuffRightP1;

	// along the X2 axis
	Arr2d* sendBuffDownU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownU3 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownRo = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownRon = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownT = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownTn = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownU12 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownU22 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownU32 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownRo2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffDownP2 = new Arr2d(n1 - 1, n3 - 1);

	Arr2d* sendBuffUpU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpU3 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpRo = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpRon = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpT = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpTn = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpU12 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpU22 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpU32 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpRo2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* sendBuffUpP2 = new Arr2d(n1 - 1, n3 - 1);

	Arr2d* recvBuffDownU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownU3 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownRo = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownRon = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownT = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownTn = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownU12 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownU22 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownU32 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownRo2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffDownP2 = new Arr2d(n1 - 1, n3 - 1);

	Arr2d* recvBuffUpU1 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpU2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpU3 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpRo = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpRon = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpT = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpTn = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpU12 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpU22 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpU32 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpRo2 = new Arr2d(n1 - 1, n3 - 1);
	Arr2d* recvBuffUpP2 = new Arr2d(n1 - 1, n3 - 1);

	for (int i = 0; i < n1 - 1; ++i) {
		for (int k = 0; k < n3 - 1; ++k) {
			sendBuffDownU1->elem(i, k) = u1nCon->elem(i + 1, 1, k + 1);
			sendBuffDownU2->elem(i, k) = u2nCon->elem(i + 1, 1, k + 1);
			sendBuffDownU3->elem(i, k) = u3nCon->elem(i + 1, 1, k + 1);
			sendBuffDownRo->elem(i, k) = roCon->elem(i + 1, 1, k + 1);
			sendBuffDownRon->elem(i, k) = ronCon->elem(i + 1, 1, k + 1);
			sendBuffDownT->elem(i, k) = tCon->elem(i + 1, 1, k + 1);
			sendBuffDownTn->elem(i, k) = tnCon->elem(i + 1, 1, k + 1);
			sendBuffDownU12->elem(i, k) = u12->elem(i + 1, 2, k + 1);
			sendBuffDownU22->elem(i, k) = u22->elem(i + 1, 2, k + 1);
			sendBuffDownU32->elem(i, k) = u32->elem(i + 1, 2, k + 1);
			sendBuffDownRo2->elem(i, k) = ro2->elem(i + 1, 2, k + 1);
			sendBuffDownP2->elem(i, k) = p2->elem(i + 1, 2, k + 1);

			sendBuffUpU1->elem(i, k) = u1nCon->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpU2->elem(i, k) = u2nCon->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpU3->elem(i, k) = u3nCon->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpRo->elem(i, k) = roCon->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpRon->elem(i, k) = ronCon->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpT->elem(i, k) = tCon->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpTn->elem(i, k) = tnCon->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpU12->elem(i, k) = u12->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpU22->elem(i, k) = u22->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpU32->elem(i, k) = u32->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpRo2->elem(i, k) = ro2->elem(i + 1, n2 - 1, k + 1);
			sendBuffUpP2->elem(i, k) = p2->elem(i + 1, n2 - 1, k + 1);
		}
	}


	count = (n1 - 1) * (n3 - 1);
	MPI_Sendrecv(sendBuffUpU1->ref(), count, MPI_DOUBLE, up, 38, recvBuffDownU1->ref(), count, MPI_DOUBLE, down, 38, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpU2->ref(), count, MPI_DOUBLE, up, 39, recvBuffDownU2->ref(), count, MPI_DOUBLE, down, 39, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpU3->ref(), count, MPI_DOUBLE, up, 40, recvBuffDownU3->ref(), count, MPI_DOUBLE, down, 40, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpRo->ref(), count, MPI_DOUBLE, up, 41, recvBuffDownRo->ref(), count, MPI_DOUBLE, down, 41, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpRon->ref(), count, MPI_DOUBLE, up, 104, recvBuffDownRon->ref(), count, MPI_DOUBLE, down, 104, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpT->ref(), count, MPI_DOUBLE, up, 42, recvBuffDownT->ref(), count, MPI_DOUBLE, down, 42, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpTn->ref(), count, MPI_DOUBLE, up, 105, recvBuffDownTn->ref(), count, MPI_DOUBLE, down, 105, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpU12->ref(), count, MPI_DOUBLE, up, 43, recvBuffDownU12->ref(), count, MPI_DOUBLE, down, 43, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpU22->ref(), count, MPI_DOUBLE, up, 44, recvBuffDownU22->ref(), count, MPI_DOUBLE, down, 44, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpU32->ref(), count, MPI_DOUBLE, up, 45, recvBuffDownU32->ref(), count, MPI_DOUBLE, down, 45, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpRo2->ref(), count, MPI_DOUBLE, up, 46, recvBuffDownRo2->ref(), count, MPI_DOUBLE, down, 46, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffUpP2->ref(), count, MPI_DOUBLE, up, 47, recvBuffDownP2->ref(), count, MPI_DOUBLE, down, 47, newComm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(sendBuffDownU1->ref(), count, MPI_DOUBLE, down, 48, recvBuffUpU1->ref(), count, MPI_DOUBLE, up, 48, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownU2->ref(), count, MPI_DOUBLE, down, 49, recvBuffUpU2->ref(), count, MPI_DOUBLE, up, 49, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownU3->ref(), count, MPI_DOUBLE, down, 50, recvBuffUpU3->ref(), count, MPI_DOUBLE, up, 50, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownRo->ref(), count, MPI_DOUBLE, down, 51, recvBuffUpRo->ref(), count, MPI_DOUBLE, up, 51, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownRon->ref(), count, MPI_DOUBLE, down, 106, recvBuffUpRon->ref(), count, MPI_DOUBLE, up, 106, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownT->ref(), count, MPI_DOUBLE, down, 52, recvBuffUpT->ref(), count, MPI_DOUBLE, up, 52, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownTn->ref(), count, MPI_DOUBLE, down, 107, recvBuffUpTn->ref(), count, MPI_DOUBLE, up, 107, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownU12->ref(), count, MPI_DOUBLE, down, 53, recvBuffUpU12->ref(), count, MPI_DOUBLE, up, 53, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownU22->ref(), count, MPI_DOUBLE, down, 54, recvBuffUpU22->ref(), count, MPI_DOUBLE, up, 54, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownU32->ref(), count, MPI_DOUBLE, down, 55, recvBuffUpU32->ref(), count, MPI_DOUBLE, up, 55, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownRo2->ref(), count, MPI_DOUBLE, down, 56, recvBuffUpRo2->ref(), count, MPI_DOUBLE, up, 56, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffDownP2->ref(), count, MPI_DOUBLE, down, 57, recvBuffUpP2->ref(), count, MPI_DOUBLE, up, 57, newComm, MPI_STATUS_IGNORE);

	for (int i = 0; i < n1 - 1; ++i)
	{
		for (int k = 0; k < n3 - 1; ++k)
		{
			u1nCon->elem(i + 1, 0, k + 1) = recvBuffDownU1->elem(i, k);
			u2nCon->elem(i + 1, 0, k + 1) = recvBuffDownU2->elem(i, k);
			u3nCon->elem(i + 1, 0, k + 1) = recvBuffDownU3->elem(i, k);
			roCon->elem(i + 1, 0, k + 1) = recvBuffDownRo->elem(i, k);
			ronCon->elem(i + 1, 0, k + 1) = recvBuffDownRon->elem(i, k);
			tCon->elem(i + 1, 0, k + 1) = recvBuffDownT->elem(i, k);
			tnCon->elem(i + 1, 0, k + 1) = recvBuffDownTn->elem(i, k);
			u12->elem(i + 1, 0, k + 1) = recvBuffDownU12->elem(i, k);
			u22->elem(i + 1, 0, k + 1) = recvBuffDownU22->elem(i, k);
			u32->elem(i + 1, 0, k + 1) = recvBuffDownU32->elem(i, k);
			ro2->elem(i + 1, 0, k + 1) = recvBuffDownRo2->elem(i, k);
			p2->elem(i + 1, 0, k + 1) = recvBuffDownP2->elem(i, k);

			u1nCon->elem(i + 1, n2, k + 1) = recvBuffUpU1->elem(i, k);
			u2nCon->elem(i + 1, n2, k + 1) = recvBuffUpU2->elem(i, k);
			u3nCon->elem(i + 1, n2, k + 1) = recvBuffUpU3->elem(i, k);
			roCon->elem(i + 1, n2, k + 1) = recvBuffUpRo->elem(i, k);
			ronCon->elem(i + 1, n2, k + 1) = recvBuffUpRon->elem(i, k);
			tCon->elem(i + 1, n2, k + 1) = recvBuffUpT->elem(i, k);
			tnCon->elem(i + 1, n2, k + 1) = recvBuffUpTn->elem(i, k);
			u12->elem(i + 1, n2 + 1, k + 1) = recvBuffUpU12->elem(i, k);
			u22->elem(i + 1, n2 + 1, k + 1) = recvBuffUpU22->elem(i, k);
			u32->elem(i + 1, n2 + 1, k + 1) = recvBuffUpU32->elem(i, k);
			ro2->elem(i + 1, n2 + 1, k + 1) = recvBuffUpRo2->elem(i, k);
			p2->elem(i + 1, n2 + 1, k + 1) = recvBuffUpP2->elem(i, k);
		}
	}

	delete sendBuffDownU1;
	delete sendBuffDownU2;
	delete sendBuffDownU3;
	delete sendBuffDownRo;
	delete sendBuffDownRon;
	delete sendBuffDownT;
	delete sendBuffDownTn;
	delete sendBuffDownU12;
	delete sendBuffDownU22;
	delete sendBuffDownU32;
	delete sendBuffDownRo2;
	delete sendBuffDownP2;

	delete sendBuffUpU1;
	delete sendBuffUpU2;
	delete sendBuffUpU3;
	delete sendBuffUpRo;
	delete sendBuffUpRon;
	delete sendBuffUpT;
	delete sendBuffUpTn;
	delete sendBuffUpU12;
	delete sendBuffUpU22;
	delete sendBuffUpU32;
	delete sendBuffUpRo2;
	delete sendBuffUpP2;

	delete recvBuffDownU1;
	delete recvBuffDownU2;
	delete recvBuffDownU3;
	delete recvBuffDownRo;
	delete recvBuffDownRon;
	delete recvBuffDownT;
	delete recvBuffDownTn;
	delete recvBuffDownU12;
	delete recvBuffDownU22;
	delete recvBuffDownU32;
	delete recvBuffDownRo2;
	delete recvBuffDownP2;

	delete recvBuffUpU1;
	delete recvBuffUpU2;
	delete recvBuffUpU3;
	delete recvBuffUpRo;
	delete recvBuffUpRon;
	delete recvBuffUpT;
	delete recvBuffUpTn;
	delete recvBuffUpU12;
	delete recvBuffUpU22;
	delete recvBuffUpU32;
	delete recvBuffUpRo2;
	delete recvBuffUpP2;

	// along the X3 axis
	Arr2d* sendBuffBottomU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomU3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomRo = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomRon = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomT = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomTn = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomU13 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomU23 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomU33 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomRo3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffBottomP3 = new Arr2d(n1 - 1, n2 - 1);

	Arr2d* sendBuffTopU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopU3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopRo = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopRon = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopT = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopTn = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopU13 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopU23 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopU33 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopRo3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* sendBuffTopP3 = new Arr2d(n1 - 1, n2 - 1);

	Arr2d* recvBuffBottomU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomU3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomRo = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomRon = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomT = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomTn = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomU13 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomU23 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomU33 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomRo3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffBottomP3 = new Arr2d(n1 - 1, n2 - 1);

	Arr2d* recvBuffTopU1 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopU2 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopU3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopRo = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopRon = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopT = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopTn = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopU13 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopU23 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopU33 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopRo3 = new Arr2d(n1 - 1, n2 - 1);
	Arr2d* recvBuffTopP3 = new Arr2d(n1 - 1, n2 - 1);

	for (int i = 0; i < n1 - 1; ++i) {
		for (int j = 0; j < n2 - 1; ++j) {
			sendBuffBottomU1->elem(i, j) = u1nCon->elem(i + 1, j + 1, 1);
			sendBuffBottomU2->elem(i, j) = u2nCon->elem(i + 1, j + 1, 1);
			sendBuffBottomU3->elem(i, j) = u3nCon->elem(i + 1, j + 1, 1);
			sendBuffBottomRo->elem(i, j) = roCon->elem(i + 1, j + 1, 1);
			sendBuffBottomRon->elem(i, j) = ronCon->elem(i + 1, j + 1, 1);
			sendBuffBottomT->elem(i, j) = tCon->elem(i + 1, j + 1, 1);
			sendBuffBottomTn->elem(i, j) = tnCon->elem(i + 1, j + 1, 1);
			sendBuffBottomU13->elem(i, j) = u13->elem(i + 1, j + 1, 2);
			sendBuffBottomU23->elem(i, j) = u23->elem(i + 1, j + 1, 2);
			sendBuffBottomU33->elem(i, j) = u33->elem(i + 1, j + 1, 2);
			sendBuffBottomRo3->elem(i, j) = ro3->elem(i + 1, j + 1, 2);
			sendBuffBottomP3->elem(i, j) = p3->elem(i + 1, j + 1, 2);

			sendBuffTopU1->elem(i, j) = u1nCon->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopU2->elem(i, j) = u2nCon->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopU3->elem(i, j) = u3nCon->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopRo->elem(i, j) = roCon->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopRon->elem(i, j) = ronCon->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopT->elem(i, j) = tCon->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopTn->elem(i, j) = tnCon->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopU13->elem(i, j) = u13->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopU23->elem(i, j) = u23->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopU33->elem(i, j) = u33->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopRo3->elem(i, j) = ro3->elem(i + 1, j + 1, n3 - 1);
			sendBuffTopP3->elem(i, j) = p3->elem(i + 1, j + 1, n3 - 1);
		}
	}


	count = (n1 - 1) * (n2 - 1);
	MPI_Sendrecv(sendBuffTopU1->ref(), count, MPI_DOUBLE, top, 57, recvBuffBottomU1->ref(), count, MPI_DOUBLE, bottom, 57, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopU2->ref(), count, MPI_DOUBLE, top, 58, recvBuffBottomU2->ref(), count, MPI_DOUBLE, bottom, 58, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopU3->ref(), count, MPI_DOUBLE, top, 59, recvBuffBottomU3->ref(), count, MPI_DOUBLE, bottom, 59, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopRo->ref(), count, MPI_DOUBLE, top, 60, recvBuffBottomRo->ref(), count, MPI_DOUBLE, bottom, 60, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopRon->ref(), count, MPI_DOUBLE, top, 108, recvBuffBottomRon->ref(), count, MPI_DOUBLE, bottom, 108, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopT->ref(), count, MPI_DOUBLE, top, 61, recvBuffBottomT->ref(), count, MPI_DOUBLE, bottom, 61, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopTn->ref(), count, MPI_DOUBLE, top, 109, recvBuffBottomTn->ref(), count, MPI_DOUBLE, bottom, 109, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopU13->ref(), count, MPI_DOUBLE, top, 62, recvBuffBottomU13->ref(), count, MPI_DOUBLE, bottom, 62, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopU23->ref(), count, MPI_DOUBLE, top, 63, recvBuffBottomU23->ref(), count, MPI_DOUBLE, bottom, 63, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopU33->ref(), count, MPI_DOUBLE, top, 64, recvBuffBottomU33->ref(), count, MPI_DOUBLE, bottom, 64, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopRo3->ref(), count, MPI_DOUBLE, top, 65, recvBuffBottomRo3->ref(), count, MPI_DOUBLE, bottom, 65, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffTopP3->ref(), count, MPI_DOUBLE, top, 66, recvBuffBottomP3->ref(), count, MPI_DOUBLE, bottom, 66, newComm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(sendBuffBottomU1->ref(), count, MPI_DOUBLE, bottom, 67, recvBuffTopU1->ref(), count, MPI_DOUBLE, top, 67, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomU2->ref(), count, MPI_DOUBLE, bottom, 68, recvBuffTopU2->ref(), count, MPI_DOUBLE, top, 68, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomU3->ref(), count, MPI_DOUBLE, bottom, 69, recvBuffTopU3->ref(), count, MPI_DOUBLE, top, 69, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomRo->ref(), count, MPI_DOUBLE, bottom, 70, recvBuffTopRo->ref(), count, MPI_DOUBLE, top, 70, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomRon->ref(), count, MPI_DOUBLE, bottom, 110, recvBuffTopRon->ref(), count, MPI_DOUBLE, top, 110, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomT->ref(), count, MPI_DOUBLE, bottom, 71, recvBuffTopT->ref(), count, MPI_DOUBLE, top, 71, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomTn->ref(), count, MPI_DOUBLE, bottom, 111, recvBuffTopTn->ref(), count, MPI_DOUBLE, top, 111, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomU13->ref(), count, MPI_DOUBLE, bottom, 72, recvBuffTopU13->ref(), count, MPI_DOUBLE, top, 72, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomU23->ref(), count, MPI_DOUBLE, bottom, 73, recvBuffTopU23->ref(), count, MPI_DOUBLE, top, 73, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomU33->ref(), count, MPI_DOUBLE, bottom, 74, recvBuffTopU33->ref(), count, MPI_DOUBLE, top, 74, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomRo3->ref(), count, MPI_DOUBLE, bottom, 75, recvBuffTopRo3->ref(), count, MPI_DOUBLE, top, 75, newComm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuffBottomP3->ref(), count, MPI_DOUBLE, bottom, 76, recvBuffTopP3->ref(), count, MPI_DOUBLE, top, 76, newComm, MPI_STATUS_IGNORE);

	for (int i = 0; i < n1 - 1; ++i)
	{
		for (int j = 0; j < n2 - 1; ++j)
		{
			u1nCon->elem(i + 1, j + 1, 0) = recvBuffBottomU1->elem(i, j);
			u2nCon->elem(i + 1, j + 1, 0) = recvBuffBottomU2->elem(i, j);
			u3nCon->elem(i + 1, j + 1, 0) = recvBuffBottomU3->elem(i, j);
			roCon->elem(i + 1, j + 1, 0) = recvBuffBottomRo->elem(i, j);
			ronCon->elem(i + 1, j + 1, 0) = recvBuffBottomRon->elem(i, j);
			tCon->elem(i + 1, j + 1, 0) = recvBuffBottomT->elem(i, j);
			tnCon->elem(i + 1, j + 1, 0) = recvBuffBottomTn->elem(i, j);
			u13->elem(i + 1, j + 1, 0) = recvBuffBottomU13->elem(i, j);
			u23->elem(i + 1, j + 1, 0) = recvBuffBottomU23->elem(i, j);
			u33->elem(i + 1, j + 1, 0) = recvBuffBottomU33->elem(i, j);
			ro3->elem(i + 1, j + 1, 0) = recvBuffBottomRo3->elem(i, j);
			p3->elem(i + 1, j + 1, 0) = recvBuffBottomP3->elem(i, j);

			u1nCon->elem(i + 1, j + 1, n3) = recvBuffTopU1->elem(i, j);
			u2nCon->elem(i + 1, j + 1, n3) = recvBuffTopU2->elem(i, j);
			u3nCon->elem(i + 1, j + 1, n3) = recvBuffTopU3->elem(i, j);
			roCon->elem(i + 1, j + 1, n3) = recvBuffTopRo->elem(i, j);
			ronCon->elem(i + 1, j + 1, n3) = recvBuffTopRon->elem(i, j);
			tCon->elem(i + 1, j + 1, n3) = recvBuffTopT->elem(i, j);
			tnCon->elem(i + 1, j + 1, n3) = recvBuffTopTn->elem(i, j);
			u13->elem(i + 1, j + 1, n3 + 1) = recvBuffTopU13->elem(i, j);
			u23->elem(i + 1, j + 1, n3 + 1) = recvBuffTopU23->elem(i, j);
			u33->elem(i + 1, j + 1, n3 + 1) = recvBuffTopU33->elem(i, j);
			ro3->elem(i + 1, j + 1, n3 + 1) = recvBuffTopRo3->elem(i, j);
			p3->elem(i + 1, j + 1, n3 + 1) = recvBuffTopP3->elem(i, j);
		}
	}

	delete sendBuffBottomU1;
	delete sendBuffBottomU2;
	delete sendBuffBottomU3;
	delete sendBuffBottomRo;
	delete sendBuffBottomRon;
	delete sendBuffBottomT;
	delete sendBuffBottomTn;
	delete sendBuffBottomU13;
	delete sendBuffBottomU23;
	delete sendBuffBottomU33;
	delete sendBuffBottomRo3;
	delete sendBuffBottomP3;

	delete sendBuffTopU1;
	delete sendBuffTopU2;
	delete sendBuffTopU3;
	delete sendBuffTopRo;
	delete sendBuffTopRon;
	delete sendBuffTopT;
	delete sendBuffTopTn;
	delete sendBuffTopU13;
	delete sendBuffTopU23;
	delete sendBuffTopU33;
	delete sendBuffTopRo3;
	delete sendBuffTopP3;

	delete recvBuffBottomU1;
	delete recvBuffBottomU2;
	delete recvBuffBottomU3;
	delete recvBuffBottomRo;
	delete recvBuffBottomRon;
	delete recvBuffBottomT;
	delete recvBuffBottomTn;
	delete recvBuffBottomU13;
	delete recvBuffBottomU23;
	delete recvBuffBottomU33;
	delete recvBuffBottomRo3;
	delete recvBuffBottomP3;

	delete recvBuffTopU1;
	delete recvBuffTopU2;
	delete recvBuffTopU3;
	delete recvBuffTopRo;
	delete recvBuffTopRon;
	delete recvBuffTopT;
	delete recvBuffTopTn;
	delete recvBuffTopU13;
	delete recvBuffTopU23;
	delete recvBuffTopU33;
	delete recvBuffTopRo3;
	delete recvBuffTopP3;

	if (coordX == 0) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k < n3; ++k) {
				// on the east plane
				u1nCon->elem(0, j, k) = u1nCon->elem(1, j, k);
				u2nCon->elem(0, j, k) = u2nCon->elem(1, j, k);
				u3nCon->elem(0, j, k) = u3nCon->elem(1, j, k);
				ronCon->elem(0, j, k) = ronCon->elem(1, j, k);
				tnCon->elem(0, j, k) = tnCon->elem(1, j, k);
			}
		}
	}

	if (coordX == px - 1) {
		for (int j = 1; j < n2; ++j) {
			for (int k = 1; k < n3; ++k) {
				// on the west plane
				u1nCon->elem(n1, j, k) = u1nCon->elem(n1 - 1, j, k);
				u2nCon->elem(n1, j, k) = u2nCon->elem(n1 - 1, j, k);
				u3nCon->elem(n1, j, k) = u3nCon->elem(n1 - 1, j, k);
				ronCon->elem(n1, j, k) = ronCon->elem(n1 - 1, j, k);
				tnCon->elem(n1, j, k) = tnCon->elem(n1 - 1, j, k);
			}
		}
	}


	if (coordZ == 0)
	{	for (int i = 1; i < n1; ++i) {
			for (int j = 1; j < n2; ++j) {
				// on the top plane
				u1nCon->elem(i, j, 0) = u1nCon->elem(i, j, 1);
				u2nCon->elem(i, j, 0) = u2nCon->elem(i, j, 1);
				u3nCon->elem(i, j, 0) = u3nCon->elem(i, j, 1);
				ronCon->elem(i, j, 0) = ronCon->elem(i, j, 1);
				roCon->elem(i, j, 0) = roCon->elem(i, j, 1);
				tnCon->elem(i, j, 0) = tnCon->elem(i, j, 1);
				tCon->elem(i, j, 0) = tCon->elem(i, j, 1);
				u13->elem(i, j, 0) = u13->elem(i, j, 1);
				u23->elem(i, j, 0) = u23->elem(i, j, 1);
				u33->elem(i, j, 0) = u33->elem(i, j, 1);
				ro3->elem(i, j, 0) = ro3->elem(i, j, 1);
				p3->elem(i, j, 0) = p3->elem(i, j, 1);
			}
		}
	}

	if (coordZ == pz - 1) {
		for (int i = 1; i < n1; ++i) {
			for (int j = 1; j < n2; ++j) {
				// on the bottom plane
				u1nCon->elem(i, j, n3) = u1nCon->elem(i, j, n3 - 1);
				u2nCon->elem(i, j, n3) = u2nCon->elem(i, j, n3 - 1);
				u3nCon->elem(i, j, n3) = u3nCon->elem(i, j, n3 - 1);
				ronCon->elem(i, j, n3) = ronCon->elem(i, j, n3 - 1);
				roCon->elem(i, j, n3) = roCon->elem(i, j, n3 - 1);
				tnCon->elem(i, j, n3) = tnCon->elem(i, j, n3 - 1);
				tCon->elem(i, j, n3) = tCon->elem(i, j, n3 - 1);
				u13->elem(i, j, n3 + 1) = u13->elem(i, j, n3);
				u23->elem(i, j, n3 + 1) = u23->elem(i, j, n3);
				u33->elem(i, j, n3 + 1) = u33->elem(i, j, n3);
				ro3->elem(i, j, n3 + 1) = ro3->elem(i, j, n3);
				p3->elem(i, j, n3 + 1) = p3->elem(i, j, n3);
			}
		}
	}

	double  u1_f, u1_b, u1_cn, u1_c,
	        u2_f, u2_fn, u2_b, u2_bn, u2_cn, u2_c,
	        u3_f, u3_fn, u3_b, u3_bn, u3_cn, u3_c,
	        ro_cn, ro_c,
	        t_f, t_fn, t_b, t_bn, t_cn, t_c,
	        p_f, p_b, p_cn, p_c,
	        r_f, r_fn, r_b, r_cn, r_c,
	        q_f, q_b, q_bn, q_cn, q_c;

	double cond_f, cond_b;

	double gr, gt, gu2, gu3, gq;

	double rmax, rmin, qmax, qmin, tmax, tmin, u2_max, u2_min, u3_max, u3_min;

	double qn, pn, rn, ro_n, tn, un, u2_n, u3_n, ucf, ucb;

	// first local invariants for the interior faces puts in the buffer arrays, bypass on the center of the cell
	// then by taking into account the boundary condition calculates extreme elements of the buffers
	// and only then calculates the flow variables

	// flow variables calculation on DS1 faces orthogonal X1 axis

	// bypassing along the X1 axis

	// only interior faces !
	for (int k = 1; k < n3; k++)
	{
		for (int j = 1; j < n2; j++)
		{
			double _i = 0, _n1 = n1;
			if (coordX == 0) {
				_i++;
				rBuf[1] = (condition->elem(1, j, k) < 0.5) ? 1. : 0.;
				tfBuf[1] = t0;
				u2fBuf[1] = 0.;
				u3fBuf[1] = 0.;
			}

			if (coordX == px - 1) {
				_n1--;
				qBuf[n1] = 5. / 6.;
				tbBuf[n1] = t0;
				u2bBuf[n1] = 0.;
				u3bBuf[n1] = 0.;
			}

			for (int i = _i; i <= _n1; i++)
			{
				u1_f = u11->elem(i + 1, j, k);
				u1_b = u11->elem(i, j, k);
				u1_cn = u1nCon->elem(i, j, k);
				u1_c = u1Con->elem(i, j, k);

				u2_f = u21->elem(i + 1, j, k);
				u2_b = u21->elem(i, j, k);
				u2_cn = u2nCon->elem(i, j, k);
				u2_c = u2Con->elem(i, j, k);

				u3_f = u31->elem(i + 1, j, k);
				u3_b = u31->elem(i, j, k);
				u3_cn = u3nCon->elem(i, j, k);
				u3_c = u3Con->elem(i, j, k);

				ro_cn = ronCon->elem(i, j, k);
				ro_c = roCon->elem(i, j, k);

				t_f = t1->elem(i + 1, j, k);
				t_b = t1->elem(i, j, k);
				t_cn = tnCon->elem(i, j, k);
				t_c = tCon->elem(i, j, k);

				p_f = p1->elem(i + 1, j, k);
				p_b = p1->elem(i, j, k);
				p_cn = sound * sound * (ro_cn - ro0_g);
				p_c = sound * sound * (ro_c - ro0_g);

				// invariant calculation

				r_f = u1_f + p_f / (ro0_g * sound);
				r_b = u1_b + p_b / (ro0_g * sound);
				r_cn = u1_cn + p_cn / (ro0_g * sound);
				r_c = u1_c + p_c / (ro0_g * sound);

				r_fn = 2 * r_cn - r_b;

				q_f = u1_f - p_f / (ro0_g * sound);
				q_b = u1_b - p_b / (ro0_g * sound);
				q_cn = u1_cn - p_cn / (ro0_g * sound);
				q_c = u1_c - p_c / (ro0_g * sound);

				q_bn = 2 * q_cn - q_f;

				t_fn = 2 * t_cn - t_b;
				t_bn = 2 * t_cn - t_f;

				u2_fn = 2 * u2_cn - u2_b;
				u2_bn = 2 * u2_cn - u2_f;

				u3_fn = 2 * u3_cn - u3_b;
				u3_bn = 2 * u3_cn - u3_f;

				// the permissible range of changes
				gr = 2 * (r_cn - r_c) / dt + (u1_cn + sound) * (r_f - r_b) / dx1;
				gq = 2 * (r_cn - r_c) / dt + (u1_cn - sound) * (q_f - q_b) / dx1;

				gt = 2 * (t_cn - t_c) / dt + u1_cn * (t_f - t_b) / dx1;
				gu2 = 2 * (u2_cn - u2_c) / dt + u1_cn * (u2_f - u2_b) / dx1;
				gu3 = 2 * (u3_cn - u3_c) / dt + u1_cn * (u3_f - u3_b) / dx1;

				// RMAX=MAX(RF,RC,RB) +dt*GR
				rmax = max3d(r_f, r_c, r_b) + dt * gr;

				// RMIN=MIN(RF,RC,RB) +dt*GR
				rmin = min3d(r_f, r_c, r_b) + dt * gr;

				// QMAX=MAX(QF,QC,QB) +dt*GQ
				qmax = max3d(q_f, q_c, q_b) + dt * gq;

				// QMIN=MIN(QF,QC,QB) +dt*GQ
				qmin = min3d(q_f, q_c, q_b) + dt * gq;

				// TMAX=MAX(TF,TC,TB) +dt*GT
				tmax = max3d(t_f, t_c, t_b) + dt * gt;

				// TMIN=MIN(TF,TC,TB) +dt*GT
				tmin = min3d(t_f, t_c, t_b) + dt * gt;

				// U2MAX=MAX(U2F,U2C,U2B) +dt*GU2
				u2_max = max3d(u2_f, u2_c, u2_b) + dt * gu2;

				// U2MIN=MIN(U2F,U2C,U2B) +dt*GU2
				u2_min = min3d(u2_f, u2_c, u2_b) + dt * gu2;

				// U3MAX=MAX(U3F,U3C,U3B) +dt*GU3
				u3_max = max3d(u3_f, u3_c, u3_b) + dt * gu3;

				// U3MIN=MIN(U3F,U3C,U3B) +dt*GU3
				u3_min = min3d(u3_f, u3_c, u3_b) + dt * gu3;

				// invariants correction
				if (r_fn > rmax) r_fn = rmax;
				if (r_fn < rmin) r_fn = rmin;

				if (q_bn > qmax) q_bn = qmax;
				if (q_bn < qmin) q_bn = qmin;

				if (t_fn > tmax) t_fn = tmax;
				if (t_fn < tmin) t_fn = tmin;

				if (t_bn > tmax) t_bn = tmax;
				if (t_bn < tmin) t_bn = tmin;

				if (u2_fn > u2_max) u2_fn = u2_max;
				if (u2_fn < u2_min) u2_fn = u2_min;

				if (u2_bn > u2_max) u2_bn = u2_max;
				if (u2_bn < u2_min) u2_bn = u2_min;

				if (u3_fn > u3_max) u3_fn = u3_max;
				if (u3_fn < u3_min) u3_fn = u3_min;

				if (u3_bn > u3_max) u3_bn = u3_max;
				if (u3_bn < u3_min) u3_bn = u3_min;

				// put invariants to buffers
				rBuf[i + 1] = r_fn;
				qBuf[i] = q_bn;
				tfBuf[i + 1] = t_fn;
				tbBuf[i] = t_bn;
				u2fBuf[i + 1] = u2_fn;
				u2bBuf[i] = u2_bn;
				u3fBuf[i + 1] = u3_fn;
				u3bBuf[i] = u3_bn;
			}

			// boundary conditions along the X1 axis
			// assignment of boundary invatiants and add them to the buffer arrays

			// no-slip conditions



// 			if (u1nCon->elem(n1 - 1, j , k) >= 0) {
// 				un = 2 * u1nCon->elem(n1 - 1, j, k) - u11->elem(n1, j, k);
// 				u2_n = u2fBuf[n1];
// 				u3_n = u2fBuf[n1];

// 				rn = rBuf[n1];
// 				pn = (rn - un) * sound * ro0_g;
// 				ro_n = ro0_g + pn / (sound * sound);
// 				tn = tfBuf[n1];
// 			} else {
// 				rn = rBuf[n1];
// 				qn = 1.;

// 				un = (rn + qn) / 2;
// 				pn = (rn - qn) * sound * ro0_g / 2;
// 				ro_n = ro0_g + pn / (sound * sound);

// 				u2_n = 0.;
// 				u3_n = 0.;
// 				tn = t0;

// //                un = u1nCon->elem(n1 - 1, j, k);
// //                u2_n = u2nCon->elem(n1 - 1, j, k);
// //                u3_n = u3nCon->elem(n1 - 1, j, k);
// //                ro_n = ronCon->elem(n1 - 1, j, k);
// //                pn = (ro_n - ro0_g) * (sound*sound);
// //                tn = tnCon->elem(n1 - 1, j, k);
// 			}

// 			p1->elem(n1, j, k) = pn;
// 			u11->elem(n1, j, k) = un;
// 			ro1->elem(n1, j, k) = ro_n;
// 			t1->elem(n1, j, k) = tn;
// 			u21->elem(n1, j, k) = u2_n;
// 			u31->elem(n1, j, k) = u3_n;

			// the flow variables calculations
			for (int i = 1; i <= n1; i++) {
				cond_b = condition->elem(i - 1, j, k);
				cond_f = condition->elem(i, j, k);

				if (cond_b < 0.5 && cond_f < 0.5) {
					rn = rBuf[i];
					qn = qBuf[i];

					un = (rn + qn) / 2;
					pn = (rn - qn) * sound * ro0_g / 2;
					ro_n = ro0_g + pn / (sound * sound);

					ucf = u1nCon->elem(i, j, k);
					ucb = u1nCon->elem(i - 1, j, k);

					if (ucf >= 0 && ucb >= 0) {
						tn = tfBuf[i];
						u2_n = u2fBuf[i];
						u3_n = u3fBuf[i];
					}
					else if (ucf <= 0 && ucb <= 0) {
						tn = tbBuf[i];
						u2_n = u2bBuf[i];
						u3_n = u3bBuf[i];
					}
					else if (ucb >= 0 && ucf <= 0) {
						if (ucb > -ucf) {
							tn = tfBuf[i];
							u2_n = u2fBuf[i];
							u3_n = u3fBuf[i];
						}
						else {
							tn = tbBuf[i];
							u2_n = u2bBuf[i];
							u3_n = u3bBuf[i];
						}
					}
					else {
						tn = tnCon->elem(i, j, k) + tnCon->elem(i - 1, j, k) - t1->elem(i, j, k);
						u2_n = u2nCon->elem(i, j, k) + u2nCon->elem(i - 1, j, k) - u21->elem(i, j, k);
						u3_n = u3nCon->elem(i, j, k) + u3nCon->elem(i - 1, j, k) - u31->elem(i, j, k);
					}
				} else {
					un = 0.;
					if (cond_b > 0.5 && cond_f < 0.5) {
						if (u1nCon->elem(i, j, k) <= 0) {
							qn = qBuf[i];
							u2_n = /*0.;*/u2bBuf[i];
							u3_n = /*0.;*/u3bBuf[i];
							pn = -qn * sound * ro0_g;
							ro_n = ro0_g + pn / (sound * sound);
							tn = tbBuf[i];
						} else {
							u2_n = /*0.;*/u2nCon->elem(i, j, k);
							u3_n = /*0.;*/u3nCon->elem(i, j, k);
							ro_n = ronCon->elem(i, j, k);
							pn = (ro_n - ro0_g) * (sound * sound);
							tn = tnCon->elem(i, j, k);
						}
					} else if (cond_b < 0.5 && cond_f > 0.5) {
						if (u1nCon->elem(i - 1, j, k) >= 0) {
							rn = rBuf[i];
							u2_n = /*0.;*/u2fBuf[i];
							u3_n = /*0.;*/u3fBuf[i];
							pn = rn * sound * ro0_g;
							ro_n = ro0_g + pn / (sound * sound);
							tn = tfBuf[i];
						} else {
							u2_n = /*0.;*/u2nCon->elem(i - 1, j, k);
							u3_n = /*0.;*/u3nCon->elem(i - 1, j, k);
							ro_n = ronCon->elem(i - 1, j, k);
							pn = (ro_n - ro0_g) * (sound * sound);
							tn = tnCon->elem(i - 1, j, k);
						}
					} else {
						u2_n = 0.;
						u3_n = 0.;
						pn = 0.;
						ro_n = ro0_g;
						tn = t0;
					}
				}

				p1->elem(i, j, k) = pn;
				u11->elem(i, j, k) = un;
				ro1->elem(i, j, k) = ro_n;
				t1->elem(i, j, k) = tn;
				u21->elem(i, j, k) = u2_n;
				u31->elem(i, j, k) = u3_n;
			}
		}
	}

	// flow variables calculation on DS2 faces orthogonal X2 axis

	// bypassing along the X2 axis

	double u1_fn, u1_bn;

	double gu1;

	double u1_max, u1_min, u1_n;

	for (int k = 1; k < n3; k++)
	{
		for (int i = 1; i < n1; i++)
		{
			for (int j = 0; j <= n2; j++)
			{
				u2_f = u22->elem(i, j + 1, k);
				u2_b = u22->elem(i, j, k);
				u2_cn = u2nCon->elem(i, j, k);
				u2_c = u2Con->elem(i, j, k);

				u1_f = u12->elem(i, j + 1, k);
				u1_b = u12->elem(i, j, k);
				u1_cn = u1nCon->elem(i, j, k);
				u1_c = u1Con->elem(i, j, k);

				u3_f = u32->elem(i, j + 1, k);
				u3_b = u32->elem(i, j, k);
				u3_cn = u3nCon->elem(i, j, k);
				u3_c = u3Con->elem(i, j, k);

				ro_cn = ronCon->elem(i, j, k);
				ro_c = roCon->elem(i, j, k);

				t_f = t2->elem(i, j + 1, k);
				t_b = t2->elem(i, j, k);
				t_cn = tnCon->elem(i, j, k);
				t_c = tCon->elem(i, j, k);

				p_f = p2->elem(i, j + 1, k);
				p_b = p2->elem(i, j, k);
				p_cn = sound * sound * (ro_cn - ro0_g);
				p_c = sound * sound * (ro_c - ro0_g);

				// invariant calculation
				r_f = u2_f + p_f / (ro0_g * sound);
				r_b = u2_b + p_b / (ro0_g * sound);
				r_cn = u2_cn + p_cn / (ro0_g * sound);
				r_c = u2_c + p_c / (ro0_g * sound);

				r_fn = 2 * r_cn - r_b;

				q_f = u2_f - p_f / (ro0_g * sound);
				q_b = u2_b - p_b / (ro0_g * sound);
				q_cn = u2_cn - p_cn / (ro0_g * sound);
				q_c = u2_c - p_c / (ro0_g * sound);

				q_bn = 2 * q_cn - q_f;

				t_fn = 2 * t_cn - t_b;
				t_bn = 2 * t_cn - t_f;

				u1_fn = 2 * u1_cn - u1_b;
				u1_bn = 2 * u1_cn - u1_f;

				u3_fn = 2 * u3_cn - u3_b;
				u3_bn = 2 * u3_cn - u3_f;

				// the permissible range of changes
				gr = 2 * (r_cn - r_c) / dt + (u2_cn + sound) * (r_f - r_b) / dx2;
				gq = 2 * (q_cn - q_c) / dt + (u2_cn - sound) * (q_f - q_b) / dx2;
				gt = 2 * (t_cn - t_c) / dt + u2_cn * (t_f - t_b) / dx2;
				gu1 = 2 * (u1_cn - u1_c) / dt + u2_cn * (u1_f - u1_b) / dx2;
				gu3 = 2 * (u3_cn - u3_c) / dt + u2_cn * (u3_f - u3_b) / dx2;


				// RMAX=MAX(RF,RC,RB) +dt*GR
				rmax = max3d(r_f, r_c, r_b) + dt * gr;

				// RMIN=MIN(RF,RC,RB) +dt*GR
				rmin = min3d(r_f, r_c, r_b) + dt * gr;

				// QMAX=MAX(QF,QC,QB) +dt*GQ
				qmax = max3d(q_f, q_c, q_b) + dt * gq;

				// QMIN=MIN(QF,QC,QB) +dt*GQ
				qmin = min3d(q_f, q_c, q_b) + dt * gq;

				// TMAX=MAX(TF,TC,TB) +dt*GT
				tmax = max3d(t_f, t_c, t_b) + dt * gt;

				// TMIN=MIN(TF,TC,TB) +dt*GT
				tmin = min3d(t_f, t_c, t_b) + dt * gt;

				// U1MAX=MAX(U1F,U1C,U1B) +dt*GU1
				u1_max = max3d(u1_f, u1_c, u1_b) + dt * gu1;

				// U1MIN=MIN(U1F,U1C,U1B) +dt*GU1
				u1_min = min3d(u1_f, u1_c, u1_b) + dt * gu1;

				// U3MAX=MAX(U3F,U3C,U3B) +dt*GU3
				u3_max = max3d(u3_f, u3_c, u3_b) + dt * gu3;

				// U3MIN=MIN(U3F,U3C,U3B) +dt*GU3
				u3_min = min3d(u3_f, u3_c, u3_b) + dt * gu3;

				// invariants correction
				if (r_fn > rmax) r_fn = rmax;
				if (r_fn < rmin) r_fn = rmin;

				if (q_bn > qmax) q_bn = qmax;
				if (q_bn < qmin) q_bn = qmin;

				if (t_fn > tmax) t_fn = tmax;
				if (t_fn < tmin) t_fn = tmin;

				if (t_bn > tmax) t_bn = tmax;
				if (t_bn < tmin) t_bn = tmin;

				if (u1_fn > u1_max) u1_fn = u1_max;
				if (u1_fn < u1_min) u1_fn = u1_min;

				if (u1_bn > u1_max) u1_bn = u1_max;
				if (u1_bn < u1_min) u1_bn = u1_min;

				if (u3_fn > u3_max) u3_fn = u3_max;
				if (u3_fn < u3_min) u3_fn = u3_min;

				if (u3_bn > u3_max) u3_bn = u3_max;
				if (u3_bn < u3_min) u3_bn = u3_min;

				// put invariants to buffers
				// ==================================================
				// !!! IMPORTANT !!!
				// ==================================================
				// u2fBuf and u2bBuf are actially the u1fBuf and u1bBuf
				// It's not an error. We do it to save dynamic memory
				// ==================================================
				rBuf[j + 1] = r_fn;
				qBuf[j] = q_bn;
				tfBuf[j + 1] = t_fn;
				tbBuf[j] = t_bn;
				u2fBuf[j + 1] = u1_fn;
				u2bBuf[j] = u1_bn;
				u3fBuf[j + 1] = u3_fn;
				u3bBuf[j] = u3_bn;
			}

			// the flow variables calculations
			for (int j = 1; j <= n2; j++)
			{
				cond_b = condition->elem(i, j - 1, k);
				cond_f = condition->elem(i, j, k);

				if (cond_b < 0.5 && cond_f < 0.5) {
					rn = rBuf[j];
					qn = qBuf[j];

					un = (rn + qn) / 2;

					pn = (rn - qn) * sound * ro0_g / 2;
					ro_n = ro0_g + pn / (sound * sound);

					ucf = u2nCon->elem(i, j, k);
					ucb = u2nCon->elem(i, j - 1, k);

					if (ucf >= 0 && ucb >= 0) {
						tn = tfBuf[j];
						u1_n = u2fBuf[j];
						u3_n = u3fBuf[j];
					}
					else if (ucf <= 0 && ucb <= 0) {
						tn = tbBuf[j];
						u1_n = u2bBuf[j];
						u3_n = u3bBuf[j];
					}
					else if (ucb >= 0 && ucf <= 0) {
						if (ucb > -ucf) {
							tn = tfBuf[j];
							u1_n = u2fBuf[j];
							u3_n = u3fBuf[j];
						}
						else {
							tn = tbBuf[j];
							u1_n = u2bBuf[j];
							u3_n = u3bBuf[j];
						}
					}
					else {
						tn = tnCon->elem(i, j, k) + tnCon->elem(i, j - 1, k) - t2->elem(i, j, k);
						u1_n = u1nCon->elem(i, j, k) + u1nCon->elem(i, j - 1, k) - u12->elem(i, j, k);
						u3_n = u3nCon->elem(i, j, k) + u3nCon->elem(i, j - 1, k) - u32->elem(i, j, k);
					}
				} else {
					un = 0.;
					if (cond_b > 0.5 && cond_f < 0.5) {
						if (u2nCon->elem(i, j, k) <= 0) {
							qn = qBuf[j];
							pn = -qn * sound * ro0_g;
							ro_n = ro0_g + pn / (sound * sound);
							u1_n = /*0.;*/u2bBuf[j];
							u3_n = /*0.;*/u3bBuf[j];
							tn = tbBuf[i];
						} else {
							ro_n = ronCon->elem(i, j, k);
							pn = (ro_n - ro0_g) * (sound * sound);
							u1_n = /*0.;*/u2nCon->elem(i, j, k);
							u3_n = /*0.;*/u3nCon->elem(i, j, k);
							tn = tnCon->elem(i, j, k);
						}
					} else if (cond_b < 0.5 && cond_f > 0.5) {
						if (u2nCon->elem(i, j - 1, k) >= 0) {
							rn = rBuf[j];
							pn = rn * sound * ro0_g;
							ro_n = ro0_g + pn / (sound * sound);
							u1_n = /*0.;*/u2fBuf[j];
							u3_n = /*0.;*/u3fBuf[j];
							tn = tfBuf[j];
						} else {
							ro_n = ronCon->elem(i, j - 1, k);
							pn = (ro_n - ro0_g) * (sound * sound);
							u1_n = /*0.;*/u2nCon->elem(i, j - 1, k);
							u3_n = /*0.;*/u3nCon->elem(i, j - 1, k);
							tn = tnCon->elem(i, j - 1, k);
						}
					} else {
						pn = 0.;
						ro_n = ro0_g;
						u1_n = 0.;
						u3_n = 0.;
						tn = t0;
					}
				}

				p2->elem(i, j, k) = pn;
				u22->elem(i, j, k) = un;
				ro2->elem(i, j, k) = ro_n;
				t2->elem(i, j, k) = tn;
				u12->elem(i, j, k) = u1_n;
				u32->elem(i, j, k) = u3_n;
			}
		}
	}

	// flow variables calculation on DS3 faces orthogonal X3 axis

	// bypassing along the X3 axis

	for (int i = 1; i < n1; i++)
	{
		for (int j = 1; j < n2; j++)
		{
			for (int k = 0; k <= n3; k++)
			{
				u3_f = u33->elem(i, j, k + 1);
				u3_b = u33->elem(i, j, k);
				u3_cn = u3nCon->elem(i, j, k);
				u3_c = u3Con->elem(i, j, k);

				u1_f = u13->elem(i, j, k + 1);
				u1_b = u13->elem(i, j, k);
				u1_cn = u1nCon->elem(i, j, k);
				u1_c = u1Con->elem(i, j, k);

				u2_f = u23->elem(i, j, k + 1);
				u2_b = u23->elem(i, j, k);
				u2_cn = u2nCon->elem(i, j, k);
				u2_c = u2Con->elem(i, j, k);

				ro_cn = ronCon->elem(i, j, k);
				ro_c = roCon->elem(i, j, k);

				t_f = t3->elem(i, j, k + 1);
				t_b = t3->elem(i, j, k);
				t_cn = tnCon->elem(i, j, k);
				t_c = tCon->elem(i, j, k);

				p_f = p3->elem(i, j, k + 1);
				p_b = p3->elem(i, j, k);
				p_cn = sound * sound * (ro_cn - ro0_g);
				p_c = sound * sound * (ro_c - ro0_g);

				// invariant calculation
				r_f = u3_f + p_f / (ro0_g * sound);
				r_b = u3_b + p_b / (ro0_g * sound);
				r_cn = u3_cn + p_cn / (ro0_g * sound);
				r_c = u3_c + p_c / (ro0_g * sound);

				r_fn = 2 * r_cn - r_b;

				q_f = u3_f - p_f / (ro0_g * sound);
				q_b = u3_b - p_b / (ro0_g * sound);
				q_cn = u3_cn - p_cn / (ro0_g * sound);
				q_c = u3_c - p_c / (ro0_g * sound);

				q_bn = 2 * q_cn - q_f;

				t_fn = 2 * t_cn - t_b;
				t_bn = 2 * t_cn - t_f;

				u2_fn = 2 * u2_cn - u2_b;
				u2_bn = 2 * u2_cn - u2_f;

				u1_fn = 2 * u1_cn - u1_b;
				u1_bn = 2 * u1_cn - u1_f;

				// the permissible range of changes
				gr = 2 * (r_cn - r_c) / dt + (u3_cn + sound) * (r_f - r_b) / dx3;
				gq = 2 * (r_cn - r_c) / dt + (u3_cn - sound) * (q_f - q_b) / dx3;

				gt = 2 * (t_cn - t_c) / dt + u3_cn * (t_f - t_b) / dx3;
				gu1 = 2 * (u1_cn - u1_c) / dt + u3_cn * (u1_f - u1_b) / dx3;
				gu2 = 2 * (u2_cn - u2_c) / dt + u3_cn * (u2_f - u2_b) / dx3;

				// RMAX=MAX(RF,RC,RB) +dt*GR
				rmax = max3d(r_f, r_c, r_b) + dt * gr;

				// RMIN=MIN(RF,RC,RB) +dt*GR
				rmin = min3d(r_f, r_c, r_b) + dt * gr;

				// QMAX=MAX(QF,QC,QB) +dt*GQ
				qmax = max3d(q_f, q_c, q_b) + dt * gq;

				// QMIN=MIN(QF,QC,QB) +dt*GQ
				qmin = min3d(q_f, q_c, q_b) + dt * gq;

				// TMAX=MAX(TF,TC,TB) +dt*GT
				tmax = max3d(t_f, t_c, t_b) + dt * gt;

				// TMIN=MIN(TF,TC,TB) +dt*GT
				tmin = min3d(t_f, t_c, t_b) + dt * gt;

				// U1MAX=MAX(U1F,U1C,U1B) +dt*GU1
				u1_max = max3d(u1_f, u1_c, u1_b) + dt * gu1;

				// U1MIN=MIN(U1F,U1C,U1B) +dt*GU1
				u1_min = min3d(u1_f, u1_c, u1_b) + dt * gu1;

				// U2MAX=MAX(U2F,U2C,U2B) +dt*GU2
				u2_max = max3d(u2_f, u2_c, u2_b) + dt * gu2;

				// U2MIN=MIN(U2F,U2C,U2B) +dt*GU2
				u2_min = min3d(u2_f, u2_c, u2_b) + dt * gu2;

				// invariants correction
				if (r_fn > rmax) r_fn = rmax;
				if (r_fn < rmin) r_fn = rmin;

				if (q_bn > qmax) q_bn = qmax;
				if (q_bn < qmin) q_bn = qmin;

				if (t_fn > tmax) t_fn = tmax;
				if (t_fn < tmin) t_fn = tmin;

				if (t_bn > tmax) t_bn = tmax;
				if (t_bn < tmin) t_bn = tmin;

				if (u1_fn > u1_max) u1_fn = u1_max;
				if (u1_fn < u1_min) u1_fn = u1_min;

				if (u1_bn > u1_max) u1_bn = u1_max;
				if (u1_bn < u1_min) u1_bn = u1_min;

				if (u2_fn > u2_max) u2_fn = u2_max;
				if (u2_fn < u2_min) u2_fn = u2_min;

				if (u2_bn > u2_max) u2_bn = u2_max;
				if (u2_bn < u2_min) u2_bn = u2_min;

				// put invariants to buffers
				// ====================================================
				// !!! IMPORTANT !!!
				// ====================================================
				// u2fBuf and u2bBuf are actially the u1fBuf and u1bBuf
				// u3fBuf and u3bBuf are actially the u2fBuf and u2bBuf
				// It's not an error. We do it to save dynamic memory
				// ====================================================
				rBuf[k + 1] = r_fn;
				qBuf[k] = q_bn;
				tfBuf[k + 1] = t_fn;
				tbBuf[k] = t_bn;
				u2fBuf[k + 1] = u1_fn;
				u2bBuf[k] = u1_bn;
				u3fBuf[k + 1] = u2_fn;
				u3bBuf[k] = u2_bn;
			}

			// the flow variables calculations
			for (int k = 1; k <= n3; k++)
			{
				cond_b = condition->elem(i, j, k - 1);
				cond_f = condition->elem(i, j, k);

				if (cond_b < 0.5 && cond_f < 0.5) {
					rn = rBuf[k];
					qn = qBuf[k];

					pn = (rn - qn) * sound * ro0_g / 2;
					un = (rn + qn) / 2;

					ro_n = ro0_g + pn / (sound * sound);

					ucf = u3nCon->elem(i, j, k);
					ucb = u3nCon->elem(i, j, k - 1);

					if (ucf >= 0 && ucb >= 0) {
						tn = tfBuf[k];
						u1_n = u2fBuf[k];
						u2_n = u3fBuf[k];
					}
					else if (ucf <= 0 && ucb <= 0) {
						tn = tbBuf[k];
						u1_n = u2bBuf[k];
						u2_n = u3bBuf[k];
					}
					else if (ucb >= 0 && ucf <= 0) {
						if (ucb > -ucf) {
							tn = tfBuf[k];
							u1_n = u2fBuf[k];
							u2_n = u3fBuf[k];
						}
						else {
							tn = tbBuf[k];
							u1_n = u2bBuf[k];
							u2_n = u3bBuf[k];
						}
					}
					else {
						tn = tnCon->elem(i, j, k) + tnCon->elem(i, j, k - 1) - t3->elem(i, j, k);
						u1_n = u1nCon->elem(i, j, k) + u1nCon->elem(i, j, k - 1) - u13->elem(i, j, k);
						u2_n = u2nCon->elem(i, j, k) + u2nCon->elem(i, j, k - 1) - u23->elem(i, j, k);
					}
				} else {
					un = 0.;
					if (cond_b > 0.5 && cond_f < 0.5) {
						if (u3nCon->elem(i, j, k) <= 0) {
							qn = qBuf[k];
							pn = -qn * sound * ro0_g;
							ro_n = ro0_g + pn / (sound * sound);
							u1_n = /*0.;*/u2bBuf[k];
							u2_n = /*0.;*/u3bBuf[k];
							tn = tbBuf[k];
						} else {
							ro_n = ronCon->elem(i, j, k);
							pn = (ro_n - ro0_g) * (sound * sound);
							u1_n = /*0.;*/u1nCon->elem(i, j, k);
							u2_n = /*0.;*/u2nCon->elem(i, j, k);
							tn = tnCon->elem(i, j, k);
						}
					} else if (cond_b < 0.5 && cond_f > 0.5) {
						if (u3nCon->elem(i, j, k - 1) >= 0) {
							rn = rBuf[k];
							pn = rn * sound * ro0_g;
							ro_n = ro0_g + pn / (sound * sound);
							u1_n = /*0.;*/u2fBuf[k];
							u2_n = /*0.;*/u3fBuf[k];
							tn = tfBuf[k];
						} else {
							ro_n = ronCon->elem(i, j, k - 1);
							pn = (ro_n - ro0_g) * (sound * sound);
							u1_n = /*0.;*/u1nCon->elem(i, j, k - 1);
							u2_n = /*0.;*/u2nCon->elem(i, j, k - 1);
							tn = tnCon->elem(i, j, k - 1);
						}
					} else {
						pn = 0.;
						ro_n = ro0_g;
						u1_n = 0.;
						u2_n = 0.;
						tn = t0;
					}
				}

				p3->elem(i, j, k) = pn;
				u33->elem(i, j, k) = un;
				ro3->elem(i, j, k) = ro_n;
				t3->elem(i, j, k) = tn;
				u13->elem(i, j, k) = u1_n;
				u23->elem(i, j, k) = u2_n;
			}
		}
	}
}

// void WriteDataTecplot() {
// 	char filename[100];
// 	sprintf(filename, "%saverage.vtk", dirPath);

// 	char header[200], v[15];
// 	sprintf(header, "TITLE=\"OUT\"\nVARIABLES = \"X\" \"Y\" \"U1\" \"U2\"\n
// 	        ZONE I=%d, J=%d, DATAPACKING=BLOCK, VARLOCATION=([3-4]=CELLCENTERED)\n", (n1_g - 1) , (n2_g - 1));

// 	MPI_Offset offset = 0, localOffset;

// 	MPI_File fwh;
// 	MPI_File_open(newComm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fwh);

// 	if (rank == 0) {
// 		MPI_Offset off = 0;

// 		MPI_File_write_at(fwh, off, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
// 		off += strlen(header) * sizeof(char);

// 		// X1
// 		for (int j = 1; j <= n2_g; ++j) {
// 			for (int i = 1; i <= n1_g; ++i) {
// 				sprintf(v, "%15E", x1[i]);
// 				MPI_File_write_at(fwh, off, v, 15, MPI_CHAR, MPI_STATUS_IGNORE);
// 				off += 15;
// 			}
// 			fprintf(fd, "\n");
// 			++off;
// 		}
// 		fprintf(fd, "\n");
// 		++off;


// 		// X2
// 		for (int j = 1; j <= n2_g; ++j) {
// 			for (int i = 1; i <= n1_g; ++i) {
// 				sprintf(v, "%15E", x2[i]);
// 				MPI_File_write_at(fwh, off, v, 15, MPI_CHAR, MPI_STATUS_IGNORE);
// 				off += 15;
// 			}
// 			fprintf(fd, "\n");
// 			++off;
// 		}
// 		fprintf(fd, "\n");
// 		++off;

// 		offset = off;
// 	}

// 	MPI_Bcast(&offset, 1, MPI_Offset, 0, newComm);

// 	// U1
// 	for (int j = 1; j < n2; ++j) {
// 		for (int i = 1; i < n1; ++i) {
// 			sprintf(v, "%15E", u1nCon->elem(i, j, k));
// 			MPI_File_write_at(fwh, off, v, 15, MPI_CHAR, MPI_STATUS_IGNORE);
// 			off += 15;
// 		}
// 		fprintf(fd, "\n");
// 		++off;
// 	}
// 	fprintf(fd, "\n");
// 	++off;

// 	// U2
// 	for (int j = 1; j < n2; ++j) {
// 		for (int i = 1; i < n1; ++i) {
// 			sprintf(v, "%15E", u2nCon->elem(i, j, k));
// 			MPI_File_write_at(fwh, off, v, 15, MPI_CHAR, MPI_STATUS_IGNORE);
// 			off += 15;
// 		}
// 		fprintf(fd, "\n");
// 		++off;
// 	}
// 	fprintf(fd, "\n");
// 	++off;


// }

void WriteDataParaView() {
	char filename[100];
	double v;

	sprintf(filename, "%sout_%d.vtk", dirPath, nStep);

	char header1[50], header2[200], headerX[100], headerY[100], headerZ[100], headerU[200],
	     headerPC[50], headerR[50];

	sprintf(header1, "# vtk DataFile Version 3.0\nvtk output\nBINARY\n");
	sprintf(header2, "DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d", n1_g, n2_g, n3_g);
	sprintf(headerX, "\nX_COORDINATES %d double\n", n1_g);
	sprintf(headerY, "\nY_COORDINATES %d double\n", n2_g);
	sprintf(headerZ, "\nZ_COORDINATES %d double\n", n3_g);
	sprintf(headerU, "\nCELL_DATA %d\nVECTORS U double\n", (n1_g - 1) * (n2_g - 1) * (n3_g - 1));
	sprintf(headerPC, "\nscalars PC double\nLOOKUP_TABLE default\n");
	sprintf(headerR, "\nVECTORS R double\n");


	MPI_Offset offset =
	    (strlen(header1) + strlen(header2) + strlen(headerX) + strlen(headerY) + strlen(headerZ) + strlen(headerU)) * sizeof(char) +
	    (n1_g + n2_g + n3_g) * sizeof(double), localOffset;

	MPI_File fwh;

	MPI_File_open(newComm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fwh);

	if (rank == 0) {
		MPI_Offset off = 0;

		int nmax = (int) max3d(n1_g, n2_g, n3_g);
		double *xbuff = new double[nmax];

		MPI_File_write_at(fwh, off, header1, strlen(header1), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(header1) * sizeof(char);

		MPI_File_write_at(fwh, off, header2, strlen(header2), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(header2) * sizeof(char);

		MPI_File_write_at(fwh, off, headerX, strlen(headerX), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(headerX) * sizeof(char);
		for (int i = 1; i <= n1_g; ++i) {
			v = (i - 1) * dx1;
			swap8(&v);
			xbuff[i - 1] = v;
		}
		MPI_File_write_at(fwh, off, xbuff, n1_g, MPI_DOUBLE, MPI_STATUS_IGNORE);
		off += n1_g * sizeof(double);

		MPI_File_write_at(fwh, off, headerY, strlen(headerY), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(headerY) * sizeof(char);
		for (int j = 1; j <= n2_g; ++j) {
			v = (j - 1) * dx2;
			swap8(&v);
			xbuff[j - 1] = v;
		}
		MPI_File_write_at(fwh, off, xbuff, n2_g, MPI_DOUBLE, MPI_STATUS_IGNORE);
		off += n2_g * sizeof(double);

		MPI_File_write_at(fwh, off, headerZ, strlen(headerZ), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(headerZ) * sizeof(char);
		for (int k = 1; k <= n3_g; ++k) {
			v = (k - 1) * dx3;
			swap8(&v);
			xbuff[k - 1] = v;
		}
		MPI_File_write_at(fwh, off, xbuff, n3_g, MPI_DOUBLE, MPI_STATUS_IGNORE);
		off += n3_g * sizeof(double);

		delete[] xbuff;

		MPI_File_write_at(fwh, off, headerU, strlen(headerU), MPI_CHAR, MPI_STATUS_IGNORE);
	}

	// U
	double *buff = new double[(n1 - 1) * 3];

	for (int k = 1; k < n3; ++k)
	{
		for (int j = 1; j < n2; ++j)
		{
			for (int i = 1; i < n1; ++i)
			{
				v = u1nCon->elem(i, j, k);
				swap8(&v);
				buff[(i - 1) * 3] = v;

				v = u2nCon->elem(i, j, k);
				swap8(&v);
				buff[(i - 1) * 3 + 1] = v;

				v = u3nCon->elem(i, j, k);
				swap8(&v);
				buff[(i - 1) * 3 + 2] = v;
			}
			localOffset =
			    offset + 3 * ((n1_g - 1) * ((n2 - 1) * coordY + j - 1) + coordX * (n1 - 1)) * sizeof(double);
			MPI_File_write_at(fwh, localOffset, buff, (n1 - 1) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}
	// ------------------------------------------------------------------------------------------------------------

	offset += 3 * ((n1_g - 1) * (n2_g - 1) * (n3_g - 1)) * sizeof(double);

	// PC
	if (rank == 0) MPI_File_write_at(fwh, offset, headerPC, strlen(headerPC), MPI_CHAR, MPI_STATUS_IGNORE);

	offset += strlen(headerPC) * sizeof(char);

	for (int k = 1; k < n3; ++k)
	{
		for (int j = 1; j < n2; ++j)
		{
			for (int i = 1; i < n1; ++i)
			{
				v = ronCon->elem(i, j, k);
				swap8(&v);
				buff[i - 1] = v;
			}
			localOffset =
			    offset + ((n1_g - 1) * ((n2 - 1) * coordY + j - 1) + coordX * (n1 - 1)) * sizeof(double);
			MPI_File_write_at(fwh, localOffset, buff, n1 - 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}
	// ------------------------------------------------------------------------------------------------------------

	offset += ((n1_g - 1) * (n2_g - 1) * (n3_g - 1)) * sizeof(double);

	// R
	if (rank == 0) MPI_File_write_at(fwh, offset, headerR, strlen(headerR), MPI_CHAR, MPI_STATUS_IGNORE);

	offset += strlen(headerR) * sizeof(char);

	for (int k = 1; k < n3; k++)
	{
		for (int j = 1; j < n2; j++)
		{
			for (int i = 1; i < n1; i++)
			{
				double d2u3 = (u32->elem(i, j + 1, k) - u32->elem(i, j, k)) / dx2,
				       d3u2 = (u23->elem(i, j, k + 1) - u23->elem(i, j, k)) / dx3,
				       d3u1 = (u13->elem(i, j, k + 1) - u13->elem(i, j, k)) / dx3,
				       d1u3 = (u31->elem(i + 1, j, k) - u31->elem(i, j, k)) / dx1,
				       d1u2 = (u21->elem(i + 1, j, k) - u21->elem(i, j, k)) / dx1,
				       d2u1 = (u12->elem(i, j + 1, k) - u12->elem(i, j, k)) / dx2;

				v = d2u3 - d3u2;
				swap8(&v);
				buff[(i - 1) * 3] = v;

				v = d3u1 - d1u3;
				swap8(&v);
				buff[(i - 1) * 3 + 1] = v;

				v = d1u2 - d2u1;
				swap8(&v);
				buff[(i - 1) * 3 + 2] = v;
			}
			localOffset =
			    offset + 3 * ((n1_g - 1) * ((n2 - 1) * coordY + j - 1) + coordX * (n1 - 1)) * sizeof(double);
			MPI_File_write_at(fwh, localOffset, buff, (n1 - 1) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}

	MPI_File_close(&fwh);
	// ------------------------------------------------------------------------------------------------------------

	delete buff;
}

void WriteAverageValues(double time) {
	char filename[100];
	double v;

	FILE *fd;
	sprintf(filename, "%saverage.vtk", dirPath);

	char header1[50], header2[200], headerX[100], headerY[100], headerZ[100], headerU[200];

	sprintf(header1, "# vtk DataFile Version 3.0\nvtk output\nBINARY\n");
	sprintf(header2, "DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d", n1_g, n2_g, n3_g);
	sprintf(headerX, "\nX_COORDINATES %d double\n", n1_g);
	sprintf(headerY, "\nY_COORDINATES %d double\n", n2_g);
	sprintf(headerZ, "\nZ_COORDINATES %d double\n", n3_g);
	sprintf(headerU, "\nCELL_DATA %d\nVECTORS U double\n", (n1_g - 1) * (n2_g - 1) * (n3_g - 1));


	MPI_Offset offset =
	    (strlen(header1) + strlen(header2) + strlen(headerX) + strlen(headerY) + strlen(headerZ) + strlen(headerU)) * sizeof(char) +
	    (n1_g + n2_g + n3_g) * sizeof(double), localOffset;

	MPI_File fwh;

	MPI_File_open(newComm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fwh);

	if (rank == 0) {
		MPI_Offset off = 0;

		int nmax = (int) max3d(n1_g, n2_g, n3_g);
		double *xbuff = new double[nmax];

		MPI_File_write_at(fwh, off, header1, strlen(header1), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(header1) * sizeof(char);

		MPI_File_write_at(fwh, off, header2, strlen(header2), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(header2) * sizeof(char);

		MPI_File_write_at(fwh, off, headerX, strlen(headerX), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(headerX) * sizeof(char);
		for (int i = 1; i <= n1_g; ++i) {
			v = (i - 1) * dx1;
			swap8(&v);
			xbuff[i - 1] = v;
		}
		MPI_File_write_at(fwh, off, xbuff, n1_g, MPI_DOUBLE, MPI_STATUS_IGNORE);
		off += n1_g * sizeof(double);

		MPI_File_write_at(fwh, off, headerY, strlen(headerY), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(headerY) * sizeof(char);
		for (int j = 1; j <= n2_g; ++j) {
			v = (j - 1) * dx2;
			swap8(&v);
			xbuff[j - 1] = v;
		}
		MPI_File_write_at(fwh, off, xbuff, n2_g, MPI_DOUBLE, MPI_STATUS_IGNORE);
		off += n2_g * sizeof(double);

		MPI_File_write_at(fwh, off, headerZ, strlen(headerZ), MPI_CHAR, MPI_STATUS_IGNORE);
		off += strlen(headerZ) * sizeof(char);
		for (int k = 1; k <= n3_g; ++k) {
			v = (k - 1) * dx3;
			swap8(&v);
			xbuff[k - 1] = v;
		}
		MPI_File_write_at(fwh, off, xbuff, n3_g, MPI_DOUBLE, MPI_STATUS_IGNORE);
		off += n3_g * sizeof(double);

		delete[] xbuff;

		MPI_File_write_at(fwh, off, headerU, strlen(headerU), MPI_CHAR, MPI_STATUS_IGNORE);
	}

	// U
	double *buff = new double[(n1 - 1) * 3];

	for (int k = 1; k < n3; ++k)
	{
		for (int j = 1; j < n2; ++j)
		{
			for (int i = 1; i < n1; ++i)
			{
				v = averageU1->elem(i, j, k) / (TIME - time);
				swap8(&v);
				buff[(i - 1) * 3] = v;

				v = averageU2->elem(i, j, k) / (TIME - time);
				swap8(&v);
				buff[(i - 1) * 3 + 1] = v;

				v = averageU3->elem(i, j, k) / (TIME - time);
				swap8(&v);
				buff[(i - 1) * 3 + 2] = v;
			}
			localOffset =
			    offset + 3 * ((n1_g - 1) * ((n2 - 1) * coordY + j - 1) + coordX * (n1 - 1)) * sizeof(double);
			MPI_File_write_at(fwh, localOffset, buff, (n1 - 1) * 3, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}
	MPI_File_close(&fwh);
	// ------------------------------------------------------------------------------------------------------------

	delete buff;
}

void WriteEnergy() {
	double energy = 0;
	for (int i = 1; i < n1; i++) {
		for (int j = 1; j < n2; j++) {
			for (int k = 1; k < n3; k++) {
				energy += (u1nCon->elem(i, j, k) * u1nCon->elem(i, j, k) +
				           u2nCon->elem(i, j, k) * u2nCon->elem(i, j, k) +
				           u3nCon->elem(i, j, k) * u3nCon->elem(i, j, k)) / 2;
			}
		}
	}

	double allEnergy = 0.;
	MPI_Reduce(&energy, &allEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, newComm);

	if (rank == 0) {
		char filename[100];
		sprintf(filename, "%senergy.txt", dirPath);

		double volume = (x3_t - x3_b) * (x2_n - x2_s) * (x1_e - x1_w);
		allEnergy *= volume;
		FILE *fd = fopen(filename, "a");
		fprintf(fd, "(%f; %f) ", TIME, allEnergy);
		fclose(fd);
	}
}

void AverageValues() {
	for (int k = 1; k < n3; k++)
	{
		for (int j = 1; j < n2; j++)
		{
			for (int i = 1; i < n1; i++)
			{
				averageU1->elem(i, j, k) += u1nCon->elem(i, j, k) * dt;
				averageU2->elem(i, j, k) += u2nCon->elem(i, j, k) * dt;
				averageU3->elem(i, j, k) += u3nCon->elem(i, j, k) * dt;
			}
		}
	}
}

void FreeMemory() {
	delete[] x1;
	delete[] x2;
	delete[] x3;

	delete roCon;
	delete u1Con;
	delete u2Con;
	delete u3Con;
	delete tCon;

	delete ronCon;
	delete u1nCon;
	delete u2nCon;
	delete u3nCon;
	delete tnCon;

	delete condition;

	delete averageU1;
	delete averageU2;
	delete averageU3;

	delete ro1;
	delete t1;
	delete u11;
	delete u21;
	delete u31;
	delete p1;

	delete ro2;
	delete t2;
	delete u12;
	delete u22;
	delete u32;
	delete p2;

	delete ro3;
	delete t3;
	delete u13;
	delete u23;
	delete u33;
	delete p3;

	delete f1;
	delete f2;
	delete f3;

	delete sigm11;
	delete sigm21;
	delete sigm31;

	delete sigm12;
	delete sigm22;
	delete sigm32;

	delete sigm13;
	delete sigm23;
	delete sigm33;

	delete[] rBuf;
	delete[] qBuf;
	delete[] tfBuf;
	delete[] tbBuf;
	delete[] u2fBuf;
	delete[] u2bBuf;
	delete[] u3fBuf;
	delete[] u3bBuf;
}

double min3d(double x1, double x2, double x3) {
	x1 = x1 > x2 ? x2 : x1;
	x1 = x1 > x3 ? x3 : x1;
	return x1;
}

double max3d(double x1, double x2, double x3) {
	x1 = x1 < x2 ? x2 : x1;
	x1 = x1 < x3 ? x3 : x1;
	return x1;
}

static void swap8(void *v) {
	if (needSwap) {
		char    in[8], out[8];
		memcpy(in, v, 8);
		out[0] = in[7];
		out[1] = in[6];
		out[2] = in[5];
		out[3] = in[4];
		out[4] = in[3];
		out[5] = in[2];
		out[6] = in[1];
		out[7] = in[0];
		memcpy(v, out, 8);
	}
}