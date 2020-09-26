#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#define	X	105	// number of voxels in x direction
#define	Y	115	// number of voxels in y direction
#define Z	115	// number of voxels in z direction

#define	dx	0.1	// voxel size in x direction (mm)
#define	dy	0.1	// voxel size in y direction (mm)
#define dz	0.1	// voxel size in z direction (mm)

#define xa	50	// apical x coordinate of centroid
#define ya	72	// apical y coordinate of centroid
#define za	9 // apical z coordinate of centroid

#define xb	55	// basal x coordinate of centroid
#define yb	65	// basal y coordinate of centroid
#define zb	104	// basal z coordinate of centroid

#define geo_threshold   0.1     // threshold for geometry (values below this in the geometry input file will be regarded as outside the tissue)

int main () {

	FILE *geofile;
	FILE *ev1xfile;
	FILE *ev1yfile;
	FILE *ev1zfile;
	FILE *evtkfile;
	FILE *hvtkfile;

	double geo;					// geometry flag
	double e1x;					// x component of the primary eigenvector
	double e1y;					// y component of the primary eigenvector
	double e1z;					// z component of the primary eigenvector

	int centroid_x[X+Y+Z];		// list of centroid locations
	int centroid_y[X+Y+Z];		// list of centroid locations
	int centroid_z[X+Y+Z];		// list of centroid locations
	int centroid_cnt = 0;		// centroid location counter
	int cnt;					// centroid location counter
	int centroid_flag;			// centroid location flag
	int x, y, z;				// x, y and z locations of current voxel
	int xc, yc, zc;				// x, y and z location of centroid relating to the current voxel
	double dradius;				// distance from voxel to centroid

	double long_x, long_y, long_z;			// components of longitudinal vector
	double radial_x, radial_y, radial_z;	// components of radial vector
	double tangent_x, tangent_y, tangent_z;	// components of tangential vector

	double e1s_x, e1s_y, e1s_z;	// components of the primary eigenvector projected onto short axis plane
	double e1t_x, e1t_y, e1t_z; // components of the primary eigenvector projected onto tangential plane

	double angle_e1t_l;			// angle between the projection of the primary eigenvector onto the tangential plane and the longitudinal vector
	double angle_e1t_t;			// angle between the projection of the primary eigenvector onto the tangential plane and the tangential vector
	double angle_e1s_t;			// angle between the projection of the primary eigenvector onto the short axis plane and the tangential vector
	double angle_e1s_r;			// angle between the projection of the primary eigenvector onto the short axis plane and the radial vector

	double e_length;			// vector length, for normalisation

	double helix;				// angles
	double trans;				// angles

	double pi;					// pi

	geofile = fopen("geo_hybrid_het.txt", "rt");
	ev1xfile = fopen("e1x_MH1_crop.txt", "rt");
	ev1yfile = fopen("e1y_MH1_crop.txt", "rt");
	ev1zfile = fopen("e1z_MH1_crop.txt", "rt");
	evtkfile = fopen("eigenvec_MH1.vtk", "w");
	hvtkfile = fopen("helix_MH1.vtk", "w");

	// print VTK file headers
	printf("printing headers...\n");

	fprintf(evtkfile, "# vtk DataFile Version 3.0\n");
	fprintf(evtkfile, "vtk output\n");
	fprintf(evtkfile, "ASCII\n");
	fprintf(evtkfile, "DATASET STRUCTURED_POINTS\n");
	fprintf(evtkfile, "DIMENSIONS %d %d %d\n", X, Y, Z);
	fprintf(evtkfile, "SPACING %f %f %f\n", 1.0*dx, 1.0*dy, 1.0*dz);
	fprintf(evtkfile, "ORIGIN 0 0 0\n");
	fprintf(evtkfile, "POINT_DATA %d\n", X*Y*Z);
	fprintf(evtkfile, "SCALARS ImageFile float 3\n");
	fprintf(evtkfile, "LOOKUP_TABLE default\n");

	fprintf(hvtkfile, "# vtk DataFile Version 3.0\n");
	fprintf(hvtkfile, "vtk output\n");
	fprintf(hvtkfile, "ASCII\n");
	fprintf(hvtkfile, "DATASET STRUCTURED_POINTS\n");
	fprintf(hvtkfile, "DIMENSIONS %d %d %d\n", X, Y, Z);
	fprintf(hvtkfile, "SPACING %f %f %f\n", 1.0*dx, 1.0*dy, 1.0*dz);
	fprintf(hvtkfile, "ORIGIN 0 0 0\n");
	fprintf(hvtkfile, "POINT_DATA %d\n", X*Y*Z);
	fprintf(hvtkfile, "SCALARS ImageFile float 1\n");
	fprintf(hvtkfile, "LOOKUP_TABLE default\n");

	printf("printing headers complete...\n");

	pi = 3.141592654;

	printf("calculating central axis...\n");

	// Calculate longitudinal unit vector (i.e. pointing along the long axis from apex to base)
	e_length = sqrt( (double)((double)xb-(double)xa)*dx*((double)xb-(double)xa)*dx + ((double)yb-(double)ya)*dy*((double)yb-(double)ya)*dy + ((double)zb-(double)za)*dz*((double)zb-(double)za)*dz );
	long_x = ((double)xb-(double)xa)*dx / e_length;
	long_y = ((double)yb-(double)ya)*dy / e_length;
	long_z = ((double)zb-(double)za)*dz / e_length;

	if ( sqrt(long_x*long_x + long_y*long_y + long_z*long_z) < 0.99 || sqrt(long_x*long_x + long_y*long_y + long_z*long_z) > 1.01 ) {printf("unit vector out of range\n"); exit(0);}

	// calculate long axis (centroid) locations
	for (x=1; x<=X; x++) {
		yc = yb-(((xb-x)*(yb-ya))/(double)(xb-xa));
		zc = zb-(((xb-x)*(zb-za))/(double)(xb-xa));
		if (yc>=0 && yc<=Y && zc>=0 && zc<=Z) {
	        	centroid_cnt++;
	        	centroid_x[centroid_cnt] = x;
	        	centroid_y[centroid_cnt] = yc;
	        	centroid_z[centroid_cnt] = zc;	
		}
	}
	for (y=1; y<=Y; y++) {
		xc = xb-(((yb-y)*(xb-xa))/(double)(yb-ya));
		zc = zb-(((yb-y)*(zb-za))/(double)(yb-ya));
		if (xc>=0 && xc<=X && zc>=0 && zc<=Z) {
	        	centroid_cnt++;
	        	centroid_x[centroid_cnt] = xc;
	        	centroid_y[centroid_cnt] = y;
	        	centroid_z[centroid_cnt] = zc;	
		}
	}
	for (z=1; z<=Z; z++) {
		xc = xb-(((zb-z)*(xb-xa))/(double)(zb-za));
		yc = yb-(((zb-z)*(yb-ya))/(double)(zb-za));
		if (xc>=0 && xc<=X && yc>=0 && yc<=Y) {
	        	centroid_cnt++;
	        	centroid_x[centroid_cnt] = xc;
	        	centroid_y[centroid_cnt] = yc;
	        	centroid_z[centroid_cnt] = z;	
		}
	}

	printf("calculating central axis complete...\n");

	float (*vec)[Y][Z][3] = malloc(sizeof(*vec)*X);

	printf("calculating angles...\n");

	// calculate angles
	for (z=1; z<=Z; z++) {
	for (y=1; y<=Y; y++) {
	for (x=1; x<=X; x++) {

		// read in eigenvector components and geometry
		fscanf(geofile, "%lf ", &geo);
		fscanf(ev1xfile, "%lf ", &e1x);
		fscanf(ev1yfile, "%lf ", &e1y);
		fscanf(ev1zfile, "%lf ", &e1z);

		centroid_flag = 0;				// check whether this voxel is a centroid
		for (cnt=1; cnt<=centroid_cnt; cnt++)
			if (x==centroid_x[cnt] && y==centroid_y[cnt] && z==centroid_z[cnt]) centroid_flag = 1;
		if (centroid_flag == 1) {			// this voxel is a centroid, so print 1000 to output files and move on
			fprintf(evtkfile, "%f %f %f ", 0.0,0.0,0.0);
			fprintf(hvtkfile, "1000.0 ");
		}

		else if (geo<geo_threshold) {			// this voxel is a centroid, so print -500 to output files and move on
			fprintf(evtkfile, "%f %f %f ", 0.0,0.0,0.0);
			fprintf(hvtkfile, "-500.0 ");
		}

		else {						// this voxel is tissue, so calculate angles and print to output files

			// CALCULATE LOCAL COORDINATE SYSTEM

	        // ensure eigenvectors are normalised
	        e_length = sqrt(e1x*e1x + e1y*e1y + e1z*e1z);
	        e1x /= e_length;
	        e1y /= e_length;
	        e1z /= e_length;

			// find this voxel's centroid (i.e. the closest point to the voxel on the long axis)
			dradius = 1000000.0;
			for (cnt=1; cnt<=centroid_cnt; cnt++) {
				e_length = sqrt( (double)(x-centroid_x[cnt])*dx*(x-centroid_x[cnt])*dx + (y-centroid_y[cnt])*dy*(y-centroid_y[cnt])*dy + (z-centroid_z[cnt])*dz*(z-centroid_z[cnt])*dz );
				if (e_length<dradius) {
					dradius = e_length;
					xc = centroid_x[cnt];
					yc = centroid_y[cnt];
					zc = centroid_z[cnt];
				}
			}

	        // calculate radial unit vector (i.e. pointing from centroid to voxel)
	        e_length = sqrt( ((double)x-(double)xc)*dx*((double)x-(double)xc)*dx + ((double)y-(double)yc)*dy*((double)y-(double)yc)*dy + ((double)z-(double)zc)*dz*((double)z-(double)zc)*dz );
	        radial_x = ((double)x-(double)xc)*dx / e_length;
	        radial_y = ((double)y-(double)yc)*dy / e_length;
	        radial_z = ((double)z-(double)zc)*dz / e_length;

	        // calculate tangential unit vector (i.e. the cross product of the longitudinal and radial unit vectors).
	        // The tangential vector lies in the short axis plane and, from the right hand rule, points anti-clockwise when looking from base to apex.
	        tangent_x = long_y * radial_z - long_z * radial_y;
	        tangent_y = long_z * radial_x - long_x * radial_z;
	        tangent_z = long_x * radial_y - long_y * radial_x;

			// CALCULATE FIBRE (HELIX AND TRANSVERSE) ANGLES

	        // DETERMINE IF PRIMARY (FIBRE) EIGENVECTOR NEEDS TO BE FLIPPED
	        // calculate projection of primary eigenvector onto the short axis plane then normalise
	        e1s_x = e1x - (e1x*long_x + e1y*long_y + e1z*long_z) * long_x;
	        e1s_y = e1y - (e1x*long_x + e1y*long_y + e1z*long_z) * long_y;
	        e1s_z = e1z - (e1x*long_x + e1y*long_y + e1z*long_z) * long_z;
	        e_length = sqrt(e1s_x*e1s_x + e1s_y*e1s_y + e1s_z*e1s_z);

	        e1s_x /= e_length;
	        e1s_y /= e_length;
	        e1s_z /= e_length;
	        // calculate angle between the projection of primary eigenvector onto short axis plane and the tangential vector
	        angle_e1s_t = acos(e1s_x*tangent_x + e1s_y*tangent_y + e1s_z*tangent_z);
	        // flip the primary eigenvector if necessary
			if (angle_e1s_t > (pi/2)) {
				e1x = -e1x;
				e1y = -e1y;
				e1z = -e1z;
				// re-calculate projections of this new primary eigenvector onto the short axis plane then normalise
				e1s_x = -e1s_x;
				e1s_y = -e1s_y;
				e1s_z = -e1s_z;
				// re-calculate angle between the projection of primary eigenvector onto short axis plane and the tangential vector
				angle_e1s_t = pi - angle_e1s_t;
			}

			// calculate projection of primary eigenvector onto tangential plane then normalise
	        e1t_x = e1x - (e1x*radial_x + e1y*radial_y + e1z*radial_z) * radial_x;
	        e1t_y = e1y - (e1x*radial_x + e1y*radial_y + e1z*radial_z) * radial_y;
	        e1t_z = e1z - (e1x*radial_x + e1y*radial_y + e1z*radial_z) * radial_z;
	        e_length = sqrt(e1t_x*e1t_x + e1t_y*e1t_y + e1t_z*e1t_z);
	        e1t_x /= e_length;
	        e1t_y /= e_length;
	        e1t_z /= e_length;

	        // calculate angle between the projection of primary eigenvector onto tangential plane and the longitudinal vector
	        angle_e1t_l = acos(e1t_x*long_x + e1t_y*long_y + e1t_z*long_z);

	        // calculate angle between the projection of primary eigenvector onto tangential plane and the tangential vector
	        angle_e1t_t = acos(e1t_x*tangent_x + e1t_y*tangent_y + e1t_z*tangent_z);

	        // calculate angle between the projection of primary eigenvector onto short axis plane and the radial vector
	        angle_e1s_r = acos(e1s_x*radial_x + e1s_y*radial_y + e1s_z*radial_z);

			// calculate fibre helix (elevation) angle
			if (angle_e1t_l > (pi/2)) helix = -angle_e1t_t;
			else helix = angle_e1t_t;

			// calculate fibre transverse (rotation) angle
			if (angle_e1s_r > (pi/2)) trans = -angle_e1s_t;
			else trans = angle_e1s_t;

			// convert radians to degrees
			helix *= 57.29577951;
			trans *= 57.29577951;

			vec[x][y][z][0] = e1x;
			vec[x][y][z][1] = e1y;
			vec[x][y][z][2] = e1z;

			// print angles to output files
			fprintf(evtkfile, "%.1f %.1f %.1f ", vec[x][y][z][0], vec[x][y][z][1], vec[x][y][z][2]);
			fprintf(hvtkfile, "%.1f ", helix);

		} // end of calculating and printing angles

	}
	fscanf(geofile, "\n");
	fscanf(ev1xfile, "\n");
	fscanf(ev1yfile, "\n");
	fscanf(ev1zfile, "\n");
	fprintf(evtkfile, "\n");
	fprintf(hvtkfile, "\n");
	}
	printf("slice %d of %d finished...\n", z, Z);
	}

	printf("FINISHED\n");

	fclose (geofile);
	fclose (ev1xfile);
	fclose (ev1yfile);
	fclose (ev1zfile);
	fclose (evtkfile);
	fclose (hvtkfile);

	free(vec);

} // end of main
